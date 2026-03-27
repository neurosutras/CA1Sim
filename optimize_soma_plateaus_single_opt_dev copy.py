import os
import sys
import yaml
import numpy as np
import random
import time
import datetime
import copy
import matplotlib.pyplot as plt
from scipy.optimize import minimize, basinhopping
from scipy.signal import find_peaks
from neuron import h, gui
from specify_cells import CA1_Pyr, QuickSim
import csv

# ==================================================================================================
# 1. SETUP & TARGETS
# ==================================================================================================

target_plateau_per_cycle = [3.337, 3.5297, 4.356, 5.2235, 5.6729]
target_spikes_per_cycle  = [6, 6, 7, 7, 7]
target_trough_per_cycle  = [1.2946, 3.1996, 8.0606, 17.4871, 22.5451]

# FIX: Centralized burst timing constants — were duplicated across calculate_cycle_features()
# and objective(), which would silently diverge if the protocol ever changed.
BURST_STARTS_MS   = [500.0, 650.0, 800.0, 950.0, 1100.0]
BURST_DURATION_MS = 150.0
TROUGH_WINDOW_MS  = 50.0
STIM_WIN_START_MS = 500.0
STIM_WIN_END_MS   = 1250.0
 

# Optimization Configuration
#config_path = 'optimization_config_BasinHopping_20260322_4th_Run_Overnight.yaml'
config_path = 'optimization_config_BasinHopping_20260324_Run_Overnight.yaml'
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)

# Override mech_filename to use the best-so-far biophysics as the starting point
#config['simulation']['mech_filename'] = 'recovered_biophysics_Best_Model_so_Far_20260322_2311.yaml'
config['simulation']['mech_filename'] = 'recovered_biophysics_eval_43_20260324_2113.yaml'

# Simulation parameters
duration  = config['simulation']['equilibrate'] + config['simulation']['sim_duration']
dt        = config['simulation']['dt']
v_init    = config['simulation']['v_init']
tbs_times = config['simulation']['tbs_times']

# Search budgets
# MAX_EVALS_BH = 40    # 40 jumps but kill switch fires at 40 total evals
# MAX_ITERS_BH = 3     # Limit L-BFGS-B to 3 iters per jump — keeps it moving
# MAX_EVALS_NM = 10     # ~30 mins polish


#IF YOU CHANGE THE NUMBER OF PARAMETERS for MAX_EVALS_BH, YOU MUST CHANGE:
#minimizer_kwargs = {
    # "method": "Nelder-Mead",
    # "options": {
    #     "maxiter": MAX_ITERS_BH,
    #     "maxfev":  MAX_ITERS_BH * 6,  <------------
# MAX_JUMPS_BH  = 8                              # number of basin-hopping jumps # 8 × 3 × 6 = 144 evals → ~8.5 hrs
# MAX_ITERS_BH  = 3                                # Nelder-Mead iters per jump
# MAX_EVALS_BH  = MAX_JUMPS_BH * MAX_ITERS_BH * 6 # safety ceiling: 6 = num params (simplex size)
# MAX_EVALS_NM  = 15                               # Nelder-Mead polish evals

MAX_JUMPS_BH  = 25                               # 25 × 3 × 6 = 450 evals → ~8 hrs
MAX_ITERS_BH  = 3
MAX_EVALS_BH  = MAX_JUMPS_BH * MAX_ITERS_BH * 6
MAX_EVALS_NM  = 30                               # ~30 mins polish at the end                        # Nelder-Mead polish evals

# Optimizer seed — randomized each run for diverse starting points.
# NOTE: This is NOT the SYNAPSE_SEED (which stays fixed at 0 for reproducible synapse placement).
# This controls: (a) the random starting parameters, (b) numpy's RNG for basinhopping jumps.
OPT_SEED = int.from_bytes(os.urandom(4), 'big') % 100000  # 5-digit seed for readability
np.random.seed(OPT_SEED)
print(f"\n--- Optimizer Seed: {OPT_SEED} (use this to reproduce this run) ---")

# Global history for plotting and tracking
error_history    = []
evaluation_count = 0
best_error       = float('inf')
best_trace       = None
best_x           = None  # FIX: track best params explicitly for kill-switch recovery

# ==================================================================================================
# FIX: Custom kill-switch exception.
# Previously used StopIteration, which Python 3.7+ converts to RuntimeError inside generators,
# and scipy may silently swallow it — leaving the optimizer in an undefined state.
# A custom exception gives clean, predictable termination and allows state recovery.
# ==================================================================================================
class EvalLimitReached(Exception):
    pass

# ==================================================================================================
# CSV LOG INITIALIZATION (With Unique Timestamp)
# ==================================================================================================
start_time_str = datetime.datetime.now().strftime("%Y%m%d_%H%M")
log_filename = f'optimization_log_{start_time_str}.csv'

headers = [
    'Seed',
    'Eval', 'Total_Error',
    'cal', 'calH', 'car', 'cat', 'mykca', 'kca',
    'Area_C1', 'Area_C2', 'Area_C3', 'Area_C4', 'Area_C5',
    'Spks_C1', 'Spks_C2', 'Spks_C3', 'Spks_C4', 'Spks_C5',
    'Trough_C1', 'Trough_C2', 'Trough_C3', 'Trough_C4', 'Trough_C5'
]

with open(log_filename, 'w', newline='') as f:
    csv.writer(f).writerow(headers)
 
# FIX: Persistent line-buffered CSV writer — previously opened and closed on every evaluation,
# adding I/O overhead and risking corruption if the process was killed mid-write.
# buffering=1 gives line-buffered writes: safe against crashes, no open/close cost per eval.
_log_fh     = open(log_filename, 'a', newline='', buffering=1)
_log_writer = csv.writer(_log_fh)
 
print(f"--- Logging progress to: {log_filename} ---")

# ==================================================================================================
# 2. UTILITY FUNCTIONS
# ==================================================================================================

def calculate_cycle_features(trace, dt):
    """
    Returns areas (list), spikes (list), and troughs (list) for all 5 theta cycles.
    - Areas:   AUC (mV*s) within each 150ms burst window
    - Spikes:  AP count within each 150ms burst window
    - Troughs: minimum mV in the 50ms window just before the next burst onset
 
    FIX: Now returns three plain lists instead of (list, dict-of-lists, dict-of-lists).
    The dict-of-single-element-lists was unnecessarily awkward and forced [0] unwrapping
    at every call site.
 
    FIX: Converts trace to ndarray once at the top, then slices into views — previously
    called np.array(..., copy=True) on every cycle, creating multiple redundant copies.
    """
    # Single conversion; zero-copy if trace is already an ndarray
    t = np.asarray(trace) if trace is not None else None
 
    areas   = []
    spikes  = []
    troughs = []
 
    for cycle_num, start in enumerate(BURST_STARTS_MS, start=1):
        s_idx = int(start / dt)
        e_idx = int((start + BURST_DURATION_MS) / dt)
 
        # --- AREA ---
        if t is None or len(t) <= s_idx:
            areas.append(0.0)
        else:
            segment      = t[s_idx:min(len(t), e_idx)]          # view, no copy
            segment_rect = np.clip(segment, 0, None)
            areas.append(np.trapz(segment_rect, dx=(dt / 1000.0)))
 
        # --- SPIKES ---
        if t is None or len(t) <= s_idx:
            spikes.append(0)
        else:
            segment  = t[s_idx:min(len(t), e_idx)]              # view, no copy
            AP_peaks = find_peaks(segment, height=-50, width=(10, 1000), distance=50, prominence=20)[0]
            spikes.append(len(AP_peaks))
 
        # --- TROUGH ---
        if cycle_num < 5:
            next_burst = BURST_STARTS_MS[cycle_num]
            t_end      = int(next_burst / dt)
            t_start    = int((next_burst - TROUGH_WINDOW_MS) / dt)
        else:
            t_start = int(1100.0 / dt)
            t_end   = int(1250.0 / dt)
 
        if t is None or len(t) <= t_start:
            troughs.append(0.0)
        else:
            trough_seg = t[t_start:min(len(t), t_end)]          # view, no copy
            troughs.append(float(np.min(trough_seg)))
 
    return areas, spikes, troughs


# ==================================================================================================
# 3. CELL INITIALIZATION
# ==================================================================================================
print("--- Initializing Cell and Synaptic Infrastructure ---")
cell = CA1_Pyr(morph_filename=config['simulation']['morph_filename'],
               mech_filename=config['simulation']['mech_filename'],
               full_spines=False)

# --- FORCE GRADIENT ACROSS THE DENDRITE ---
# This forces the soma to be the most depolarized region, mimicking the experimental data.
cell.set_terminal_branch_na_gradient()

exc_syn_locs_by_sec_type = cell.get_excitatory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
inh_syn_locs_by_sec_type = cell.get_inhibitory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
exc_syn_types = ['AMPA_KIN', 'NMDA_KIN5']
inh_syn_types = ['GABAb']
SYNAPSE_SEED   = 0          # fixed seed — must match notebook Cell 11
local_random   = random.Random(SYNAPSE_SEED)

tbs_vec = h.Vector(tbs_times)

sim = QuickSim(duration, cvode=False, dt=dt, verbose=0)
sim.append_rec(cell, cell.tree.root, description='soma_v')
spike_vec = h.Vector()
cell.spike_detector.record(spike_vec)

stim_exc_syns = {'CA3': [], 'ECIII': []}
stim_inh_syns = {'CA3': [], 'ECIII': []}

# ==================================================================================================
# 4. SYNAPSE ASSIGNMENT LOGIC
# ==================================================================================================

def assign_exc_and_inh_synapse_stims(num_exc_syns, excitatory_stochastic, num_inh_syns):
    # Dictionaries to store the actual synapse objects
    stim_exc_syns = {'CA3': [], 'ECIII': []}
    stim_inh_syns = {'CA3': [], 'ECIII': []}

    # 1. Place Excitatory Synapses
    for pathway, num_to_insert in num_exc_syns.items():
        list_of_valid_syn_locs = []
        if pathway == 'ECIII':
            list_of_valid_syn_locs.extend(exc_syn_locs_by_sec_type['tuft'])
        else:
            for sec_type in ['trunk', 'apical']:
                list_of_valid_syn_locs.extend(exc_syn_locs_by_sec_type[sec_type])

        actual_num = min(int(num_to_insert), len(list_of_valid_syn_locs))
        selected_locs = local_random.sample(list_of_valid_syn_locs, actual_num)
        
        syn_list = cell.insert_synapses_at_syn_locs(selected_locs, exc_syn_types, stochastic=excitatory_stochastic)
        stim_exc_syns[pathway].extend(syn_list)
    
    # 2. Place Inhibitory Synapses
    for pathway, num_to_insert in num_inh_syns.items():
        list_of_valid_syn_locs = []
        if pathway == 'ECIII':
            list_of_valid_syn_locs.extend(inh_syn_locs_by_sec_type['tuft'])
        else:
            for sec_type in ['trunk', 'apical']:
                list_of_valid_syn_locs.extend(inh_syn_locs_by_sec_type[sec_type])

        actual_num = min(int(num_to_insert), len(list_of_valid_syn_locs))
        selected_locs = local_random.sample(list_of_valid_syn_locs, actual_num)
        
        syn_list = cell.insert_synapses_at_syn_locs(selected_locs, inh_syn_types, stochastic=False)
        stim_inh_syns[pathway].extend(syn_list)

    return stim_exc_syns, stim_inh_syns 


#The optimizer reuses the same cell object across every evaluation — it only changes conductances via modify_mech_param(). 
#It does NOT rebuild the cell from scratch each eval. So Without clearing: the synapse numbers will accumulate across evals.
# you clear synapses because the cell persists but synapses must be re-randomized identically each eval (via the fixed SYNAPSE_SEED), 
# and adding on top of existing ones would corrupt the simulation.
def clear_all_synapses():
    for node in cell.tree:
        node.content['synapses'] = []
    for pathway in stim_exc_syns: stim_exc_syns[pathway] = []
    for pathway in stim_inh_syns: stim_inh_syns[pathway] = []

# ==================================================================================================
# 5. PARAMETER UPDATE LOGIC
# ==================================================================================================

def update_calcium_and_potassium(m_cal, m_calH, m_car, m_cat, m_mykca, m_kca):
    cell.modify_mech_param('soma', 'cal', 'gcalbar', value=0.007 * m_cal) #optimize soma cal
    for sec in ['trunk', 'basal']:
        cell.modify_mech_param(sec, 'cal', 'gcalbar', origin='soma')
    for sec in ['apical', 'tuft']:
        cell.modify_mech_param(sec, 'cal', 'gcalbar', origin='trunk')

    cell.modify_mech_param('trunk', 'calH', 'gcalbar', value=0.000031635 * m_calH, max_loc=50.0, origin='soma') #optimize soma calH
    cell.modify_mech_param('trunk', 'calH', 'gcalbar', value=0.001455 * m_calH, min_loc=50.0, origin='soma', replace=False) #optimize trunk calH
    for sec in ['apical', 'tuft']:
        cell.modify_mech_param(sec, 'calH', 'gcalbar', origin='trunk')

    cell.modify_mech_param('soma', 'car', 'gcabar', value=0.003 * m_car) #optimize soma car
    cell.modify_mech_param('trunk', 'car', 'gcabar', value=0.00003 * m_car) #optimize trunk car
    cell.modify_mech_param('basal', 'car', 'gcabar', origin='soma')
    for sec in ['apical', 'tuft']:
        cell.modify_mech_param(sec, 'car', 'gcabar', origin='trunk')

    #The max value of the trunk is the free param - then compute the slope from max value
    cell.modify_mech_param('soma', 'cat', 'gcatbar', value=0.00005 * m_cat)
    cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0 * m_cat, max_loc=100.0, origin='soma') 
    cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0, slope=1.143e-06 * m_cat, min_loc=100.0, max_loc=350.0, origin='soma', replace=False) 
    cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0004 * m_cat, min_loc=350.0, origin='soma', replace=False) 
    cell.modify_mech_param('basal', 'cat', 'gcatbar', origin='soma')
    for sec in ['apical', 'tuft']:
        cell.modify_mech_param(sec, 'cat', 'gcatbar', origin='trunk')

    #The max value of the trunk is the free param - then compute the slope from max value
    cell.modify_mech_param('soma', 'mykca', 'gkbar', value=0.09075 * m_mykca) #optimize soma mykca
    cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.09075 * m_mykca, max_loc=50.0, origin='soma')
    cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.09075 * m_mykca, slope=-0.0005543 * m_mykca, min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
    cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.007585 * m_mykca, min_loc=200.0, origin='soma', replace=False)
    cell.modify_mech_param('tuft', 'mykca', 'gkbar', value=0.007585 * m_mykca)
    cell.modify_mech_param('basal', 'mykca', 'gkbar', origin='soma')
    cell.modify_mech_param('apical', 'mykca', 'gkbar', origin='trunk')

    #The max value of the trunk is the free param - then compute the slope from max value
    cell.modify_mech_param('soma', 'kca', 'gbar', value=0.0005 * m_kca)
    cell.modify_mech_param('trunk', 'kca', 'gbar', value=0.0005 * m_kca, max_loc=50.0, origin='soma')
    cell.modify_mech_param('trunk', 'kca', 'gbar', value=0.0005 * m_kca, slope=-3.056e-6 * m_kca, min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
    cell.modify_mech_param('trunk', 'kca', 'gbar', value=4.167e-05 * m_kca, min_loc=200.0, origin='soma', replace=False)
    cell.modify_mech_param('tuft', 'kca', 'gbar', value=4.167e-05 * m_kca)
    cell.modify_mech_param('basal', 'kca', 'gbar', origin='soma')
    cell.modify_mech_param('apical', 'kca', 'gbar', origin='trunk')

    all_secs  = ['soma', 'trunk', 'basal', 'apical', 'tuft']
    all_mechs = ['cal', 'calH', 'car', 'cat', 'mykca', 'kca', 'pas']
    for sec in all_secs:
        for mech in all_mechs:
            cell.reinitialize_subset_mechanisms(sec, mech)

#update model will assign the parameters to the cell, clear synapses, and reassign synapses each time it is called - for every evaluation
def update_model(x_vec):
    # FIX: Declare globals before assignment — previously update_model() created local
    # variables that shadowed the module-level stim_*_syns, leaving the globals empty
    # after clear_all_synapses() wiped them. The sim worked by accident (cell.tree is the
    # real source of truth) but the tracking state was silently wrong.
    global stim_exc_syns, stim_inh_syns
 
    m_cal, m_calH, m_car, m_cat, m_mykca, m_kca = x_vec
 
    try:
        n_exc_CA3   = config['parameters']['exc_CA3']
        n_exc_ECIII = config['parameters']['exc_ECIII']
        n_inh_CA3   = config['parameters']['inh_CA3']
        n_inh_ECIII = config['parameters']['inh_ECIII']
    except KeyError as e:
        print(f"YAML ERROR: Could not find key {e}. Check your 'parameters' block.")
        sys.exit()
 
    update_calcium_and_potassium(m_cal, m_calH, m_car, m_cat, m_mykca, m_kca)
 
    clear_all_synapses()
    # CRITICAL: reset the RNG to the fixed seed before every synapse assignment so
    # that synapse locations are IDENTICAL on every evaluation, independent of how
    # many prior evaluations have been run or whether the script was restarted.
    local_random.seed(SYNAPSE_SEED)
    num_exc = {'CA3': int(n_exc_CA3), 'ECIII': int(n_exc_ECIII)}
    num_inh = {'CA3': int(n_inh_CA3), 'ECIII': int(n_inh_ECIII)}
    stim_exc_syns, stim_inh_syns = assign_exc_and_inh_synapse_stims(num_exc, True, num_inh)
 
    cell.init_synaptic_mechanisms()
 
    for stim_dict in (stim_exc_syns, stim_inh_syns):
        for pathway in stim_dict:
            for syn in stim_dict[pathway]:
                syn.source.play(tbs_vec)

# ==================================================================================================
# 6. OBJECTIVE FUNCTION
# ==================================================================================================

def objective(x, kill_limit=None):
    global best_error, best_trace, best_x, evaluation_count
    evaluation_count += 1

    W_SPKS         = 1000.0
    W_AREA         =  800.0
    W_TROUGH       = 1200.0
    TROUGH_ASYM_HARD  = 3.0   # penalty multiplier when too hyperpolarized
    ROUGH_ASYM_SOFT  = 1.0   # penalty multiplier when too depolarized
 
    # FIX: Use custom EvalLimitReached instead of StopIteration.
    # In Python 3.7+ StopIteration raised inside a generator is silently converted to
    # RuntimeError, and scipy's internals may catch and swallow it entirely.
    if kill_limit is not None and evaluation_count > kill_limit:
        print("\n" + "!"*60 + f"\n!!! EVAL LIMIT ({kill_limit}) REACHED - TERMINATING !!!\n" + "!"*60)
        raise EvalLimitReached(f"Evaluation limit {kill_limit} reached")
 
    evals_left = (kill_limit - evaluation_count) if kill_limit else '?'
    mins_left  = round(evals_left * 1.07, 1) if isinstance(evals_left, int) else '?'
    print(f"\n--- Eval #{evaluation_count} | ~{mins_left} mins left ---")
    print(f"PARAMS (6 Chans): {[round(v, 4) for v in x]}")
 
    # Run simulation
    update_model(x)
    sim.run(v_init=v_init)
    trace_graph.erase(); curr_v_vec.line(trace_graph, h.dt, 2, 1); trace_graph.flush()
 
    # Signal processing
    soma_v = np.array(sim.get_rec('soma_v')['vec'])
    #mean subtrace from baseline window (100ms before TBS)
    v_norm = soma_v - np.mean(soma_v[:int(100.0/dt)])
 
    # Feature extraction
    # Snapshot spike times from NEURON's NetCon.record buffer — needed for:
    # (a) ISIs (timing-based penalties) and (b) Activity gating (hard-limit checks)
    sim_spikes = np.array(spike_vec.to_python())
 
    # FIX: calculate_cycle_features now returns three plain lists; was (list, dict-of-lists,
    # dict-of-lists) which required awkward [0] unwrapping at every call site.
    plateau_areas, spikes_per_cycle, troughs_per_cycle = calculate_cycle_features(v_norm, dt)
 
    # ========================================================================
    # ERROR CALCULATION - All normalized to comparable scales
    # ========================================================================
 
    # 1. AREA ERROR - Plateau potential magnitude
    # Normalize by target to make errors scale-invariant
    err_area = 0.0
    for a, t in zip(plateau_areas, target_plateau_per_cycle):
        if abs(t) < 1e-6:
            err_area += (a - t)**2
        else:
            err_area += ((a - t) / abs(t))**2
 
    # 2. SPIKE COUNT ERROR - Must match experimental burst pattern
    err_spks = 0.0
    for s, t in zip(spikes_per_cycle, target_spikes_per_cycle):
        if t < 1:
            err_spks += (s - t)**2
        else:
            err_spks += ((s - t) / t)**2
 
    # 3. TROUGH ERROR - Progressive depolarization between bursts
    # CRITICAL: Target shows troughs getting progressively MORE depolarized (less negative).
    # diff < 0: measured is MORE HYPERPOLARIZED than target — penalize hard
    # diff > 0: measured is LESS HYPERPOLARIZED than target — gentle penalty
    err_trough = 0.0
    for tr, t in zip(troughs_per_cycle, target_trough_per_cycle):
        diff = tr - t
        norm_diff = diff / abs(t) if abs(t) >= 1e-6 else diff
        multiplier = TROUGH_ASYM_HARD if diff < 0 else TROUGH_ASYM_SOFT
        err_trough += (norm_diff**2) * multiplier
 
    # ========================================================================
    # FINAL ERROR COMPILATION
    # ========================================================================
    total_error = (
        (err_spks   * W_SPKS) 
        + (err_area   * W_AREA)
        + (err_trough * W_TROUGH)) 
    
    error_history.append(total_error)
 
    # FIX: Use persistent writer — was opening/closing the file on every evaluation
    _log_writer.writerow(
        [OPT_SEED, evaluation_count, round(total_error, 6)]
        + [round(v, 6) for v in x]
        + [round(a, 6) for a in plateau_areas]
        + spikes_per_cycle
        + [round(tr, 6) for tr in troughs_per_cycle]
    )
 
    if total_error < best_error:
        best_error = total_error
        best_x     = x.copy()   # FIX: track best params for kill-switch recovery
        # FIX: soma_v.copy() instead of copy.deepcopy(soma_v) — deepcopy walks the
        # full Python object graph; .copy() is a single C-level memcpy (~50-100x faster)
        best_trace = soma_v.copy()
        print(f"*** NEW BEST FOUND: {total_error:.4f} ***")
 
    print(f"\n--- Cycle Features ---")
    print(f"  Plateau AUC (mV·s)     | Model:  {[round(a, 3) for a in plateau_areas]}")
    print(f"  [Ca2+ plateau strength] | Target: {[round(t, 3) for t in target_plateau_per_cycle]}")
    print(f"  AP count / burst        | Model:  {spikes_per_cycle}")
    print(f"  [theta-burst firing]    | Target: {target_spikes_per_cycle}")
    print(f"  Inter-burst Vm (mV)     | Model:  {[round(tr, 3) for tr in troughs_per_cycle]}")
    print(f"  [AHP / slow depol.]     | Target: {[round(t, 3) for t in target_trough_per_cycle]}")

    print(f"\n--- Normalized Errors ---")
    print(f"  Burst AP count error    : {err_spks:.3f}   (0 = correct spikes each theta cycle)")
    print(f"  Plateau AUC error       : {err_area:.3f}   (0 = correct Ca2+ plateau magnitude)")
    print(f"  Inter-burst Vm error    : {err_trough:.3f}   (0 = correct progressive depolarization)")

 
    return total_error


# Wrappers for each optimization stage with appropriate kill limits
def objective_bh(x):
    return objective(x, kill_limit=MAX_EVALS_BH)

def objective_nm(x):
    return objective(x, kill_limit=None)


# ==================================================================================================
# LIVE MONITORING SETUP
# ==================================================================================================
trace_graph = h.Graph(0)
trace_graph.size(0, duration, -80, 40)
trace_graph.view(0, -80, duration, 120, 100, 100, 400, 300)

curr_v_vec = sim.get_rec('soma_v')['vec']

# ==================================================================================================
# 7. OPTIMIZATION - Basin-Hopping + Nelder-Mead polish
# ==================================================================================================

# p_bounds = [
#     (0.4, 0.8),   # cal
#     (0.9, 1.1),   # calH
#     (0.9, 1.4),   # car
#     (0.9, 1.4),   # cat
#     (1.0, 1.8),   # mykca — narrowed to intermediate range
#     (1.0, 1.8)    # kca   — narrowed to intermediate range
# ]

p_bounds = [
    (0.35, 0.55),  # cal   — best was 0.4, give it room either side
    (1.0,  1.2),   # calH  — best was 1.1
    (1.3,  1.5),   # car   — best was 1.4
    (1.1,  1.35),  # cat   — best was 1.22
    (1.0,  1.25),  # mykca — best was 1.12
    (1.35, 1.65),  # kca   — best was 1.49
]

# Starting parameters back-calculated from:
# These are the multipliers relative to the base conductances in update_calcium_and_potassium()
# x_start is fixed to the known-good initial values; OPT_SEED only randomizes basinhopping jump directions.
# x_start = np.array([0.599473, 1.033540, 1.119427, 1.136134, 1.300000, 1.300000])
x_start = np.array([0.4, 1.155, 1.4, 1.222622, 1.121581, 1.492279])
print(f"x_start (fixed): {[round(v, 4) for v in x_start]}, OPT_SEED={OPT_SEED}")

# p_bounds = [
#     (0.1, 1.5),  # cal
#     (0.5, 1.1),  # calH (Strict ceiling)
#     (0.5, 1.5),  # car
#     (0.5, 1.5),  # cat
#     (0.1, 2.0),  # mykca
#     (0.5, 2.0)   # kca
# ]

# def callback(x, f, accept):
#     if f < best_error:
#         checkpoint = {
#             'channels':  {k: float(v) for k, v in zip(['cal', 'calH', 'car', 'cat', 'mykca', 'kca'], x)},
#             'error':     float(f),
#             'timestamp': datetime.datetime.now().strftime("%H:%M:%S")
#         }
#         with open('current_best_checkpoint.yaml', 'w') as file:
#             yaml.dump(checkpoint, file)
#         print(f"--> Checkpoint Saved: Error {f:.4f}")

#Auto-save the best trace to a file
# The 'Auto-Saver': Records the best parameters and plot to disk whenever a new best is found.
# This ensures you never lose the best result if the script crashes or is interrupted.
# This checks if the current trial (f) is better than the best_error. If it is, it immediately over-writes current_best_checkpoint.yaml and current_best_trace.png.
def callback(x, f, accept):
    global best_error
    if f < best_error:
        checkpoint = {
            'opt_seed': OPT_SEED,
            'channels': {k: float(v) for k, v in zip(['cal', 'calH', 'car', 'cat', 'mykca', 'kca'], x)},
            'error': float(f),
            'timestamp': datetime.datetime.now().strftime("%H:%M:%S")
        }
        with open('current_best_checkpoint.yaml', 'w') as file:
            yaml.dump(checkpoint, file)

        if best_trace is not None:
            time_axis = np.arange(0, duration, dt)
            plt.figure(figsize=(12, 4))
            plt.plot(time_axis, best_trace[:len(time_axis)], color='black', linewidth=1.0)
            plt.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
            plt.title(f'Best Trace So Far | Error: {f:.1f} | Eval #{evaluation_count}')
            plt.xlabel('Time (ms)'); plt.ylabel('Voltage (mV)')
            plt.tight_layout()
            plt.savefig('current_best_trace.png', dpi=150)
            plt.close()

        print(f"--> Checkpoint Saved: Error {f:.4f}")

# The 'Guardrail': Forces the optimizer to stay within biological bounds during random jumps.
# Built-in basinhopping is 'blind'; this prevents it from trying negative or extreme conductances.
# This implements the 'BoundStep' class, which acts as a 'Guardrail' for the optimizer.
# During the random 'jumps' of the Basin-Hopping algorithm, this class ensures the parameter values never go outside the biologically realistic bounds defined in p_bounds.
class BoundStep:
    def __init__(self, stepsize=0.2, bounds=None):
        self.stepsize = stepsize
        self.bounds   = np.array(bounds)
    def __call__(self, x):
        # Take a random step within 'stepsize'
        step = np.random.uniform(-self.stepsize, self.stepsize, size=x.shape)
        # Clip the result so we never land outside the [min, max] allowed for each channel
        return np.clip(x + step, self.bounds[:, 0], self.bounds[:, 1])

if __name__ == "__main__":

    # ── STAGE 1: BASIN-HOPPING (Global Search) ────────────────────────────────
    # Basin-hopping = global explorer that "jumps" to random spots in parameter space.
    # At each jump, it runs a local Nelder-Mead minimizer to find the best nearby basin floor.
    # NOTE: Nelder-Mead doesn't support bounds natively — BoundStep enforces them during jumps.
    minimizer_kwargs = {
        "method": "Nelder-Mead",
        "options": {
            "maxiter": MAX_ITERS_BH,
            "maxfev":  MAX_ITERS_BH * 6,  # must match the 6 used in MAX_EVALS_BH
            "xatol":    1e-3,
            "fatol":    1e-3,
            "adaptive": True
        }
    }
 
    print(f"\n--- Starting Global Search (Basin-Hopping, Budget: {MAX_EVALS_BH} evals) ---")
    print(f"    Local minimizer: Nelder-Mead | OPT_SEED: {OPT_SEED}")
 
    # Catch EvalLimitReached (custom exception) instead of relying on StopIteration.
    # If the limit fires, recover best_x from the module-level tracker.

    try:
        res = basinhopping(
            objective_bh,
            x_start,
            niter=MAX_JUMPS_BH,
            T=2.0,
            stepsize=0.3,
            take_step=BoundStep(stepsize=0.3, bounds=p_bounds),
            minimizer_kwargs=minimizer_kwargs,
            callback=callback,
            seed=OPT_SEED,
            disp=True
        )
    except EvalLimitReached:
        print(f"\nEval limit reached. Recovering best result (error: {best_error:.4f})")
        # Reconstruct a minimal result-like object so the rest of the script runs unchanged
        class _Result:
            pass
        res     = _Result()
        res.x   = best_x if best_x is not None else x_start
        res.fun = best_error
 
    print(f"\nBasin-Hopping Complete! Best Error: {res.fun:.4f}")

    # ── STAGE 2: NELDER-MEAD POLISH (Fine-Tuning) ─────────────────────────────
    # Takes the best solution from basin-hopping and fine-tunes it.
    # This is a "warm start" — starts from res.x (the best point found so far)
    # and does a careful local search with tighter tolerances.
    print(f"\n--- Starting Nelder-Mead Polish (warm-starting from Basin-Hopping best) ---")
    res_nm = minimize(
        objective_nm,                                  # objective function (no kill_limit)
        res.x,                                         # warm-start from basin-hopping best
        method='Nelder-Mead',
        options={
            'maxiter':  MAX_EVALS_NM,                  # max iterations for polish
            'maxfev':   MAX_EVALS_NM,                  # max function evaluations
            'xatol':    1e-4,                          # tighter convergence than stage 1
            'fatol':    1e-4,                          # tighter convergence than stage 1
            'disp':     True,
            'adaptive': True                           # adapt simplex to parameter scales
        }
    )

    print(f"\nNelder-Mead Complete! Final Error: {res_nm.fun:.4f}")
    print(f"Improvement over Basin-Hopping:   {res.fun - res_nm.fun:.4f}")

    res = res_nm  # Use the polished NM result as the final answer


    # ==================================================================================================
    # 8. EXPORT RESULTS
    # ==================================================================================================
    print(f"\nOptimization Complete!")
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    update_model(res.x)
    sim.run(v_init=v_init)
    soma_v_final = np.array(sim.get_rec('soma_v')['vec'])
    v_norm_final = soma_v_final - np.mean(soma_v_final[:int(100.0/dt)])

    final_areas, final_spikes, final_troughs = calculate_cycle_features(v_norm_final, dt)

    print(f"\n--- Final Model Performance ---")
    print(f"Areas  (Cycle 1-5): {[round(a, 4) for a in final_areas]} mV*s")
    print(f"Target (Cycle 1-5): {target_plateau_per_cycle} mV*s")
    print(f"Spikes (Cycle 1-5): {final_spikes}")
    print(f"Target (Cycle 1-5): {target_spikes_per_cycle}")
    print(f"Trough (Cycle 1-5): {[round(t, 4) for t in final_troughs]} mV")
    print(f"Target (Cycle 1-5): {target_trough_per_cycle} mV")

    mech_save_path = f'optimized_biophysics_{timestamp}.yaml'
    cell.export_mech_dict(mech_filename=mech_save_path)

    best_params = {
        'opt_seed': OPT_SEED,
        'synapse_seed': SYNAPSE_SEED,
        'channels': {k: float(v) for k, v in zip(['cal', 'calH', 'car', 'cat', 'mykca', 'kca'], res.x)},
        'synapses': {
            'exc_CA3':   config['parameters']['exc_CA3'],
            'exc_ECIII': config['parameters']['exc_ECIII'],
            'inh_CA3':   config['parameters']['inh_CA3'],
            'inh_ECIII': config['parameters']['inh_ECIII']
        },
        'final_error': float(res.fun)
    }

    with open(f'results_summary_{timestamp}.yaml', 'w') as f:
        yaml.dump(best_params, f)

    plt.figure(figsize=(12, 6))
    time_axis = np.arange(0, duration, dt)
    plt.plot(time_axis, best_trace[:len(time_axis)], color='black', linewidth=1.5, label='Optimized Soma')
    plt.axhline(y=-70, color='gray', linestyle='--', alpha=0.5)
    plt.title(f'Final Optimized Soma Trace (Error: {best_error:.4f})')
    plt.xlabel('Time (ms)'); plt.ylabel('Voltage (mV)')
    plt.legend(); plt.tight_layout()
    plt.savefig(f'best_soma_trace_{timestamp}.png')

    _log_fh.close()
    print(f"Results saved to: {os.getcwd()}")