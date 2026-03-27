import os
import sys
import yaml
import numpy as np
import random
import time
import datetime
import copy
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from neuron import h, gui
from specify_cells import CA1_Pyr, QuickSim
import csv 

# ==================================================================================================
# 1. SETUP & TARGETS (Dendritic Calibration)
# ==================================================================================================

# --- YOUR GOLDEN SOMATIC BASELINE (From Eval 29) ---
# These stay frozen to ensure the soma remains stable
FIXED_SOMA_MULTS = {
    'm_cal': 0.4714, 'm_calH': 1.0241, 'm_car': 1.0653, 
    'm_cat': 1.0844, 'm_mykca': 1.1891, 'm_kca': 1.0622
}

# ==================================================================================================
# 1. UPDATED EXPERIMENTAL TARGETS (from Milstein et al 2015 & Takahashi and Magee 2009)
# ==================================================================================================
target_dend_peak_v = 50.0       # mV above rest (The 'Full Plateau' threshold)
target_summation_rise = 7.5     # mV baseline climb per burst (Temporal Summation)
target_plateau_duration = 150.0 # ms (Total time above 20mV threshold)

# Load Config
config_path = 'optimization_config.yaml'
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)

duration = config['simulation']['equilibrate'] + config['simulation']['sim_duration']
dt = config['simulation']['dt']
v_init = config['simulation']['v_init']
tbs_times = config['simulation']['tbs_times']

evaluation_count = 0
best_error = float('inf')
best_trace = None
error_history = []

# ==================================================================================================
# 2. LOGGING & FILENAMES
# ==================================================================================================
start_time_str = datetime.datetime.now().strftime("%Y%m%d_%H%M")
log_filename = f'dendritic_opt_log_{start_time_str}.csv'

# Optimizing 3 key dendritic drivers: 1. Apical calH, 2. Tuft calH, 3. Dendritic kca
headers = ['Eval', 'Total_Error', 'apic_calH', 'tuft_calH', 'dend_kca']
with open(log_filename, 'w', newline='') as f:
    writer = csv.writer(f); writer.writerow(headers)

# ==================================================================================================
# 3. CELL & RECORDING INITIALIZATION (Fixed for Dendrite)
# ==================================================================================================
print("--- Initializing Cell for Distal Dendritic Recording ---")
cell = CA1_Pyr(morph_filename=config['simulation']['morph_filename'], 
               mech_filename=config['simulation']['mech_filename'], 
               full_spines=False)

# POINT RECORDING TO APICAL DENDRITE (e.g., 250um from Soma to match experiment)
dend_node = cell.get_node_by_distance(250.0, sec_type='apical')
sim = QuickSim(duration, cvode=False, dt=dt, verbose=0)
sim.append_rec(cell, dend_node, description='dend_v') 

spike_vec = h.Vector()
cell.spike_detector.record(spike_vec)

# Synapse Setup
exc_syn_locs_by_sec_type = cell.get_excitatory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
inh_syn_locs_by_sec_type = cell.get_inhibitory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
exc_syn_types, inh_syn_types = ['AMPA_KIN', 'NMDA_KIN5'], ['GABAb']
tbs_vec = h.Vector(tbs_times)
stim_exc_syns, stim_inh_syns = {'CA3': [], 'ECIII': []}, {'CA3': [], 'ECIII': []}
local_random = random.Random(0)

# Live Monitor (Scaled for large dendritic voltages)
trace_graph = h.Graph(0)
trace_graph.size(0, duration, -80, 40)
trace_graph.view(0, -80, duration, 120, 100, 100, 400, 300)
curr_v_vec = sim.get_rec('dend_v')['vec']

# ==================================================================================================
# 4. UPDATE LOGIC (Fixed: Soma Frozen, Dendrites Scaled)
# ==================================================================================================

def update_model_dendritic(x_vec):
    m_apic_calH, m_tuft_calH, m_dend_kca = x_vec
    
    # A. Apply Fixed Soma Base (using the update function from your original code)
    # We use somatic baseline for ALL segments first, then override
    cell.modify_mech_param('soma', 'cal', 'gcalbar', value=0.007 * FIXED_SOMA_MULTS['m_cal'])
    cell.modify_mech_param('soma', 'car', 'gcabar', value=0.003 * FIXED_SOMA_MULTS['m_car'])
    cell.modify_mech_param('soma', 'cat', 'gcatbar', value=0.00005 * FIXED_SOMA_MULTS['m_cat'])
    cell.modify_mech_param('soma', 'mykca', 'gkbar', value=0.09075 * FIXED_SOMA_MULTS['m_mykca'])
    cell.modify_mech_param('soma', 'kca', 'gbar', value=0.0005 * FIXED_SOMA_MULTS['m_kca'])
    
    # B. Override the Apical/Tuft (This is the Dendritic Polish)
    cell.modify_mech_param('apical', 'calH', 'gcalbar', value=0.001455 * m_apic_calH)
    cell.modify_mech_param('tuft', 'calH', 'gcalbar', value=0.001455 * m_tuft_calH)
    
    for sec in ['apical', 'tuft']:
        cell.modify_mech_param(sec, 'kca', 'gbar', value=4.167e-05 * m_dend_kca)

    # C. Re-init and Synapses (Constants from YAML)
    for node in cell.tree: node.content['synapses'] = []
    for pathway in stim_exc_syns: stim_exc_syns[pathway] = []
    
    num_exc = {'CA3': config['parameters']['num_exc_CA3'], 'ECIII': config['parameters']['num_exc_ECIII']}
    for pathway, num in num_exc.items():
        locs = exc_syn_locs_by_sec_type['tuft'] if pathway == 'ECIII' else exc_syn_locs_by_sec_type['apical']
        syns = cell.insert_synapses_at_syn_locs(local_random.sample(locs, min(num, len(locs))), exc_syn_types)
        stim_exc_syns[pathway].extend(syns)
    
    cell.init_synaptic_mechanisms()
    for pathway in stim_exc_syns:
        for syn in stim_exc_syns[pathway]: syn.source.play(tbs_vec)

# ==================================================================================================
# 5. DENDRITIC OBJECTIVE FUNCTION (Fixed for Summation & Peaks)
# ==================================================================================================

def objective(x):
    global best_error, best_trace, evaluation_count
    evaluation_count += 1
    if evaluation_count > 150: sys.exit()

    update_model_dendritic(x)
    sim.run(v_init=v_init)
    trace_graph.erase(); curr_v_vec.line(trace_graph, h.dt, 2, 1); trace_graph.flush()

    dend_v = np.array(sim.get_rec('dend_v')['vec'])
    v_norm = dend_v - np.mean(dend_v[:int(100.0/dt)])

    # 1. SUMMATION ERROR (Baseline rise across 5 bursts)
    burst_starts = [500, 650, 800, 950, 1100]
    summation_err = 0
    for i, start in enumerate(burst_starts):
        pre_burst_v = v_norm[int((start - 5)/dt)] # Look 5ms before next burst
        target_v = i * target_summation_rise
        summation_err += (pre_burst_v - target_v)**2

    # 2. PEAK ERROR (Goal: ~45mV above rest)
    peak_err = (np.max(v_norm) - target_dend_peak_v)**2

    # 3. ANTI-FLATLINE (Ensures we keep dendritic spikes)
    dv = np.diff(v_norm[int(500/dt):int(1250/dt)])
    activity_err = 100.0 if np.var(dv) < 0.2 else 0.0

    total_error = (summation_err * 0.2) + (peak_err * 0.5) + activity_err
    
    with open(log_filename, 'a', newline='') as f:
        csv.writer(f).writerow([evaluation_count, total_error] + list(x))

    if total_error < best_error:                        
        best_error, best_trace = total_error, copy.deepcopy(dend_v)
        print(f"*** NEW BEST DENDRITE: {total_error:.4f} ***")
    
    print(f"Eval {evaluation_count} | Peak: {np.max(v_norm):.1f}mV | SumErr: {summation_err:.1f}")
    return total_error

# ==================================================================================================
# 6. RUN DENDRITIC OPTIMIZATION
# ==================================================================================================
# Starting: [apic_calH, tuft_calH, dend_kca]
x_start = np.array([2.5, 3.5, 0.8]) 
p_bounds = [(1.0, 15.0), (1.0, 15.0), (0.01, 2.0)]

print("\n--- Starting Final Dendritic Polish (Caffeinate Active) ---")
res = minimize(objective, x_start, method='L-BFGS-B', bounds=p_bounds, 
               options={'maxiter': 8, 'maxfun': 60, 'disp': True})

# ==================================================================================================
# 7. EXPORT RESULTS
# ==================================================================================================
print(f"\nOptimization Complete!")
out_dir = "optimized_config_dendrite"
if not os.path.exists(out_dir): os.makedirs(out_dir)
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

update_model_dendritic(res.x)
cell.export_mech_dict(mech_filename=f'./{out_dir}/dendrite_biophysics_{timestamp}.yaml')

plt.figure(figsize=(10, 5))
plt.plot(np.arange(0, duration, dt), best_trace[:int(duration/dt)], color='black')
plt.title(f'Optimized Dendritic Trace (Err: {best_error:.4f})')
plt.savefig(f'./{out_dir}/dendrite_best_trace_{timestamp}.png')