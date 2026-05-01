from neuron import h
from specify_cells import CA1_Pyr, QuickSim
import yaml
from optimization_utils import *
from plot_utils import *

import random


#fraction_exc_syns = {sec_type: float(exc_syn_count[sec_type]) / float(total_exc_syns) for sec_type in exc_syn_count}
#exc_syn_locs = local_random.sample(exc_syn_locs_by_sec_type[sec_type],
#                                           int(num_exc_syns*fraction_exc_syns[sec_type]))


#inh_syn_locs_by_sec_type = cell.get_inhibitory_syn_locs(sec_type_list=inh_syns_sec_types)
#inh_syn_locs = local_random.sample(inh_syn_locs_by_sec_type[sec_type],
#                                           int(num_inh_syns * fraction_inh_syns[sec_type]))


def assign_exc_and_inh_synapse_stims(cell, exc_syn_locs_by_sec_type, inh_syn_locs_by_sec_type,
                                     exc_syn_types, inh_syn_types, num_exc_syns, 
                                     excitatory_stochastic, num_inh_syns, local_random):
    """
    Place excitatory and inhibitory synapses on the cell according to pathway rules.
    Standardized to distribute synapses among sections (e.g., trunk vs apical vs basal) 
    based on the fraction of available locations, matching simulate_place_cell.py logic.
    
    Returns (stim_exc_syns, stim_inh_syns) – dicts keyed by 'CA3' and 'ECIII'.
    """
    stim_exc_syns = {'CA3': [], 'ECIII': []}
    stim_inh_syns = {'CA3': [], 'ECIII': []}

    # 1. Place Excitatory Synapses
    for pathway, num_to_insert in num_exc_syns.items():
        # Pathway rules: ECIII on tuft; CA3 on trunk/apical/basal
        valid_sec_types = ['tuft'] if pathway == 'ECIII' else ['trunk', 'apical']
        
        # Calculate total valid locations for this pathway to get relative weights
        total_valid_locs = sum(len(exc_syn_locs_by_sec_type.get(st, [])) for st in valid_sec_types)
        
        if total_valid_locs == 0:
            continue

        for sec_type in valid_sec_types:
            locs_in_sec = exc_syn_locs_by_sec_type.get(sec_type, [])
            if not locs_in_sec:
                continue
            
            # Relative fraction within this pathway's compartments
            rel_frac = len(locs_in_sec) / total_valid_locs
            n_at_sec = int(num_to_insert * rel_frac)
            
            # Ensure we don't try to sample more than exist
            n_at_sec = min(n_at_sec, len(locs_in_sec))
            if n_at_sec > 0:
                selected_locs = local_random.sample(locs_in_sec, n_at_sec)
                syn_list = cell.insert_synapses_at_syn_locs(selected_locs, exc_syn_types,
                                                            stochastic=excitatory_stochastic)
                stim_exc_syns[pathway].extend(syn_list)

    # 2. Place Inhibitory Synapses
    for pathway, num_to_insert in num_inh_syns.items():
        # Pathway rules: ECIII on tuft; CA3 on soma/trunk/apical
        valid_sec_types = ['tuft'] if pathway == 'ECIII' else ['soma', 'trunk', 'apical']
        
        total_valid_locs = sum(len(inh_syn_locs_by_sec_type.get(st, [])) for st in valid_sec_types)
        
        if total_valid_locs == 0:
            continue

        for sec_type in valid_sec_types:
            locs_in_sec = inh_syn_locs_by_sec_type.get(sec_type, [])
            if not locs_in_sec:
                continue

            rel_frac = len(locs_in_sec) / total_valid_locs
            n_at_sec = int(num_to_insert * rel_frac)
            
            n_at_sec = min(n_at_sec, len(locs_in_sec))
            if n_at_sec > 0:
                selected_locs = local_random.sample(locs_in_sec, n_at_sec)
                syn_list = cell.insert_synapses_at_syn_locs(selected_locs, inh_syn_types, stochastic=False)
                stim_inh_syns[pathway].extend(syn_list)

    return stim_exc_syns, stim_inh_syns


# ---------------------------------------------------------------------------
# GABA_A_KIN gmax helpers
# ---------------------------------------------------------------------------

def _zero_gabaa_gmax(stim_inh_syns):
    """
    Set gmax to 0 on every GABA_A_KIN target.
    Returns a backup dict {synapse_id: original_gmax} for later restoration.
    """
    saved = {}
    for pathway in stim_inh_syns:
        for syn in stim_inh_syns[pathway]:
            if 'GABA_A_KIN' in syn._syn:
                target = syn._syn['GABA_A_KIN']['target']
                saved[id(syn)] = target.gmax
                target.gmax = 0.0
    return saved


def _restore_gabaa_gmax(stim_inh_syns, saved):
    """Restore gmax values backed up by _zero_gabaa_gmax()."""
    for pathway in stim_inh_syns:
        for syn in stim_inh_syns[pathway]:
            if 'GABA_A_KIN' in syn._syn and id(syn) in saved:
                syn._syn['GABA_A_KIN']['target'].gmax = saved[id(syn)]


# ---------------------------------------------------------------------------
# Public experiment runner
# ---------------------------------------------------------------------------

def run_experiment(sim, experiment_type='Control'):
    """
    Run a simulation on a pre-built (cell, sim) pair.

    Control  – GABA_A_KIN and GABAb both active.
    Gabazine – GABA_A_KIN gmax zeroed before the run, restored afterwards,
               so the identical cell can be reused for the next experiment.

    Parameters
    ----------
    sim : QuickSim
        Built by build_basic_cell().
    experiment_type : str
        'Control' or 'Gabazine'.
    """
    inh_syns = sim.stim_inh_syns
    saved_gmax = {}

    if experiment_type == 'Gabazine':
        saved_gmax = _zero_gabaa_gmax(inh_syns)
        print("Gabazine: GABA_A_KIN gmax → 0  (GABAb still active).")

    sim.run(v_init=sim.v_init)

    if experiment_type == 'Gabazine':
        _restore_gabaa_gmax(inh_syns, saved_gmax)
        print("Gabazine: GABA_A_KIN gmax restored.")

    return saved_gmax


# ---------------------------------------------------------------------------
# Cell / sim builder  (call once)
# ---------------------------------------------------------------------------

def build_basic_cell(config, pathway='Both'):
    """
    Build the CA1 pyramidal cell and wire stimulation — **once**.

    Both GABA_A_KIN and GABAb are inserted at every inhibitory synapse
    location.  Use run_experiment(sim, 'Control') or
    run_experiment(sim, 'Gabazine') to execute individual conditions.

    Parameters
    ----------
    config : dict
        Loaded from default_sim_config.yaml.
    pathway : str
        'Both', 'CA3', or 'ECIII' – which pathway receives stimulation.

    Returns
    -------
    cell : CA1_Pyr
    sim  : QuickSim  (has .v_init, .stim_exc_syns, .stim_inh_syns attached)
    """
    # SYNAPSE_SEED must match optimization runs
    SYNAPSE_SEED = 0
    local_random = random.Random(SYNAPSE_SEED)

    # ------------------------------------------------------------------
    # Build cell
    # ------------------------------------------------------------------
    cell = CA1_Pyr(morph_filename=config['simulation']['morph_filename'],
                   mech_filename=config['simulation']['mech_filename'],
                   full_spines=False)

    # Clear any residual synapses
    for node in cell.tree:
        node.content['synapses'] = []

    # ------------------------------------------------------------------
    # Synapse locations
    # ------------------------------------------------------------------
    exc_syns_sec_types = ['soma', 'trunk', 'apical', 'tuft', 'basal']
    inh_syns_sec_types = ['soma', 'ais', 'trunk', 'apical', 'tuft', 'basal']
    exc_syn_locs_by_sec_type = cell.get_excitatory_syn_locs(sec_type_list=exc_syns_sec_types)
    inh_syn_locs_by_sec_type = cell.get_inhibitory_syn_locs(sec_type_list=inh_syns_sec_types)

    num_exc_syns = config['simulation']['num_exc_syns']
    num_inh_syns = config['simulation']['num_inh_syns']
    excitatory_stochastic = False

    exc_syn_types = ['AMPA_KIN', 'NMDA_KIN5']
    # Always include BOTH receptors; Gabazine zeroes GABA_A_KIN gmax at run-time
    inh_syn_types = ['GABA_A_KIN', 'GABAb']

    # ------------------------------------------------------------------
    # Assign the fraction of synapses to insert for exc and inhibitory
    # ------------------------------------------------------------------
    stim_exc_syns, stim_inh_syns = assign_exc_and_inh_synapse_stims(
        cell, exc_syn_locs_by_sec_type, inh_syn_locs_by_sec_type,
        exc_syn_types, inh_syn_types, num_exc_syns, excitatory_stochastic,
        num_inh_syns, local_random)

    cell.init_synaptic_mechanisms()

    # ------------------------------------------------------------------
    # Optional I80T / GABAb mutation
    # ------------------------------------------------------------------
    with open('data/' + config['simulation']['mech_filename'], 'r') as f:
        biophys_config = yaml.safe_load(f)

    gabab_mut = biophys_config.get('gabab_mutation', None)
    if gabab_mut:
        gmax_mult = gabab_mut['gmax_mult']
        for inh_pathway in stim_inh_syns:
            for syn in stim_inh_syns[inh_pathway]:
                if 'GABAb' in syn._syn:
                    syn._syn['GABAb']['target'].gmax *= gmax_mult
        print(f"I80T mutation applied: GABAb gmax × {gmax_mult}")
    else:
        print("No gabab_mutation field found — running as WT")

    n_exc = sum(len(v) for v in stim_exc_syns.values())
    n_inh = sum(len(v) for v in stim_inh_syns.values())
    print(f"Assigned {n_exc} Exc and {n_inh} Inh synapses "
          f"(each Inh site has GABA_A_KIN + GABAb).")

    # ------------------------------------------------------------------
    # QuickSim setup
    # ------------------------------------------------------------------
    # Timing layout (Combined Rin Test + 300 ms ISI):
    #   0–250 ms  : equilibrate
    #   270–370 ms: Rin IClamp (-50 pA, 100 ms)
    #   370–390 ms: 20 ms baseline recovery
    #   390 ms    : EPSP 1 (20ms after Rin ends)
    #   690 ms    : EPSP 2 (+300 ms ISI)
    #   990 ms    : EPSP 3 (+300 ms ISI)
    #   990–1190 ms: 200 ms GABAb / post-EPSP tail
    #   Total      = 1190 ms
    equilibrate    = config['simulation']['equilibrate']   # 250 ms
    opt_stim_times = [
        equilibrate + 170.0,   # 420 ms (50 ms after Rin ends, providing recovery + baseline)
        equilibrate + 470.0,   # 720 ms  (+300 ms ISI)
        equilibrate + 770.0,   # 1020 ms (+300 ms ISI)
    ]
    opt_duration = equilibrate + 970.0   # 1220 ms total

    dt     = config['simulation']['dt']
    v_init = config['simulation']['v_init']

    sim = QuickSim(opt_duration, cvode=False, dt=dt, verbose=0)
    sim.parameters['equilibrate'] = equilibrate
    sim.parameters['duration']    = opt_duration - equilibrate
    sim.parameters['stim_dt']     = dt
    sim.v_init         = v_init
    sim.opt_stim_times = opt_stim_times   # exposed for get_features()

    # Soma recording
    sim.append_rec(cell, cell.tree.root, description='soma', loc=0.5)

    # Input Resistance IClamp — 270 ms onset, 100 ms duration
    rin_start  = config['simulation'].get('input_resistance_test_start',    270.0)
    rin_dur    = config['simulation'].get('input_resistance_test_duration',  100.0)
    rin_amp_pa = config['simulation'].get('input_resistance_test_amp',       -50.0)
    sim.append_stim(cell, cell.tree.root, loc=0.5,
                    amp=rin_amp_pa / 1000.0,
                    delay=rin_start, dur=rin_dur,
                    description='R_in_test')

    # Single shared stim vector for both pathways (same absolute times)
    stim_vec = h.Vector(opt_stim_times)
    sim.opt_stim_vec = stim_vec

    # Connect VecStim sources to synapses
    for group in stim_exc_syns:
        if pathway == 'Both' or pathway == group:
            for syn in stim_exc_syns[group]:
                syn.source.play(stim_vec)

    for group in stim_inh_syns:
        if pathway == 'Both' or pathway == group:
            for syn in stim_inh_syns[group]:
                syn.source.play(stim_vec)

    # Backward-compat aliases so _play_pathway() still works
    sim.ec3_vec = stim_vec
    sim.ca3_vec = stim_vec

    # Attach synapse dicts to sim for run_experiment()
    sim.stim_exc_syns = stim_exc_syns
    sim.stim_inh_syns = stim_inh_syns
    sim.local_random = local_random

    print("Cell built successfully.")
    print("  → run_experiment(sim, 'Control')  : GABA_A_KIN + GABAb active")
    print("  → run_experiment(sim, 'Gabazine') : GABA_A_KIN gmax=0, GABAb active")
    return cell, sim