from neuron import h
from specify_cells import CA1_Pyr, QuickSim, data_dir
from plot_utils import *
from optimization_utils import calculate_cycle_features
import numpy as np
import os, sys, gc, yaml, random
from nested.utils import Context, str_to_bool
from nested.optimize_utils import *

context = Context()

# --- CONSTANTS (Must match simulation protocol) ---
BURST_STARTS_MS   = [500.0, 650.0, 800.0, 950.0, 1100.0]
BURST_DURATION_MS = 150.0
TROUGH_WINDOW_MS  = 50.0

def config_worker():
    """Load configuration and initialize context."""
    context.config_dict = read_from_yaml(context.model_config_file_path)
    # Default parameters that aren't optimized
    context.equilibrate = context.config_dict['simulation']['equilibrate']
    context.duration = context.equilibrate + context.config_dict['simulation']['sim_duration']
    context.dt = context.config_dict['simulation']['dt']
    context.v_init = context.config_dict['simulation']['v_init']
    context.tbs_times = context.config_dict['simulation']['tbs_times']
    context.targets = context.config_dict['targets']

def get_seeds():
    return [list(range(context.seed_start, context.seed_start + context.num_seeds))]

def update_model_parameters(cell, params):
    """Map 'actual values' from optimizer to cell mechanisms."""
    # 1. cal (L-type)
    cell.modify_mech_param('soma', 'cal', 'gcalbar', value=params['cal_gcalbar_soma'])
    cell.modify_mech_param('trunk', 'cal', 'gcalbar', value=params['cal_gcalbar_trunk'])
    cell.modify_mech_param('apical', 'cal', 'gcalbar', value=params['cal_gcalbar_apical'])
    cell.modify_mech_param('tuft', 'cal', 'gcalbar', value=params['cal_gcalbar_tuft'])
    cell.modify_mech_param('basal', 'cal', 'gcalbar', origin='soma')

    # 2. calH (High-Threshold L-type)
    cell.modify_mech_param('soma', 'calH', 'gcalbar', value=params['calH_gcalbar_soma'])
    cell.modify_mech_param('trunk', 'calH', 'gcalbar', value=params['calH_gcalbar_soma'], max_loc=50.0)
    cell.modify_mech_param('trunk', 'calH', 'gcalbar', value=params['calH_gcalbar_trunk'], min_loc=50.0, replace=False)
    cell.modify_mech_param('apical', 'calH', 'gcalbar', value=params['calH_gcalbar_apical'])
    cell.modify_mech_param('tuft', 'calH', 'gcalbar', value=params['calH_gcalbar_tuft'])

    # 3. car (R-type)
    cell.modify_mech_param('soma', 'car', 'gcabar', value=params['car_gcabar_soma'])
    cell.modify_mech_param('trunk', 'car', 'gcabar', value=params['car_gcabar_trunk'])
    cell.modify_mech_param('apical', 'car', 'gcabar', value=params['car_gcabar_apical'])
    cell.modify_mech_param('tuft', 'car', 'gcabar', value=params['car_gcabar_tuft'])
    cell.modify_mech_param('basal', 'car', 'gcabar', origin='soma')

    # 4. cat (T-type)
    cell.modify_mech_param('soma', 'cat', 'gcatbar', value=params['cat_gcatbar_soma'])
    cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0, max_loc=100.0)
    # Calculate slope: reaches cat_gcatbar_trunk by 350um
    slope_cat = params['cat_gcatbar_trunk'] / 250.0 if params['cat_gcatbar_trunk'] > 0 else 0
    cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0, slope=slope_cat, min_loc=100.0, max_loc=350.0, replace=False)
    cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=params['cat_gcatbar_trunk'], min_loc=350.0, replace=False)
    cell.modify_mech_param('apical', 'cat', 'gcatbar', value=params['cat_gcatbar_apical'])
    cell.modify_mech_param('tuft', 'cat', 'gcatbar', value=params['cat_gcatbar_tuft'])
    cell.modify_mech_param('basal', 'cat', 'gcatbar', origin='soma')

    for sec in ['soma', 'trunk', 'basal', 'apical', 'tuft']:
        for mech in ['cal', 'calH', 'car', 'cat', 'pas']:
            cell.reinitialize_subset_mechanisms(sec, mech)

def run_evaluation(params, seed, model_id, export=False, plot=False):
    """The single-evaluation loop called by the optimizer."""
    # Initialize Cell
    cell = CA1_Pyr(morph_filename=context.config_dict['simulation']['morph_filename'],
                   mech_filename=context.config_dict['simulation']['mech_filename'],
                   full_spines=False)
    
    # Map array to dictionary
    param_dict = param_array_to_dict(params, context.param_names)
    
    # Update Cell Biophysics
    update_model_parameters(cell, param_dict)
    
    # Setup Synapses
    local_random = random.Random(seed)
    exc_syn_locs = cell.get_excitatory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
    inh_syn_locs = cell.get_inhibitory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
    
    num_exc = {'CA3': int(param_dict['num_exc_syns_CA3']), 'ECIII': int(param_dict['num_exc_syns_ECIII'])}
    num_inh = {'CA3': int(param_dict['num_inh_syns_CA3']), 'ECIII': int(param_dict['num_inh_syns_ECIII'])}
    
    # Place synapses
    stim_exc_syns = {'CA3': [], 'ECIII': []}
    stim_inh_syns = {'CA3': [], 'ECIII': []}
    
    for pathway in ['CA3', 'ECIII']:
        # Exc
        locs = exc_syn_locs['tuft'] if pathway == 'ECIII' else exc_syn_locs['trunk'] + exc_syn_locs['apical']
        sel = local_random.sample(locs, min(num_exc[pathway], len(locs)))
        stim_exc_syns[pathway] = cell.insert_synapses_at_syn_locs(sel, ['AMPA_KIN', 'NMDA_KIN5'], stochastic=True)
        # Inh
        locs_inh = inh_syn_locs['tuft'] if pathway == 'ECIII' else inh_syn_locs['trunk'] + inh_syn_locs['apical']
        sel_inh = local_random.sample(locs_inh, min(num_inh[pathway], len(locs_inh)))
        stim_inh_syns[pathway] = cell.insert_synapses_at_syn_locs(sel_inh, ['GABAb'], stochastic=False)
    
    # Apply GABAb weight
    for pathway in stim_inh_syns:
        for syn in stim_inh_syns[pathway]:
            syn._syn['GABAb']['target'].gmax = param_dict['gabab_gmax']

    cell.init_synaptic_mechanisms()
    
    # Play TBS times
    tbs_vec = h.Vector(context.tbs_times)
    for group in list(stim_exc_syns.values()) + list(stim_inh_syns.values()):
        for syn in group:
            syn.source.play(tbs_vec)

    # Run Session
    sim = QuickSim(context.duration, cvode=True, dt=context.dt, verbose=0)
    sim.append_rec(cell, cell.tree.root, description='soma_v')
    spike_vec = h.Vector()
    cell.spike_detector.record(spike_vec)
    
    sim.run(v_init=context.v_init)
    
    # Calculate Features
    soma_v = np.array(sim.get_rec('soma_v')['vec'])
    v_norm = soma_v - np.mean(soma_v[:int(100.0/context.dt)])
    areas, spikes, troughs = calculate_cycle_features(v_norm, context.dt)
    
    features = {
        'plateau_areas': areas,
        'spike_counts': spikes,
        'trough_vms': troughs
    }
    
    if plot:
        plt.plot(sim.tvec, soma_v)
        plt.show()

    # Cleanup
    del cell, sim, stim_exc_syns, stim_inh_syns
    gc.collect()
    
    return features

def calculate_objective(features, model_id, export=False, plot=False):
    """Compute normalized errors for each target feature."""
    targets = context.targets
    
    # 1. Area Error
    err_area = sum(((a - t) / t)**2 if t != 0 else (a-t)**2 
                   for a, t in zip(features['plateau_areas'], targets['target_plateau_per_cycle']))
    
    # 2. Spike Error
    err_spks = sum(((s - t) / t)**2 if t != 0 else (s-t)**2 
                   for s, t in zip(features['spike_counts'], targets['target_spikes_per_cycle']))
    
    # 3. Trough Error
    err_trough = sum(((tr - t) / abs(t))**2 if t != 0 else (tr-t)**2 
                     for tr, t in zip(features['trough_vms'], targets['target_trough_per_cycle']))
    
    total_error = (err_area * 800.0) + (err_spks * 1000.0) + (err_trough * 1200.0)
    
    objectives = {
        'total_error': total_error,
        'plateau_area_error': err_area,
        'spike_count_error': err_spks,
        'trough_depol_error': err_trough
    }
    
    return features, objectives