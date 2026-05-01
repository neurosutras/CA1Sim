import os
os.environ['NRN_NO_MPI'] = '1'
import matplotlib
matplotlib.use('Agg')
import sys
import numpy as np
import yaml
import time
import datetime
import random
import gc
import pickle
from collections import defaultdict
from neuron import h
from specify_cells import CA1_Pyr, QuickSim, data_dir, Synapse
from optimization_utils import calculate_cycle_features
from nested.parallel import *
from nested.optimize_utils import *
import click
import scipy.optimize
from basic_sim_utils import assign_exc_and_inh_synapse_stims
from nested.optimize_utils import nested_optimize_init_controller_context


# Global context for the nested framework
context = Context()


def log10_fit(x, a, b, c):
    """Logarithmic fit for f-I curves: f(i) = a * log10(i - b) + c"""
    return a * np.log10(np.maximum(1e-9, x - b)) + c

def inverse_log10_fit(y, a, b, c):
    """Inverse of log10_fit to find current (i) for a given frequency (f)."""
    return 10**((y - c) / a) + b



@click.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True, ))
@click.option("--config-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False),
              default='config/Nested_optimize_CA1_tbs_e_I_spiking_config.yaml')
@click.option("--mech-filename", type=str, default='20260426_default_biophysics_with_Poirazi_Ca_Ch_and_2018_Kca.yaml')
@click.option("--output-dir", type=click.Path(exists=False, file_okay=False, dir_okay=True), default='data')
@click.option("--export", is_flag=True)
@click.option("--export-file-path", type=str, default=None)
@click.option("--label", type=str, default=None)
@click.option("--verbose", type=int, default=2)
@click.option("--plot", is_flag=True)
@click.option("--interactive", is_flag=True)
@click.option("--debug", is_flag=True)
@click.option("--framework", type=str, default='serial')
def main(config_file_path, output_dir, export, export_file_path, label, verbose, plot, interactive, debug, framework):

    """

    :param config_file_path: Path to the config file.
    :param output_dir: Directory to save the output.
    :param export: Whether to export the results.
    :param export_file_path: Path to the export file.
    :param label: Label for the optimization.
    :param verbose: Whether to print verbose output.
    :param plot: Whether to plot the results.
    :param interactive: Whether to run in interactive mode.
    :param debug: Whether to run in debug mode.
    """

    # requires a global variable context: :class:'Context'
    kwargs = get_unknown_click_arg_dict(click.get_current_context().args)
    context.update(locals())
    context.kwargs = kwargs
    context.verbose = verbose
    context.plot = plot
    context.debug = debug
    context.disp = verbose > 0

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    context.interface = get_parallel_interface(framework=framework, source_file=__file__, source_package=__package__, **kwargs)
    context.interface.start(disp=context.disp)
    nested_optimize_init_controller_context(context, config_file_path=config_file_path, 
                                            output_dir=output_dir, label=label, disp=context.disp, 
                                            **kwargs)

    # Force synchronization to ALL workers using a map operation
    # We use an explicit local function to ensure the correct global 'context' is updated
    if hasattr(context.interface, 'num_workers') and context.interface.num_workers > 0:
        context.interface.map(update_local_context, [get_picklable_context()] * context.interface.num_workers)
    
    if not debug:
        context.interface.ensure_controller()
        run_sim_tests()
        if not interactive:
            context.interface.stop()


def get_trunk_bifurcation(cell):
    """
    Locate the biological trunk bifurcation point instead of defaulting to the proximal trunk.
    """
    if getattr(cell, 'trunk', None) is None:
        return None
    trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
    if trunk_bifurcation:
        trunk_branches = [branch for branch in trunk_bifurcation[0].children if branch.type == 'trunk']
        if trunk_branches:
            # get where the thickest trunk branch gives rise to the tuft
            trunk = max(trunk_branches, key=lambda node: node.sec(0.).diam)
            trunk = next((node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type for child in node.children)), trunk_bifurcation[0])
        else:
            trunk = trunk_bifurcation[0]
    else:
        trunk_bifurcations = [node for node in cell.trunk if 'tuft' in (child.type for child in node.children)]
        trunk = trunk_bifurcations[0] if trunk_bifurcations else None
    return trunk

def update_local_context(content):
    """
    Explicitly update the global context on the worker.
    """
    global context
    if content is not None:
        context.update(content)


def plot_sim_traces(sim, target_descriptions, save_path):
    """
    Custom plotter for QuickSim that filters by recorded descriptions.
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for rec_dict in sim.rec_list:
        desc = rec_dict.get('description', '')
        if desc in target_descriptions:
            ax.plot(sim.tvec, rec_dict['vec'], label=f"{rec_dict['node'].name}({rec_dict['loc']}) - {desc}")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Vm (mV)")
    ax.legend(loc='upper right')
    plt.savefig(save_path)
    plt.close(fig)

def config_controller():
    """
    Initialization hook for the master process (controller).
    Ensures that simulation parameters and target arrays are initialized on the master
    so that mapping functions like get_args_dynamic_fi can access them.
    """
    init_context()

def config_worker():
    """

    Configure the worker for the optimization - 
    """

    if 'verbose' in context():
        context.verbose = int(context.verbose)
    if 'plot' not in context():
        context.plot = False
    if not context_has_sim_env(context):
        build_sim_env(context, **context.kwargs)
    else:
        config_sim_env(context)
    
def context_has_sim_env(context):
    """
    Check if the context has a sim env
    """
    return 'env' in context() and 'sim' in context() and 'cell' in context()

def load_yaml(file_path):
    with open(file_path, 'r') as f:
        return yaml.load(f, Loader=yaml.Loader)

def init_context(default_sim_config_path='default_sim_config.yaml'):
    context.verbose = getattr(context, 'verbose', 0)
    if getattr(context, 'verbose', 0) > 0:
        print(f"DEBUG: [Master/Worker] Initializing context from {default_sim_config_path}...")

    #Basic sim parameters
    loaded_config = load_yaml(default_sim_config_path)
    if 'simulation' in loaded_config:
        default_sim_config = loaded_config['simulation']
    else:
        default_sim_config = loaded_config
        
    equilibrate = context().get('equilibrate', default_sim_config.get('equilibrate', 250.0))
    baseline = context().get('baseline', default_sim_config.get('baseline', 10.0))
    sim_duration = context().get('sim_duration', default_sim_config.get('sim_duration', 1500.0))
    dt = context().get('dt', default_sim_config.get('dt', 0.025))
    v_init = default_sim_config.get('v_init', -70.0)
    morph_filename = context().get('morph_filename', default_sim_config.get('morph_filename', 'EB2-late-bifurcation.swc'))
    mech_filename = context().get('mech_filename', default_sim_config.get('mech_filename', '20260426_default_biophysics_with_Poirazi_Ca_Ch_and_2018_Kca.yaml'))

    #input resistance test parameters
    input_resistance_test_start = 260.0 # 250 + 10 baseline
    input_resistance_test_duration = default_sim_config.get('input_resistance_test_duration', 100.0)
    input_resistance_test_amp = default_sim_config.get('input_resistance_test_amp', -50.0)
    
    #f-I curve test parameters
    fi_test_start = 260.0
    fi_test_duration = default_sim_config.get('fi_test_duration', 500.0)
    fi_test_amps = default_sim_config.get('fi_test_amps', [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0])

    #rheobase test parameters
    rheobase_test_start = 260.0
    rheobase_test_duration = default_sim_config.get('rheobase_test_duration', 500.0)
    rheobase_test_amps = default_sim_config.get('rheobase_test_amps', [0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0])

    #synaptic stimulation parameters
    ECIII_epsp_times = [260.0, 560.0, 860.0]
    CA3_epsp_times = [260.0, 560.0, 860.0]

    #theta burst stimulation parameters
    burst_starts_ms = [260.0, 410.0, 560.0, 710.0, 860.0]
    burst_duration_ms = default_sim_config.get('burst_duration_ms', 150.0)
    trough_window_ms = default_sim_config.get('trough_window_ms', 50.0)
    theta_stim_win_start_ms = 260.0
    theta_stim_win_end_ms = 1160.0
    pulses_per_burst = context().get('pulses_per_burst', default_sim_config.get('pulses_per_burst', 5))
    intra_burst_interval_ms = context().get('intra_burst_interval_ms', default_sim_config.get('intra_burst_interval_ms', 10.0))
    num_exc_syns = default_sim_config.get('num_exc_syns', None)
    num_inh_syns = default_sim_config.get('num_inh_syns', None)
    excitatory_stochastic = default_sim_config.get('excitatory_stochastic', True)

    # Default relative amplitudes (in nA) if not determined dynamically from targets
    i_inj_increment_f_I = 0.050
    num_increments_f_I = len(fi_test_amps)
    i_inj_relative_amp_array = np.array([i_inj_increment_f_I * i for i in range(1, num_increments_f_I + 1)])

    # For fitting FI curves from experimental targets dynamically
    if 'target_val' in context():
        targets = context.target_val
        exp_rates = []
        exp_amps = []
        
        for amp in fi_test_amps:
            if amp > 0:  # Ignore 0.0 current injection for the log fit
                key = f'fi_rate_{int(amp)}'
                if key in targets:
                    exp_amps.append(float(amp))
                    exp_rates.append(float(targets[key]))
                    
        if exp_rates:
            exp_i_inj_amp_f_I_0 = np.array(exp_amps) / 1000.0  # nA
            exp_rate_f_I_0 = np.array(exp_rates)
            
            fit_params_f_I_0 = [200., 0.05, 0.0]  # Initial guess
            exp_fit_params_f_I, pcov = scipy.optimize.curve_fit(log10_fit, exp_i_inj_amp_f_I_0, exp_rate_f_I_0, fit_params_f_I_0)
            
            rate_at_rheobase = 1.0 # 1 Hz
            exp_rheobase = inverse_log10_fit(rate_at_rheobase, *exp_fit_params_f_I)

            # Determine relative increment sizes straight from the fi_test_amps spacing
            i_inj_increment_f_I = (exp_amps[1] - exp_amps[0]) / 1000.0 if len(exp_amps) > 1 else 0.050
            num_increments_f_I = len(exp_amps)
            
            i_inj_relative_amp_array = np.array([i_inj_increment_f_I * i for i in range(1, num_increments_f_I + 1)])
            exp_i_inj_amp_array = np.add(exp_rheobase, i_inj_relative_amp_array)
            exp_rate_f_I_array = log10_fit(exp_i_inj_amp_array, *exp_fit_params_f_I)

    context.update(locals())
    # Compatibility for nested framework which flattens kwargs
    if not hasattr(context, 'kwargs') or context.kwargs is None:
        context.kwargs = context()

#based from basic_sim_utils.py
def build_sim_env(context, verbose=2, cvode=False, daspk=False, load_edges=False, set_edge_delays=False, **kwargs):
    """
    :param context: :class:'Context' - The global container that stores the cell, simulation, and results.
    :param verbose: int - Level of logging (0-2). High verbosity (2) prints detailed simulation progress.
    :param cvode: bool - Enables Variable Time-Step integration. Used to significantly speed up subthreshold sims 
                 by allowing the simulator to use larger time steps when voltage changes are slow.
    :param daspk: bool - Enables the DASPK solver. Used primarily for 'Consistent Initialization' to ensure the 
                 model starts at a stable resting potential, which is critical for accurate Rin measurements.
    """
    # get parameters from context
    verbose = int(verbose)
    init_context()  # Assumes init_context is called as needed
    
    # context.env = Env(comm=context.comm, verbose=verbose > 1, **kwargs)
    # configure_hoc_env(context.env)

    # Get morphology/mechanisms from config
    config = context.kwargs
    morph_filename = config['morph_filename']
    mech_filename = config['mech_filename']

    # Initialize the CA1 Cell
    context.cell = CA1_Pyr(morph_filename=morph_filename,
                           mech_filename=mech_filename,
                           full_spines=False)

    # Setup spike detector: CA1_Pyr.spike_detector is a pre-built NetCon
    context.spike_output_vec = h.Vector()
    context.cell.spike_detector.record(context.spike_output_vec)

    # Initialize the simulation engine (QuickSim)
    context.sim = QuickSim(context.kwargs.get('sim_duration', 1500.0), 
                           cvode=cvode, 
                           dt=context.kwargs['dt'], 
                           verbose=(verbose > 1))

    # Monkey-patch QuickSim to support legacy Env wrapper queries (by string description)
    context.sim.has_rec = lambda desc: any(r.get('description') == desc for r in context.sim.rec_list)
    context.sim.has_stim = lambda desc: any(s.get('description') == desc for s in context.sim.stim_list)
    
    orig_modify_stim = context.sim.modify_stim
    def modify_stim_by_desc(desc_or_index, **kwargs):
        if isinstance(desc_or_index, str):
            for i, s in enumerate(context.sim.stim_list):
                if s.get('description') == desc_or_index:
                    return orig_modify_stim(index=i, **kwargs)
            raise KeyError(f'Stimulus parameter {desc_or_index} not found')
        else:
            return orig_modify_stim(index=desc_or_index, **kwargs)
    context.sim.modify_stim = modify_stim_by_desc

    # Recording setup (Somatic Rin always measured at the root)
    context.sim.append_rec(context.cell, context.cell.tree.root, description='soma', loc=0.5)

    # Recording setup for the trunk (Dendritic Rin at the bifurcation)
    trunk_node = get_trunk_bifurcation(context.cell)
    if trunk_node:
        context.sim.append_rec(context.cell, trunk_node, description='trunk', loc=0.5)

    # Recording setup for the AIS
    if hasattr(context.cell, 'axon') and len(context.cell.axon) > 1:
        context.sim.append_rec(context.cell, context.cell.axon[1], description='ais', loc=0.5)

    # config_sim_env(context) # Standard framework hook
    
    if verbose > 0:
        print(f"    Cell built from {morph_filename}.")
        print(f"    Sim engine initialized (duration: {context.kwargs.get('sim_duration', 1500.0)}ms, cvode: {cvode}).")
        if trunk_node:
            print(f"    Trunk recording active on node index {trunk_node.index} (bifurcation).")

def config_sim_env(context):
    """
    Configure the characterized simulation environment. This ensures that the necessary 
    recordings and stimuli (IClamps) are attached to the cell for characterization.
    """
    if 'previous_module' in context() and context.previous_module == __file__:
        return
    
    init_context()
    
    # Initialize holding current tracking if not already present
    if 'i_holding' not in context():
        context.i_holding = defaultdict(dict)

    if 'cell' not in context():
        build_sim_env(context, **context.kwargs)
    
    cell = context.cell
    sim = context.sim

    # 1. Ensure Standard Characterization Recordings
    if not sim.has_rec('soma'):
        sim.append_rec(cell=cell, node=cell.tree.root, description='soma', loc=0.5)
    
    if not sim.has_rec('trunk'):
        trunk_node = get_trunk_bifurcation(cell)
        if trunk_node:
            sim.append_rec(cell=cell, node=trunk_node, description='trunk', loc=0.5)

    # 2. Setup Standard Stimulation IClamps
    equilibrate = context.kwargs.get('equilibrate', 250.0)
    stim_dur = context.kwargs.get('input_resistance_test_duration', 100.0) # Used for Rin test
    duration = context.kwargs.get('sim_duration', 1500.0)

    # The 'step' clamp is used for the Rin test (-50pA) and FI curve steps
    if not sim.has_stim('step'):
        sim.append_stim(cell=cell, node=cell.tree.root, loc=0.5, description='step', amp=0., delay=equilibrate, dur=stim_dur)
    
    # The 'holding' clamp can be used to offset the membrane potential if needed
    if not sim.has_stim('holding'):
        sim.append_stim(cell=cell, node=cell.tree.root, loc=0.5, description='holding', amp=0., delay=0., dur=duration)

    # Standardize V_active for this model's baseline
    if 'v_active' not in context():
        context.v_active = context.kwargs.get('v_init', -70.0)
        
    if context.v_active not in context.i_holding['soma']:
        context.i_holding['soma'][context.v_active] = 0.

    # Track that this module has been configured
    context.previous_module = __file__

    verbose = getattr(context, 'verbose', 0)
    if verbose > 0:
        print("    Characterization environment configured (soma/trunk/step/holding ready).")

def update_param_names(context_dict=None):
    """
    Ensure workers have access to the list of parameter names from the optimization config.
    """
    if not hasattr(context, 'config_file_path') or context.config_file_path is None:
        # Fallback: check if it's in kwargs or use a default if we're in analysis mode
        context.config_file_path = getattr(context, 'config_file_path', 
                                           context().get('config_file_path', 
                                           'config/Nested_optimize_CA1_tbs_e_I_spiking_config.yaml'))
    
    if os.path.isfile(context.config_file_path):
        config = load_yaml(context.config_file_path)
        if 'param_names' in config:
            context.param_names = config['param_names']
            # Ensure it's in the internal dict for context().get()
            context()['param_names'] = config['param_names']

def sync_worker_context(context_dict=None):
    """
    Explicitly synchronize a worker's global context with the provided dictionary.
    """
    if context_dict is not None:
        context.update(context_dict)
        # Compatibility for nested framework which flattens kwargs
        if not hasattr(context, 'kwargs') or context.kwargs is None:
            context.kwargs = context()

def compute_features_input_resistance(x, section_name, model_id=None, export=False, plot=False):
    """
    Inject a hyperpolarizing step current into the specified section and measure local Rin.
    """
    start_time = time.time()
    config_sim_env(context)
    
    update_mechanisms_CA1(x)

    sim = context.sim
    cell = context.cell
    kwargs = context.kwargs
    
    rin_amp_pA = context.input_resistance_test_amp
    rin_start = context.input_resistance_test_start
    rin_dur = context.input_resistance_test_duration
    dt = context.dt
    v_init = context.v_init

    # Determine injection node
    if section_name == 'soma':
        inj_node = cell.tree.root
    else:
        trunk_nodes = cell.get_nodes_of_subtype('trunk')
        if not trunk_nodes:
            return {}
        inj_node = trunk_nodes[0]

    # Shorten integration time to drastically improve speed without CVODE
    sim.tstop = rin_start + rin_dur + 10.0
    h.tstop = sim.tstop

    sim.modify_stim('step', node=inj_node, loc=0.5, amp=rin_amp_pA/1000.0, delay=rin_start, dur=rin_dur)
    sim.run(v_init=v_init)

    tvec = np.array(sim.tvec.to_python())
    results = {}

    if sim.has_rec(section_name):
        v_vec = np.array(sim.get_rec(section_name)['vec'].to_python())
        bl_mask = (tvec >= rin_start - 20.0) & (tvec < rin_start)
        stim_mask = (tvec >= rin_start) & (tvec <= rin_start + rin_dur)
        
        v_rest = np.mean(v_vec[bl_mask])
        v_peak = np.min(v_vec[stim_mask])
        
        rin = (v_peak - v_rest) / (rin_amp_pA / 1000.0)
        results[f'{section_name}_rin'] = float(rin)
        
        if section_name == 'soma':
            results['soma_v_rest'] = float(v_rest)

    sim.modify_stim('step', amp=0.0)

    if (plot or context.kwargs.get('plot', False)) and (model_id == 0 or model_id == '0'):
        plot_sim_traces(sim, ['soma', 'trunk'], f'data/input_resistance_trace_{section_name}_model_{model_id}.png')
        
    return results

def get_objectives_input_resistance(features, targets=None, model_id=None):
    """
    Compute objective residuals for Rin.
    :param features: dict
    :param targets: dict
    :param model_id: int or str
    :return: dict
    """
    if hasattr(context, 'target_val'):
        targets = context.target_val
    else:
        targets = context.kwargs.get('target_val', {})
    objectives = {}
    target_range = context.kwargs.get('target_range', {})
    
    for key in ['soma_rin', 'trunk_rin']:
        if key in features and key in targets:
            norm_range = target_range.get(key, 10.0)
            objectives[key] = ((features[key] - targets[key]) / norm_range)**2
            
    return objectives

def compute_features_rheobase(x, model_id=None, export=False, plot=False):
    """
    Find the minimum somatic current injection (rheobase) required to elicit at least one spike
    using a dynamic incremental step search.
    """
    start_time = time.time()
    config_sim_env(context)
    update_mechanisms_CA1(x)

    sim = context.sim
    cell = context.cell
    
    rheo_start = context.rheobase_test_start
    rheo_dur = context.rheobase_test_duration
    v_init = context.v_init

    amps_to_test = context.rheobase_test_amps
    rheobase = max(amps_to_test) # Default if no spike

    # Shorten integration time to drastically speed up sequential search
    sim.tstop = rheo_start + rheo_dur + 10.0
    h.tstop = sim.tstop

    # Binary search for rheobase to minimize simulation calls
    low = 0
    high = len(amps_to_test) - 1
    rheo_idx = high
    found = False

    while low <= high:
        mid = (low + high) // 2
        amp = amps_to_test[mid]
        
        sim.modify_stim('step', dur=rheo_dur, amp=amp/1000.0, delay=rheo_start)
        # Clear previous spikes specifically
        context.spike_output_vec.resize(0)
        sim.run(v_init=v_init)
        
        spike_times = np.array(context.spike_output_vec.to_python())
        if np.any(spike_times > rheo_start):
            rheo_idx = mid
            found = True
            high = mid - 1 # Try lower
        else:
            low = mid + 1 # Need higher
            
    rheobase = float(amps_to_test[rheo_idx]) if found else float(max(amps_to_test))

    if getattr(context, 'verbose', 0) > 0:
        print(f"    Soma Rheobase: {rheobase:.1f} pA")

    # Re-run at the actual rheobase so the saved plot shows the spike at threshold
    if (plot or context.kwargs.get('plot', False)) and (model_id == 0 or model_id == '0'):
        sim.modify_stim('step', dur=rheo_dur, amp=rheobase/1000.0, delay=rheo_start)
        context.spike_output_vec.resize(0)
        sim.run(v_init=v_init)
        plot_sim_traces(sim, ['soma'], f'data/rheobase_trace_model_{model_id}.png')

    return {'soma_rheobase': float(rheobase)}

def get_objectives_rheobase(features, targets=None, model_id=None):
    if hasattr(context, 'target_val'):
        targets = context.target_val
    else:
        targets = context.kwargs.get('target_val', {})
    objectives = {}
    target_range = context.kwargs.get('target_range', {})
    
    if 'soma_rheobase' in features and 'soma_rheobase' in targets:
        norm_range = target_range.get('soma_rheobase', 10.0)
        objectives['soma_rheobase'] = ((features['soma_rheobase'] - targets['soma_rheobase']) / norm_range)**2
        
    return objectives

def get_args_dynamic_fi(x, features):
    """
    Parallel Worker Mapping: Dynamically pulls the completed Rheobase metric from previously executed Stage 1.
    Distributes discrete tasks across distinct parallel workers, where each worker rapidly evaluates exactly 
    one step amplitude. Thus, a 7-step FI curve natively splits cleanly across 7 isolated computational workers!
    """
    rheobase = features['soma_rheobase']
    if 'i_inj_relative_amp_array' not in context():
        raise Exception('Error: Must initialize context before calling get_args_dynamic_fi')
    relative_amps = list(context.i_inj_relative_amp_array)
    # Pass rheobase, relative amps, and the filtered context dict to ensure workers are initialized
    args = [[rheobase]*len(relative_amps), relative_amps]
    return args

def filter_features_fi_rates(primitives, features, model_id=None, plot=False):
    """
    Collect FI rate results from parallel workers into a single flat dict.
    
    This filter merges all of them into one flat dict so that get_objectives_fi
    can find all fi_rate_rel_N keys in features.

    'primitives' is the list of per-worker result dicts (one per current step).
    'features' is the dict of features already collected from previous stages.
    """
    fi_rates = {'fi_list': []}
    for result in primitives:
        if not result or 'failed' in result:
            continue
        for key, val in result.items():
            if key.startswith('fi_rate_rel_'):
                fi_rates['fi_list'].append(val)
    return fi_rates

def compute_features_fi_step(x, rheobase, relative_amp_nA, model_id=None, export=False, plot=False):
    """
    Evaluate the firing rate for a single current injection amplitude relative to Rheobase.
    """    
    config_sim_env(context)
    update_mechanisms_CA1(x)

    sim = context.sim
    cell = context.cell

    fi_start = context.fi_test_start
    fi_dur = context.fi_test_duration
    v_init = context.v_init

    # amp parameter is relative so we add it to the model's dynamic rheobase
    amp_nA = (rheobase/1000.0) + relative_amp_nA

    # Generate a unique key index so filter_features can collect it robustly
    amp_idx = list(context.i_inj_relative_amp_array).index(relative_amp_nA)

    # Shorten FI integration bounds scaling
    sim.tstop = fi_start + fi_dur + 50.0
    h.tstop = sim.tstop

    if getattr(context, 'verbose', 0) > 0:
        print(f"      [Model {model_id}] Running FI Step +{relative_amp_nA*1000:.1f} pA over Rheobase...")

    sim.modify_stim('step', amp=amp_nA, delay=fi_start, dur=fi_dur)
    sim.run(v_init=v_init)

    spike_times = np.array(context.spike_output_vec.to_python())
    
    # Firing rate in Hz (count / duration in seconds), only count spikes inside step
    indexes = np.where((spike_times > fi_start) & (spike_times < fi_start + fi_dur))[0]
    rate = len(indexes) / (fi_dur / 1000.0)
    
    if getattr(context, 'verbose', 0) > 0:
        print(f"      [Model {model_id}] FI Step completed: {rate:.1f} Hz")
    
    if (plot or context.kwargs.get('plot', False)) and (model_id == 0 or model_id == '0'):
        plot_sim_traces(sim, ['soma'], f'data/fi_trace_rel_{amp_idx}_model_{model_id}.png')

    return {f'fi_rate_rel_{amp_idx}': rate}

def get_objectives_fi(features, targets=None, model_id=None):
    """
    Compute f_I_residuals comparing model rates to the experimental fitted 
    curve (evaluated at relative increments).
    """
    if targets is None:
        if hasattr(context, 'target_val'):
            targets = context.target_val
        else:
            targets = context.kwargs.get('target_val', {})
    objectives = {}
    f_I_residuals = 0.
    exp_rates = context.exp_rate_f_I_array

    for i, this_rate in enumerate(features['fi_list']):
        tol = max(1.0, 0.1 * exp_rates[i])
        f_I_residuals += ((this_rate - exp_rates[i]) / tol) ** 2

    objectives['fi_rate_residuals'] = f_I_residuals
    return objectives


#TODO: more granularity for Calcium channels and kCadepK..
def update_mechanisms_CA1(params):
    """
    Apply optimized biophysical parameters to the CA1 cell.
    Every value passed to modify_mech_param comes from params_dict.
    x0 in the config YAML matches the biophysics YAML defaults exactly.
    Slopes are derived from endpoint params: slope = (val_max - val_min) / (dist_max - dist_min).
    """
    if isinstance(params, (list, np.ndarray)):
        if 'param_names' not in context() or context.param_names is None:
            update_param_names()
        param_names = context().get('param_names')
        if param_names is None:
            raise RuntimeError("update_mechanisms_CA1: param_names missing from context.")
        params_dict = dict(zip(param_names, params))
        params_array = np.array(params)
    else:
        params_dict = params
        param_names = context().get('param_names')
        if param_names is None:
            update_param_names()
            param_names = context().get('param_names')
        params_array = np.array([params_dict[n] for n in param_names]) if param_names else np.array(list(params_dict.values()))

    # Skip if params unchanged
    if hasattr(context, 'last_params_applied') and np.allclose(context.last_params_applied, params_array):
        return
    context.last_params_applied = params_array

    cell = context.cell

    # ------------------------------------------------------------------ #
    # 1.  Na / K  (soma)                                                  #
    # ------------------------------------------------------------------ #
    cell.modify_mech_param('soma', 'nas', 'gbar',    value=params_dict['soma.gbar_nas'])
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', value=params_dict['soma.gkdrbar'])
    cell.modify_mech_param('soma', 'kap', 'gkabar',  value=params_dict['soma.gkabar'])
    cell.modify_mech_param('soma', 'nas', 'sh',      value=params_dict['soma.sh_nas/x'])
    cell.modify_mech_param('soma', 'nas', 'sha',     value=params_dict['soma.sha_nas/x'])

    # ------------------------------------------------------------------ #
    # 2.  Soma Ca channels  (cal, car, cat — single-zone)                #
    # ------------------------------------------------------------------ #
    cell.modify_mech_param('soma', 'cal', 'gcalbar', value=params_dict['soma.cal.gcalbar'])
    cell.modify_mech_param('soma', 'car', 'gcabar',  value=params_dict['soma.car.gcabar'])
    cell.modify_mech_param('soma', 'cat', 'gcatbar', value=params_dict['soma.cat.gcatbar'])

    # ------------------------------------------------------------------ #
    # 3.  calH  (L-type, trunk/apical/tuft — two zones)                  #
    #     zone 0: max_loc=50  (proximal)  <- trunk.calH.gcalbar_prox     #
    #     zone 1: min_loc=50  (distal)    <- apical.Ca.gcalHbar / tuft   #
    # ------------------------------------------------------------------ #
    calH_prox = params_dict['trunk.calH.gcalbar_prox']

    for sec, calH_dist in [('trunk',  params_dict['apical.Ca.gcalHbar']),
                            ('apical', params_dict['apical.Ca.gcalHbar']),
                            ('tuft',   params_dict['tuft.Ca.gcalHbar'])]:
        cell.modify_mech_param(sec, 'calH', 'gcalbar',
                               value=calH_prox, max_loc=50.0, origin='soma')
        cell.modify_mech_param(sec, 'calH', 'gcalbar',
                               value=calH_dist, min_loc=50.0, origin='soma', replace=False)

    # ------------------------------------------------------------------ #
    # 4.  car  (R-type — trunk value; apical/tuft inherit)               #
    # ------------------------------------------------------------------ #
    cell.modify_mech_param('trunk',  'car', 'gcabar', value=params_dict['apical.Ca.gcabar'])
    cell.modify_mech_param('apical', 'car', 'gcabar', origin='trunk')
    cell.modify_mech_param('tuft',   'car', 'gcabar', origin='trunk')

    # ------------------------------------------------------------------ #
    # 5.  cat  (T-type — three zones)                                     #
    #     zone 0: val=0,          max_loc=100                             #
    #     zone 1: slope-ramped,   min_loc=100, max_loc=350                #
    #     zone 2: val=cat_dist,   min_loc=350                             #
    #     slope = (cat_dist - 0) / (350 - 100)  derived from endpoint    #
    # ------------------------------------------------------------------ #
    for sec, cat_dist in [('trunk',  params_dict['apical.Ca.gcatbar']),
                          ('apical', params_dict['apical.Ca.gcatbar']),
                          ('tuft',   params_dict['tuft.Ca.gcatbar'])]:
        cat_slope = cat_dist / 250.0   # (val_350 - val_100=0) / (350 - 100)
        cell.modify_mech_param(sec, 'cat', 'gcatbar',
                               value=0.0, max_loc=100.0, origin='soma')
        cell.modify_mech_param(sec, 'cat', 'gcatbar',
                               value=0.0, min_loc=100.0, max_loc=350.0,
                               origin='soma', slope=cat_slope, replace=False)
        cell.modify_mech_param(sec, 'cat', 'gcatbar',
                               value=cat_dist, min_loc=350.0, origin='soma', replace=False)

    # ------------------------------------------------------------------ #
    # 6.  Cacum tau — soma sets value; dendrites inherit                  #
    # ------------------------------------------------------------------ #
    cell.modify_mech_param('soma', 'Cacum', 'tau', value=params_dict['soma.tau_Cacum_2018'])
    for sec in ['trunk', 'apical', 'tuft', 'basal']:
        cell.modify_mech_param(sec, 'Cacum', 'tau', origin='soma')

    # ------------------------------------------------------------------ #
    # 7.  CadepK — gbkbar and gskbar                                      #
    #     Soma: single zone.  Trunk: three zones.  Apical: inherits trunk #
    #     Tuft: flat value.  All values come from params_dict directly.   #
    # ------------------------------------------------------------------ #
    # Soma
    cell.modify_mech_param('soma', 'CadepK', 'gbkbar', value=params_dict['soma.CadepK.gbkbar'])
    cell.modify_mech_param('soma', 'CadepK', 'gskbar', value=params_dict['soma.CadepK.gskbar'])

    # Trunk — three zones (proximal < 50, mid 50-200, distal > 200)
    cell.modify_mech_param('trunk', 'CadepK', 'gbkbar',
                           value=params_dict['trunk.CadepK.gbkbar_prox'],
                           max_loc=50.0, origin='soma')
    cell.modify_mech_param('trunk', 'CadepK', 'gbkbar',
                           value=params_dict['trunk.CadepK.gbkbar_mid'],
                           min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
    cell.modify_mech_param('trunk', 'CadepK', 'gbkbar',
                           value=params_dict['trunk.CadepK.gbkbar_dist'],
                           min_loc=200.0, origin='soma', replace=False)

    cell.modify_mech_param('trunk', 'CadepK', 'gskbar',
                           value=params_dict['trunk.CadepK.gskbar_prox'],
                           max_loc=50.0, origin='soma')
    cell.modify_mech_param('trunk', 'CadepK', 'gskbar',
                           value=params_dict['trunk.CadepK.gskbar_mid'],
                           min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
    cell.modify_mech_param('trunk', 'CadepK', 'gskbar',
                           value=params_dict['trunk.CadepK.gskbar_dist'],
                           min_loc=200.0, origin='soma', replace=False)

    # Apical inherits trunk three-zone profile
    cell.modify_mech_param('apical', 'CadepK', 'gbkbar', origin='trunk')
    cell.modify_mech_param('apical', 'CadepK', 'gskbar', origin='trunk')

    # Tuft — flat values
    cell.modify_mech_param('tuft', 'CadepK', 'gbkbar', value=params_dict['tuft.CadepK.gbkbar'])
    cell.modify_mech_param('tuft', 'CadepK', 'gskbar', value=params_dict['tuft.CadepK.gskbar'])

    # ------------------------------------------------------------------ #
    # 8.  Reinitialize modified mechanisms                                #
    # ------------------------------------------------------------------ #
    mech_list = ['nas', 'kdr', 'kap', 'cal', 'calH', 'car', 'cat', 'CadepK', 'Cacum']
    sec_list  = ['soma', 'trunk', 'basal', 'apical', 'tuft']
    for mech in mech_list:
        for sec in sec_list:
            if sec in cell.mech_dict and mech in cell.mech_dict[sec]:
                cell.reinitialize_subset_mechanisms(sec, mech)

    if hasattr(cell, 'stim_inh_syns') and cell.stim_inh_syns:
        if 'gabab_gmax' in params_dict:
            for group in cell.stim_inh_syns.values():
                for syn in group:
                    if hasattr(syn, '_syn') and 'GABAb' in syn._syn:
                        syn._syn['GABAb']['target'].gmax = params_dict['gabab_gmax']
        if 'gabaa_gmax' in params_dict:
            for group in cell.stim_inh_syns.values():
                for syn in group:
                    if hasattr(syn, '_syn') and 'GABA_A_KIN' in syn._syn:
                        syn._syn['GABA_A_KIN']['target'].gmax = params_dict['gabaa_gmax']


def evaluate_single_model(parameters, model_id=None):
    """
    Main evaluation loop: Rin -> Rheobase -> FI (Parallel) -> Unitary -> TBS
    Uses the persistent context.cell without rebuilding.
    """
    param_names = context.param_names
    x_dict = {param_names[i]: float(parameters[i]) for i in range(len(param_names))}
    
    # 1. Update Intrinsic Biophysics
    update_mechanisms_CA1(x_dict)

    try:
        # 2. Consolidated Features (Rin, Rheobase, FI)
        # This replaces the previous separate sequential calls for massive speedup
        features = compute_features_intrinsic(x_dict, model_id=model_id)
        if 'failed' in features:
            return None, None

        # 4. Unitary and TBS (evaluate for both genotypes)
        unit_feats = compute_features_unitary_serial(x_dict, model_id=model_id)
        tbs_wt = compute_features_tbs(x_dict, genotype='WT', model_id=model_id)
        tbs_i80t = compute_features_tbs(x_dict, genotype='I80T', model_id=model_id)

        features.update(unit_feats)
        features.update(tbs_wt)
        features.update(tbs_i80t)

        # 5. Objectives
        targets = context.kwargs['target_val']
        objectives = {}
        objectives.update(get_objectives_input_resistance(features, targets))
        objectives.update(get_objectives_fi(features, targets))
        objectives.update(get_objectives_unitary(features, targets))
        objectives.update(get_objectives_tbs(features, targets))
        
        return features, objectives

    except Exception as e:
        print(f"  Evaluation error for Model {model_id}: {e}")
        return None, None

def setup_synapses_for_sim(params, config):
    """
    Assign synapses to the persistent cell exactly once, then update optimized weights.
    """
    if isinstance(params, (list, np.ndarray)):
        params_dict = dict(zip(context.param_names, params))
    else:
        params_dict = params

    cell = context.cell
    if getattr(cell, 'stim_exc_syns', None) is None:
        seed = config.get('synapse_seed', 0)
        local_random = random.Random(seed)
        exc_syn_locs = cell.get_excitatory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
        inh_syn_locs = cell.get_inhibitory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft', 'soma'])

        num_exc = {'CA3':   int(params_dict['num_exc_syns_CA3']),
                   'ECIII': int(params_dict['num_exc_syns_ECIII'])}
        num_inh = {'CA3':   int(params_dict['num_inh_syns_CA3']),
                   'ECIII': int(params_dict['num_inh_syns_ECIII'])}

        excitatory_stochastic = config.get('excitatory_stochastic', True)
        exc, inh = assign_exc_and_inh_synapse_stims(
            cell, exc_syn_locs, inh_syn_locs, 
            ['AMPA_KIN', 'NMDA_KIN5'], ['GABA_A_KIN', 'GABAb'], 
            num_exc, excitatory_stochastic, num_inh, local_random)
        
        cell.stim_exc_syns = exc
        cell.stim_inh_syns = inh
        
        # Ensure trunk and soma have at least one synapse of each type to prevent inheritance crashes 
        # (e.g. tuft inheriting AMPA_KIN from a synapse-less trunk)
        for sec_type in ['soma', 'trunk', 'apical', 'tuft', 'basal']:
            nodes = cell.get_nodes_of_subtype(sec_type)
            if nodes:
                for syn_t in ['AMPA_KIN', 'NMDA_KIN5', 'GABA_A_KIN', 'GABAb']:
                    if not cell.node_has_synapses(nodes[0], syn_t):
                        Synapse(cell, nodes[0], type_list=[syn_t], stochastic=0)

        # Load the base weights from the biophysics YAML
        cell.init_synaptic_mechanisms()

    # Only overwrite parameters actively being optimized
    if 'gabab_gmax' in params_dict:
        gabab_gmax = params_dict['gabab_gmax']
        for group in cell.stim_inh_syns.values():
            for syn in group:
                if hasattr(syn, '_syn') and 'GABAb' in syn._syn:
                    syn._syn['GABAb']['target'].gmax = gabab_gmax
                    
    if 'gabaa_gmax' in params_dict:
        gabaa_gmax = params_dict['gabaa_gmax']
        for group in cell.stim_inh_syns.values():
            for syn in group:
                if hasattr(syn, '_syn') and 'GABA_A_KIN' in syn._syn:
                    syn._syn['GABA_A_KIN']['target'].gmax = gabaa_gmax
                
def _scale_gabab(scale_factor):
    for group in context.cell.stim_inh_syns.values():
        for syn in group:
            if hasattr(syn, '_syn') and 'GABAb' in syn._syn:
                syn._syn['GABAb']['target'].gmax *= scale_factor

def filter_features(primitives, current_features, model_id=None, export=False, plot=False):
    """Safely merge results from both serial (dict) and parallel (list) stages."""
    features = {}
    if isinstance(primitives, dict):
        if not primitives or 'failed' in primitives:
            return {'failed': True}
        features.update(primitives)
    elif isinstance(primitives, list):
        for instance_features in primitives:
            if not instance_features or 'failed' in instance_features:
                return {'failed': True}
            features.update(instance_features)
    return features

    return features

def compute_features_unitary_serial(x, model_id=None, export=False, plot=False):
    """
    Serial evaluation of all Unitary pathways and conditions on a single worker.
    Unitary tests are WT-only. I80T genotype differences are only tested in TBS.
    """
    results = {}
    for pathway in ['CA3', 'ECIII']:
        for condition in ['Control', 'Gabazine']:
            res = compute_features_unitary(x, pathway, condition, model_id=model_id, export=export, plot=plot)
            results.update(res)
    return results

def compute_features_unitary(x, pathway, condition, model_id=None, export=False, plot=False):
    """
    Measure unitary EPSP amplitude and GABAb area for a single pathway/condition.
    Always runs as WT — genotype-specific tests (I80T) happen only in the TBS stage.
    """
    if context_dict is not None:
        if 'param_names' not in context() or not context.param_names:
            context.update(context_dict)
    
    config_sim_env(context)
    update_mechanisms_CA1(x)
    setup_synapses_for_sim(x, context.kwargs)

    sim = context.sim
    cell = context.cell
    config = context.kwargs

    dt = context.dt
    if pathway == 'CA3':
        stim_times = context.CA3_epsp_times
    else:
        stim_times = context.ECIII_epsp_times
    
    # Sim only needs to run 300ms past the last stim (the GABAb analysis window) + small buffer
    sim.parameters['duration'] = max(stim_times) + 310.0
    sim.tstop = sim.parameters['duration']

    # Turn off current injections
    sim.modify_stim('step', amp=0.0)

    # Gabazine removes GABAA
    gabaa_saved = {}
    if condition == 'Gabazine':
        for group in cell.stim_inh_syns.values():
            for syn in group:
                if 'GABA_A_KIN' in syn._syn:
                    t = syn._syn['GABA_A_KIN']['target']
                    gabaa_saved[id(syn)] = t.gmax
                    t.gmax = 0.0

    # Only activate the specific pathway
    play_vec = h.Vector(stim_times)
    active_exc = cell.stim_exc_syns.get(pathway, [])
    active_inh = cell.stim_inh_syns.get(pathway, [])

    for syn in active_exc + active_inh:
        syn.source.play(play_vec)

    if getattr(context, 'verbose', 0) > 0:
        print(f"      [Model {model_id}] Running Unitary test ({pathway} | {condition})...")

    # Use CVODE for unitary subthreshold simulations
    was_cvode = sim.cvode.active() if sim.cvode else False
    if sim.cvode:
        sim.cvode.active(1)

    sim.run(v_init=context.v_init)

    if sim.cvode and not was_cvode:
        sim.cvode.active(0)

    # Stop playing
    for syn in active_exc + active_inh:
        syn.source.play(h.Vector())

    # Restore GABAA
    if condition == 'Gabazine':
        for group in cell.stim_inh_syns.values():
            for syn in group:
                if 'GABA_A_KIN' in syn._syn and id(syn) in gabaa_saved:
                    syn._syn['GABA_A_KIN']['target'].gmax = gabaa_saved[id(syn)]

    # No genotype-specific cleanup needed — unitary is always WT

    soma_v = np.array(sim.get_rec('soma')['vec'].to_python())
    
    # Typically 3 pulses at 300ms ISI. 
    # To capture the full GABAb area, we analyze 300ms of the trace starting from each stimulus
    window_steps = int(300.0 / dt)
    
    segments = []
    for t_start in stim_times:
        idx = int(t_start / dt)
        if idx + window_steps <= len(soma_v):
            # Baseline subtract relative to the exact moment of each stimulus onset
            v_seg = soma_v[idx:idx + window_steps] - soma_v[idx]
            segments.append(v_seg)
            
    avg_v = np.mean(segments, axis=0) if len(segments) > 0 else np.zeros(window_steps)
    
    # For unitary tests, return specific features based on pathway and condition
    prefix = 'ca3' if pathway == 'CA3' else 'ec3'
    cond_prefix = 'gab' if condition == 'Gabazine' else 'ctrl'
    suffix = ''  # Unitary is always WT, no genotype suffix needed
    
    # Peak Amplitude (mirrors find_peaks height limit mathematically assuming no noise)
    peak = float(np.max(avg_v))
    
    res = {}
    res[f'{prefix}_epsp_{cond_prefix}{suffix}'] = peak
    
    if condition == 'Gabazine':
        # Calculate GABAb Area isolating the negative trace using trapz matches analysis_utils
        negative_trace = np.where(avg_v < 0, avg_v, 0)
        area = float(np.trapz(negative_trace, dx=(dt / 1000.0))) # mV * s
        res[f'{prefix}_gab_area{suffix}'] = area

    if (plot or context.kwargs.get('plot', False)) and (model_id == 0 or model_id == '0'):
        plot_sim_traces(sim, ['soma'], f'data/unitary_{prefix}_{cond_prefix}{suffix}_model_{model_id}.png')

    return res

def get_args_static_unitary():
    """
    Returns the pathways and conditions for the unitary synaptic tests.
    """
    pathways = ['CA3', 'CA3', 'ECIII', 'ECIII']
    conditions = ['Control', 'Gabazine', 'Control', 'Gabazine']
    return pathways, conditions

def get_args_static_tbs():
    """
    Parallel Worker Mapping: Executes identical TBS stimulation protocols simultaneously across
    2 isolated processing workers: one for wild type ('WT'), and perfectly identical stimuli into 'I80T'.
    """
    return [['WT', 'I80T']]

def filter_features_tbs(primitives, features, model_id=None, plot=False):
    """
    Merge TBS results from 2 parallel workers (WT and I80T) into a single features dict.

    'primitives' is a LIST of dicts — one per worker:
        primitives[0] = {'plateau_area_cycle_WT': [...], 'spike_count_cycle_WT': [...], 'trough_vm_cycle_WT': [...]}

    We simply merge all keys from each worker result into one flat dict.
    """
    features = {}
    for result in primitives:
        features.update(result)
    return features

def compute_features_tbs(x, genotype, model_id=None, export=False, plot=False):
    pathway = 'CA3' # TBS is restricted exclusively to CA3
    config_sim_env(context)
    update_mechanisms_CA1(x)
    setup_synapses_for_sim(x, context.kwargs)

    sim = context.sim
    cell = context.cell
    config = context.kwargs

    if genotype == 'I80T':
        _scale_gabab(0.5)

    dt = context.dt
    burst_starts = context.burst_starts_ms
    
    # Read from config with fallbacks to context or defaults
    pulses = config.get('pulses_per_burst', getattr(context, 'pulses_per_burst', 5))
    isi = config.get('intra_burst_interval_ms', getattr(context, 'intra_burst_interval_ms', 10.0))
    
    expanded = []
    for s in burst_starts:
        for i in range(pulses):
            expanded.append(s + i * isi)
    play_vec = h.Vector(expanded)

    sim.parameters['duration'] = max(burst_starts) + context.burst_duration_ms + 250.0
    sim.tstop = sim.parameters['duration']

    sim.modify_stim('step', amp=0.0)

    # WT with NO GABAA (Gabazine) is used for TBS in the GNB1 paper
    gabaa_saved = {}
    for group in cell.stim_inh_syns.values():
        for syn in group:
            if 'GABA_A_KIN' in syn._syn:
                t = syn._syn['GABA_A_KIN']['target']
                gabaa_saved[id(syn)] = t.gmax
                t.gmax = 0.0

    active_exc = cell.stim_exc_syns.get(pathway, [])
    active_inh = cell.stim_inh_syns.get(pathway, [])

    for syn in active_exc + active_inh:
        syn.source.play(play_vec)

    if getattr(context, 'verbose', 0) > 0:
        print(f"      [Model {model_id}] Running TBS simulation for {genotype} genotype...")

    sim.run(v_init=context.v_init)

    for syn in active_exc + active_inh:
        syn.source.play(h.Vector())

    # Restore GABAA
    for group in cell.stim_inh_syns.values():
        for syn in group:
            if 'GABA_A_KIN' in syn._syn and id(syn) in gabaa_saved:
                syn._syn['GABA_A_KIN']['target'].gmax = gabaa_saved[id(syn)]

    if genotype == 'I80T':
        _scale_gabab(2.0)

    soma_v = np.array(sim.get_rec('soma')['vec'].to_python())
    bl_end = int(min(burst_starts) / dt)
    v_norm = soma_v - np.mean(soma_v[max(0, bl_end - 4000):bl_end])
    
    areas, spikes, troughs = calculate_cycle_features(v_norm, dt, burst_starts=burst_starts, 
                                    burst_duration=context.burst_duration_ms)
    
    results_dict = {}

    results_dict[f'plateau_area_cycle_{genotype}'] = areas
    results_dict[f'spike_count_cycle_{genotype}'] = spikes
    results_dict[f'trough_vm_cycle_{genotype}'] = troughs

    return results_dict

def get_objectives_unitary(features, targets=None, model_id=None):
    if targets is None:
        if hasattr(context, 'target_val'):
            targets = context.target_val
        else:
            targets = context.kwargs.get('target_val', {})
    objectives = {}
    target_range = context.kwargs.get('target_range', {})
    epsp_range = target_range.get('epsp_amplitude', 0.2)
    gabab_range = target_range.get('gabab_area', 0.05)
    
    for k, v in targets.items():
        if k in features and ('epsp' in k or 'gab_area' in k):
            norm = epsp_range if 'epsp' in k else gabab_range
            objectives[k] = ((features[k] - v) / norm)**2
    return objectives

def get_objectives_tbs(features, targets=None, model_id=None):
    """
    Compute SSE objectives for TBS plateau area, spike count, and trough Vm
    for both WT and I80T conditions.

    features keys (from compute_features_tbs / filter_features_tbs):
        'plateau_area_cycle_WT'  -> list of 5 floats
        'spike_count_cycle_WT'   -> list of 5 floats
        'trough_vm_cycle_WT'     -> list of 5 floats
        'plateau_area_cycle_I80T' -> list of 5 floats  (etc.)

    target_val keys (from YAML):
        plateau_area_wt, spike_count_wt, trough_vm_wt  -> lists of 5
        plateau_area_i80t, spike_count_i80t, trough_vm_i80t -> lists of 5
    """
    if targets is None:
        if hasattr(context, 'target_val'):
            targets = context.target_val
        else:
            targets = context.kwargs.get('target_val', {})

    objectives = {}
    target_range = context.kwargs.get('target_range', {})
    area_range   = target_range['plateau_area']
    spk_range    = target_range['spike_count']
    trough_range = target_range['trough_vm']

    def _sse(model_list, target_list, norm):
        """Normalised sum of squared errors across cycles."""
        sse = 0.0
        for m, t in zip(model_list, target_list):
            sse += ((m - t) / norm) ** 2
        return sse

    # WT
    for feat_key, tgt_key, norm, obj_key in [
        ('plateau_area_cycle_WT',  'plateau_area_wt',  area_range,   'plateau_area_res_WT'),
        ('spike_count_cycle_WT',   'spike_count_wt',   spk_range,    'spike_count_res_WT'),
        ('trough_vm_cycle_WT',     'trough_vm_wt',     trough_range, 'trough_vm_res_WT'),
    ]:
        if feat_key in features and tgt_key in targets:
            objectives[obj_key] = _sse(features[feat_key], targets[tgt_key], norm)

    # I80T
    for feat_key, tgt_key, norm, obj_key in [
        ('plateau_area_cycle_I80T', 'plateau_area_i80t', area_range,   'plateau_area_res_i80t'),
        ('spike_count_cycle_I80T',  'spike_count_i80t',  spk_range,    'spike_count_res_i80t'),
        ('trough_vm_cycle_I80T',    'trough_vm_i80t',    trough_range, 'trough_vm_res_i80t'),
    ]:
        if feat_key in features and tgt_key in targets:
            objectives[obj_key] = _sse(features[feat_key], targets[tgt_key], norm)

    return objectives

def run_sim_tests():
    """
    Standard characterization suite: Rin -> Rheobase -> FI Curve Mapping.
    """
    model_id = context.gid if hasattr(context, 'gid') else 0
    x0 = context.x0_array
    model_label = context.model_key if hasattr(context, 'model_key') else 'test'
    
    print(f"\n--- Starting CA1 Characterization tests for [ {model_label} ] ---")
    
    # 1. Initialize Cell and Environment
    build_sim_env(context)
    config_sim_env(context)
    
    features = {}

    # Stage 0: Input Resistance (Peak)
    print("  Stage 0: Computing Soma & Trunk Input Resistance...")
    for section_name in ['soma', 'trunk']:
        rin_features = compute_features_input_resistance(x0, section_name=section_name, model_id=model_id)
        features.update(rin_features)

    # Stage 1: Rheobase
    print("  Stage 1: Searching for Rheobase...")
    rheo_features = compute_features_rheobase(x0, model_id=model_id)
    features.update(rheo_features)

    # Stage 2: F-I Curve Mapping
    print("  Stage 2: Evaluating F-I Curve (Parallel Mapping)...")
    fi_relative_amps = context.i_inj_relative_amp_array
    group_size = len(fi_relative_amps)
    rheobase = features['soma_rheobase']
    
    # Use context interface to map across workers for FI
    # We pass a filtered context to each worker to ensure they are fully initialized without NEURON objects
    sequences = [[x0] * group_size, [rheobase] * group_size, fi_relative_amps, [get_picklable_context()] * group_size, [model_id] * group_size]
    fi_rates = context.interface.map(compute_features_fi_step, *sequences)
    
    for i, rate_dict in enumerate(fi_rates):
        features.update(rate_dict)
        rate = list(rate_dict.values())[0]
        print(f"    +{fi_relative_amps[i]*1000.0:.0f} pA above Rheobase: {rate:.1f} Hz")

    # Stage 3: Unitary CA3/ECIII
    print("  Stage 3: Unitary Synaptic Responses (WT only)...")
    pathways, conditions = get_args_static_unitary()
    for p, c in zip(pathways, conditions):
        unit = compute_features_unitary(x0, pathway=p, condition=c, model_id=model_id)
        features.update(unit)
        print(f"    {p} {c} EPSP Extracted")
        
    # Stage 4: TBS WT and I80T/+
    print("  Stage 4: Theta Burst Stimulation (TBS) Responses...")
    for geno in ['WT', 'I80T']:
        tbs = compute_features_tbs(x0, genotype=geno, model_id=model_id)
        features.update(tbs)
        print(f"    TBS CA3 {geno} simulated.")

    # Score against targets
    targets = context.kwargs.get('target_val', {})
    objectives = {}
    if targets:
        objectives.update(get_objectives_input_resistance(features, targets))
        objectives.update(get_objectives_fi(features, targets))
        objectives.update(get_objectives_unitary(features, targets))
        objectives.update(get_objectives_tbs(features, targets))
        
        total_err = sum(objectives.values())
        print(f"\n--- Characterization Complete. Total Objective Error: {total_err:.4f} ---\n")
    else:
        print(f"\n--- Characterization Complete. No objective targets found. ---\n")

    if context.export:
        merge_exported_data(context, param_arrays=[x0],
                            model_ids=[model_id], model_labels=[model_label], features=[features],
                            objectives=[objectives], export_file_path=context.export_file_path, 
                            verbose=getattr(context, 'verbose', 0) > 1)
        sys.stdout.flush()

if __name__ == '__main__':
    # When running with mpirun, only Rank 0 should act as the controller.
    # Non-master ranks should wait for the framework to manage them or exit early
    # to avoid conflicting with the nested framework's MPI initialization.
    try:
        from mpi4py import MPI
        if MPI.COMM_WORLD.rank == 0:
            main()
    except ImportError:
        main()


