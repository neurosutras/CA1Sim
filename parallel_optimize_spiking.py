__author__ = 'Grace Ng'
import click
from ipyparallel import interactive
from ipyparallel import Client
# from IPython.display import clear_output
from specify_cells3 import *
from moopgen import *
from plot_results import *

#Is the cell on each engine re-setting properly before the next test?
"""
Aims for spike initiation at initial segment by increasing nax density and decreasing activation V1/2 relative to soma,
axon_hill, and axon compartments. Extend linear kap gradient into basals and obliques, aim for 60% spike attenuation
at bifurcation of trunk and tuft.

Optimizes gbar_nax/nas/sh/sha, gkabar_kap/d, gkdrbar for target na spike threshold, AHP amp, and vm stability

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass


script_filename = 'parallel_optimize_spiking.py'

equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.02
amp = 0.3
th_dvdt = 10.
v_init = -77.
v_active = -77.
i_holding = {'soma': 0., 'dend': 0., 'distal_dend': 0.}
i_th = {'soma': 0.1}
soma_ek = -77.
soma_na_gbar = 0.04

default_mech_file_path = data_dir + '042717 GC optimizing spike stability.pkl'
default_neurotree_file_path = morph_dir + '121516_DGC_trees.pkl'
default_param_gen = 'BGen'
default_get_features = ('get_stability_features', 'get_fI_features')
default_process_features = ('process_stability_features', 'process_fI_features')
default_get_objectives = 'get_spiking_objectives'

default_x0_dict = {'soma.gbar_nas': 0.03, 'dend.gbar_nas': 0.03, 'axon.gbar_nax': 0.06, 'ais.gbar_nax': 1.681E-01,
                   'soma.gkabar': 2.108E-02, 'dend.gkabar': 2.709E-03, 'soma.gkdrbar': 4.299E-02,
                   'axon.gkbar': 5.266E-02, 'soma.sh_nas/x': 1.219E+00, 'ais.sha_nas': -2.659E+00,
                   'soma.gCa factor': 3.364E-01, 'soma.gCadepK factor': 4.096E+00, 'soma.gkmbar': 4.286E-03,
                   'ais.gkmbar': 1.286E-02}
default_param_names = ['soma.gbar_nas', 'dend.gbar_nas', 'axon.gbar_nax', 'ais.gbar_nax', 'soma.gkabar', 'dend.gkabar',
                       'soma.gkdrbar', 'axon.gkbar', 'soma.sh_nas/x', 'ais.sha_nas', 'soma.gCa factor',
                       'soma.gCadepK factor', 'soma.gkmbar', 'ais.gkmbar']
default_bounds_dict = {'soma.gbar_nas': (0.01, 0.05), 'dend.gbar_nas': (0.01, 0.05), 'axon.gbar_nax': (0.02, 0.1),
                       'ais.gbar_nax': (0.02, 0.5), 'soma.gkabar': (0.01, 0.05), 'dend.gkabar': (0.001, 0.25),
                       'soma.gkdrbar': (0.01, 0.06), 'axon.gkbar': (0.01, 0.18), 'soma.sh_nas/x': (0.1, 6.),
                       'ais.sha_nas': (-5., -1.), 'soma.gCa factor': (0.1, 5.), 'soma.gCadepK factor': (0.1, 5.),
                       'soma.gkmbar': (0.0005, 0.005), 'ais.gkmbar': (0.0005, 0.015)}
default_feature_names = ['v_rest', 'v_th', 'ADP', 'AHP', 'stability', 'ais_delay', 'slow_depo', 'dend_amp', 'rheobase',
                         'spike_count', 'rate', 'amp', 'adi', 'exp_adi', 'f_I', 'soma_peak']
default_objective_names = ['f_I stability', 'v_th', 'ADP', 'AHP', 'stability', 'slow_depo', 'dend_amp', 'ais_delay',
                           'f_I adaptation', 'adi']
default_target_val = {'v_rest': v_init, 'v_th': -48., 'soma_peak': 40., 'ADP': 0., 'AHP': 4., 'stability': 0.,
                      'ais_delay': 0., 'slow_depo': 20., 'dend_amp': 0.3}
default_target_range = {'v_rest': 0.25, 'v_th': .01, 'soma_peak': 2., 'ADP': 0.01, 'AHP': .005, 'stability': 1.,
                        'ais_delay': 0.0005, 'slow_depo': 0.5, 'dend_amp': 0.0002}
default_optimization_title = 'spiking'

default_param_file_path = None


@click.command()
@click.option("--cluster-id", type=str, default=None)
@click.option("--spines", is_flag=True)
@click.option("--mech-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None)
@click.option("--neurotree-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None)
@click.option("--neurotree-index", type=int, default=0)
@click.option("--param-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None)
@click.option("--param-gen", type=str, default=None)
@click.option("--get-features", nargs=2, type=str, default=())
@click.option("--process-features", nargs=2, type=str, default=())
@click.option("--get-objectives", type=str, default=None)
@click.option("--group-sizes", nargs=2, type=int, default=(1, 10))
@click.option("--pop-size", type=int, default=100)
@click.option("--seed", type=int, default=None)
@click.option("--max-iter", type=int, default=30)
@click.option("--max-gens", type=int, default=None)
@click.option("--path-length", type=int, default=1)
@click.option("--adaptive-step-factor", type=float, default=0.9)
@click.option("--niter-success", type=int, default=None)
@click.option("--survival-rate", type=float, default=0.2)
@click.option("--optimize", is_flag=True)
@click.option("--storage-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None)
@click.option("--export", type=int, default=None)
@click.option("--export-file-name", type=str, default=None)
@click.option("--disp", is_flag=True)
def main(cluster_id, spines, mech_file_path, neurotree_file_path, neurotree_index, param_file_path, param_gen,
         get_features, process_features, get_objectives, group_sizes, pop_size, seed, max_iter, max_gens, path_length,
         adaptive_step_factor, niter_success, survival_rate, optimize, storage_file_path, export, export_file_name, disp):
    """

    :param cluster_id: str
    :param spines: bool
    :param mech_file_path: str (path)
    :param neurotree_file_path: str (path)
    :param neurotree_index: int
    :param param_file_path: str (path)
    :param param_gen: str (must refer to callable in globals())
    :param get_features: tuple of two str (must refer to callable in globals())
    :param process_features: tuple of two str (must refer to callable in globals())
    :param get_objectives: str (must refer to callable in globals())
    :param group_sizes: tuple of int
    :param pop_size: int
    :param seed: int
    :param max_iter: int
    :param max_gens: int
    :param path_length: int
    :param adaptive_step_factor: float in [0., 1.]
    :param niter_success: int
    :param survival_rate: float
    :param optimize: bool
    :param storage_file_path: str (path)
    :param export: int
    :param export_file_name str
    :param disp: bool
    """
    global c

    if cluster_id is not None:
        c = Client(cluster_id=cluster_id)
    else:
        c = Client()

    global num_procs
    num_procs = len(c)

    if mech_file_path is None:
        mech_file_path = default_mech_file_path

    if neurotree_file_path is None:
        neurotree_file_path = default_neurotree_file_path

    global x0
    global param_names
    global bounds
    global feature_names
    global objective_names
    global target_val
    global target_range
    global optimization_title

    if param_file_path is not None:
        params_dict = read_from_pkl(param_file_path)
        param_names = params_dict['param_names']
        x0 = np.array(make_param_arr(params_dict['x0'], param_names))
        bounds = [params_dict['bounds'][key] for key in param_names]
        feature_names = params_dict['feature_names']
        objective_names = params_dict['objective_names']
        target_val = params_dict['target_val']
        target_range = params_dict['target_range']
        optimization_title = params_dict['optimization_title']
    else:
        param_names = default_param_names
        x0 = np.array(make_param_arr(default_x0_dict, param_names))
        bounds = [default_bounds_dict[key] for key in param_names]
        feature_names = default_feature_names
        objective_names = default_objective_names
        target_val = default_target_val
        target_range = default_target_range
        optimization_title = default_optimization_title

    globals()['path_length'] = path_length

    if param_gen is None:
        param_gen = default_param_gen
    if param_gen not in globals():
        raise NameError('Multi-Objective Optimization: %s has not been imported, or is not a valid class of parameter '
                        'generator.' % param_gen)

    global history_filename
    history_filename = '%s %s %s optimization history' % \
                       (datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)

    if not get_features:
        get_features = default_get_features
    if get_features[0] not in globals() or not callable(globals()[get_features[0]]):
        raise NameError('Multi-Objective Optimization: get_features: %s has not been imported, or is not a callable '
                        'function.' % get_features[0])
    if get_features[1] not in globals() or not callable(globals()[get_features[1]]):
        raise NameError('Multi-Objective Optimization: get_features: %s has not been imported, or is not a callable '
                        'function.' % get_features[1])
    if not process_features:
        process_features = default_process_features
    if process_features[0] not in globals() or not callable(globals()[process_features[0]]):
        raise NameError('Multi-Objective Optimization: process_features: %s has not been imported, or is not a callable '
                        'function.' % process_features[0])
    if process_features[1] not in globals() or not callable(globals()[process_features[1]]):
        raise NameError('Multi-Objective Optimization: process_features: %s has not been imported, or is not a callable '
                        'function.' % process_features[1])
    if get_objectives is None:
        get_objectives = default_get_objectives
    if get_objectives not in globals() or not callable(globals()[get_objectives]):
        raise NameError('Multi-Objective Optimization: get_objectives: %s has not been imported, or is not a callable '
                        'function.' % get_objectives)

    globals()['pop_size'] = pop_size
    new_group_sizes = [group_sizes[i] for i in range(len(group_sizes))]
    blocks = range(len(group_sizes))
    un_utilized = range(len(group_sizes))
    for ind in range(len(group_sizes)):
        if new_group_sizes[ind] > num_procs:
            new_group_sizes[ind] = num_procs
            print 'Multi-Objective Optimization: group_sizes index %i adjusted to not exceed num_processes: %i' % ind, num_procs
        un_utilized[ind] = num_procs % group_sizes[ind]
        blocks[ind] = pop_size / (num_procs / new_group_sizes[ind])
        if blocks[ind] * (num_procs / new_group_sizes[ind]) < pop_size:
            blocks[ind] += 1

    globals()['group_sizes'] = new_group_sizes

    print 'Multi-Objective Optimization: %s; Total processes: %i; Population size: %i; Group sizes: %s; ' \
          'Feature calculator: %s; process features: %s; Objective calculator: %s; Blocks / generation: %s' % \
          (param_gen, num_procs, pop_size, (','.join(str(x) for x in new_group_sizes)), get_features, process_features,
           get_objectives, (','.join(str(x) for x in new_group_sizes)))
    if un_utilized > 0:
        print 'Multi-Objective Optimization: %s processes are unutilized' % (','.join(str(x) for x in un_utilized))
    sys.stdout.flush()

    param_gen = globals()[param_gen]
    globals()['param_gen'] = param_gen
    get_features = (globals()[get_features[0]], globals()[get_features[1]])
    globals()['get_features'] = get_features
    process_features = (globals()[process_features[0]], globals()[process_features[1]])
    globals()['process_features'] = process_features
    get_objectives = globals()[get_objectives]
    globals()['get_objectives'] = get_objectives

    c[:].execute('from parallel_optimize_spiking import *', block=True)
    c[:].map_sync(init_engine, [spines] * num_procs, [mech_file_path] * num_procs, [neurotree_file_path] * num_procs,
                  [neurotree_index] * num_procs, [param_file_path] * num_procs, [disp] * num_procs)

    global local_param_gen

    if optimize:
        local_param_gen = param_gen(param_names, feature_names, objective_names, pop_size, x0=x0, bounds=bounds, seed=seed,
                                    max_iter=max_iter, max_gens=max_gens, path_length=path_length,
                                    adaptive_step_factor=adaptive_step_factor, niter_success=niter_success,
                                    survival_rate=survival_rate, disp=disp)
        run_optimization(new_group_sizes, path_length, disp)
        storage = local_param_gen.storage
    elif storage_file_path is not None:
        storage = PopulationStorage(file_path=storage_file_path)
    if export is not None:
        try:
            storage
        except NameError:
            print 'Storage object has not been defined.'
        best_inds = storage.get_best(n=export, iterations=1)
        best_x_val = [make_param_dict(ind.x, param_names) for ind in best_inds]
        for x_val in best_x_val:
            get_voltage_traces(x_val, group_sizes, combined_rec_filename=export_file_name)


@interactive
def run_optimization(group_sizes, path_length, disp):
    for ind, generation in enumerate(local_param_gen()):
        if (ind > 0) and (ind % path_length == 0):
            local_param_gen.storage.save_gen(data_dir + history_filename, ind - 1, path_length)
        features, objectives = compute_features(generation, group_sizes=group_sizes, disp=disp)
        local_param_gen.update_population(features, objectives)
    local_param_gen.storage.save_gen(data_dir + history_filename, ind, path_length)
    # local_param_gen.storage.plot()

@interactive
def get_adaptation_index(spike_times):
    """
    A large value indicates large degree of spike adaptation (large increases in interspike intervals during a train)
    :param spike_times: list of float
    :return: float
    """
    import numpy as np
    if len(spike_times) < 3:
        return None
    isi = []
    adi = []
    for i in range(len(spike_times) - 1):
        isi.append(spike_times[i + 1] - spike_times[i])
    for i in range(len(isi) - 1):
        adi.append((isi[i + 1] - isi[i]) / (isi[i + 1] + isi[i]))
    return np.mean(adi)

# GC experimental spike adaptation data from Brenner...Aldrich, Nat. Neurosci., 2005
experimental_spike_times = [0., 8.57331572, 21.79656539, 39.24702774, 60.92470277, 83.34214003, 109.5640687,
                            137.1598415, 165.7067371, 199.8546896, 236.2219287, 274.3857332, 314.2404227, 355.2575958,
                            395.8520476, 436.7635403]
experimental_adaptation_indexes = []
for i in range(3, len(experimental_spike_times)+1):
    experimental_adaptation_indexes.append(get_adaptation_index(experimental_spike_times[:i]))
experimental_f_I_slope = 53.  # Hz/ln(pA); rate = slope * ln(current - rheobase)
# GC experimental f-I data from Kowalski J...Pernia-Andrade AJ, Hippocampus, 2016
i_inj_increment = 0.05
num_increments = 10

@interactive
def init_engine(spines=False, mech_file_path=None, neurotree_file_path=None, neurotree_index=0, param_file_path=None,
                disp=False):
    """

    :param spines: bool
    :param mech_file_path: str
    :param neurotree_file_path: str
    :param neurotree_index: int
    :param param_file_path: str (path)
    :param disp: bool
    """
    global x0
    global param_names
    global param_indexes

    if param_file_path is not None:
        params_dict = read_from_pkl(param_file_path)
        param_names = params_dict['param_names']
        x0 = np.array(make_param_arr(params_dict['x0'], param_names))
    else:
        param_names = default_param_names
        x0 = np.array(make_param_arr(default_x0_dict, param_names))

    param_indexes = {param_name: i for i, param_name in enumerate(param_names)}

    if mech_file_path is None:
        mech_file_path = default_mech_file_path

    globals()['spines'] = spines

    if neurotree_file_path is None:
        neurotree_file_path = default_neurotree_file_path
    neurotree_dict = read_from_pkl(neurotree_file_path)[neurotree_index]

    globals()['disp'] = disp

    global cell
    cell = DG_GC(neurotree_dict=neurotree_dict, mech_file_path=mech_file_path, full_spines=spines)
    # in order to use a single engine to compute many different objectives, we're going to need a way to reset the cell
    # to the original state by re-initializing the mech_dict from the mech_file_path

    # get the thickest apical dendrite ~200 um from the soma
    candidate_branches = []
    candidate_diams = []
    candidate_locs = []
    for branch in cell.apical:
        if ((cell.get_distance_to_node(cell.tree.root, branch, 0.) >= 200.) &
                (cell.get_distance_to_node(cell.tree.root, branch, 1.) > 300.) & (not cell.is_terminal(branch))):
            candidate_branches.append(branch)
            for seg in branch.sec:
                loc = seg.x
                if cell.get_distance_to_node(cell.tree.root, branch, loc) > 250.:
                    candidate_diams.append(branch.sec(loc).diam)
                    candidate_locs.append(loc)
                    break
    index = candidate_diams.index(max(candidate_diams))
    dend = candidate_branches[index]
    dend_loc = candidate_locs[index]

    axon_seg_locs = [seg.x for seg in cell.axon[2].sec]

    global rec_locs
    rec_locs = {'soma': 0., 'dend': dend_loc, 'ais': 1., 'axon': axon_seg_locs[0]}
    global rec_nodes
    rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'ais': cell.axon[1], 'axon': cell.axon[2]}
    global rec_filename
    rec_filename = 'sim_output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'_pid'+str(os.getpid())

    global sim
    sim = QuickSim(duration, cvode=False, dt=dt, verbose=False)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    sim.parameters['spines'] = spines

    global spike_output_vec
    spike_output_vec = h.Vector()
    cell.spike_detector.record(spike_output_vec)


@interactive
def compute_features(generation, group_sizes=(1, 10), disp=False, export=False):
    """

    :param generation: list of arr
    :param group_sizes: list of int
    :param disp: bool
    :param export: False (for exporting voltage traces)
    :return: tuple of list of dict
    """
    pop_ids = range(generation)
    results = []
    curr_generation = {pop_id: generation[pop_id] for pop_id in pop_ids}
    final_features = {pop_id: {} for pop_id in pop_ids}
    for ind in range(len(get_features)):
        next_generation = {}
        usable_procs = num_procs - (num_procs % group_sizes[ind])
        client_ranges = [range(start, start + group_sizes[ind]) for start in range(0, usable_procs, group_sizes[ind])]
        feature_function = get_features[ind]
        indivs = [{'pop_id': pop_id, 'x': curr_generation[pop_id],
                   'features': final_features[pop_id]} for pop_id in pop_ids]
        while len(pop_ids) > 0 or len(results) > 0:
            num_groups = min(len(client_ranges), len(pop_ids))
            if num_groups > 0:
                results.extend(map(feature_function, [indivs.pop(0) for i in range(num_groups)],
                                   [client_ranges.pop(0) for i in range(num_groups)], [export] * num_groups))
            if np.any([this_result['async_result'].ready() for this_result in results]):
                for this_result in results:
                    if this_result['async_result'].ready():
                        client_ranges.append(this_result['client_range'])
                        if disp:
                            flush_engine_buffer(this_result['async_result'])
                        get_result = this_result['async_result'].get()
                        if get_result is None:
                            final_features[this_result['pop_id']] = None
                        else:
                            next_generation[this_result['pop_id']] = generation[this_result['pop_id']]
                            new_features = process_features[ind](get_result, final_features[this_result['pop_id']])
                            final_features[this_result['pop_id']].update(new_features)
                            if disp:
                                print 'Individual: %i, computing features took %.2f s' % \
                                      (this_result['pop_id'], this_result['async_result'].wall_time)
                        results.remove(this_result)
                    sys.stdout.flush()
            else:
                time.sleep(1.)
        curr_generation = next_generation
    features = [final_features[pop_id] for pop_id in pop_ids]
    if export is False:
        objectives = map(get_objectives, features)
        return features, objectives
    else:
        return features

@interactive
def get_stability_features(indiv, client_range, export=False):
    """
    Distribute simulations across available engines for testing spike stability.
    :param indiv: dict {'pop_id': pop_id, 'x': x arr, 'features': features dict}
    :param client_range: list of ints
    :param export: False (for exporting voltage traces)
    :return: dict
    """
    x = indiv['x']
    dv = c[client_range]
    result = dv.map_async(spike_shape_features, [x])
    if export is True:
        global rec_file_list
        dv.execute('export_sim_results()')
        rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir + filename + '.hdf5')]
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result}

@interactive
def process_stability_features(get_result, old_features):
    """

    :param get_result: list with one dict
    :return: dict
    """
    return get_result[0]

@interactive
def get_fI_features(indiv, client_range, export=False):
    """
    Distribute simulations across available engines for testing f-I features.
    :param indiv: dict {'pop_id': pop_id, 'x': x arr, 'features': features dict}
    :param client_range: list of ints
    :param export: False (for exporting voltage traces)
    :return: dict
    """
    dv = c[client_range]
    x = indiv['x']
    rheobase = indiv['features']['rheobase']
    # Calculate firing rates for a range of I_inj amplitudes using a stim duration of 500 ms
    result = dv.map_async(sim_f_I, [rheobase + i_inj_increment * (i + 1) for i in range(num_increments)],
                          [x] * num_increments, [False] * (num_increments-1) + [True])
    if export is True:
        global rec_file_list
        dv.execute('export_sim_results()')
        rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir + filename + '.hdf5')]
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result}

@interactive
def process_fI_features(get_result, old_features):
    """

    :param get_result: list of dict (each dict has the results from a particular simulation)
    :param old_features: dict
    :return: dict
    """
    new_features = {}
    temp_dict = {}
    temp_dict['amp'] = []
    temp_dict['rate'] = []
    new_features['adi'] = []
    new_features['exp_adi'] = []
    new_features['f_I'] = []
    for i, this_dict in enumerate(get_result):
        temp_dict['amp'].append(this_dict['amp'])
        temp_dict['rate'].append(this_dict['rate'])
        if 'stability' not in new_features:
            new_features['stability'] = this_dict['stability']
        else:
            new_features['stability'] += this_dict['stability']
        if 'slow_depo' not in new_features:
            new_features['slow_depo'] = this_dict['v_min_late'] - old_features['v_th']
        else:
            new_features['slow_depo'] += this_dict['v_min_late'] - old_features['v_th']

        spike_times = this_dict['spike_times']
        if len(spike_times) < 3:
            adi = None
            exp_adi = None
        elif len(spike_times) > len(experimental_spike_times):
            adi = get_adaptation_index(spike_times[:len(experimental_spike_times)])
            exp_adi = experimental_adaptation_indexes[len(experimental_spike_times) - 3]
        else:
            adi = get_adaptation_index(spike_times)
            exp_adi = experimental_adaptation_indexes[len(spike_times) - 3]
        new_features['adi'].append(adi)
        new_features['exp_adi'].append(exp_adi)
        this_rate = len(spike_times) / stim_dur * 1000.
        new_features['f_I'].append(this_rate)

    stability_ind = range(len(temp_dict['rate']))
    stability_ind.sort(key=temp_dict['amp'].__getitem__)
    temp_dict['amp'] = map(temp_dict['amp'].__getitem__, stability_ind)
    temp_dict['rate'] = map(temp_dict['rate'].__getitem__, stability_ind)
    new_features['rate'] = temp_dict['rate'][0]
    new_features['amp'] = temp_dict['amp'][0] #This is used later to calculate the target firing rate for this
                                                  #value of current

    adapt_ind = range(len(new_features['f_I']))
    adapt_ind.sort(key=temp_dict['amp'].__getitem__)
    new_features['adi'] = map(new_features['adi'].__getitem__, adapt_ind)
    new_features['exp_adi'] = map(new_features['exp_adi'].__getitem__, adapt_ind)
    new_features['f_I'] = map(new_features['f_I'].__getitem__, adapt_ind)
    return new_features

@interactive
def get_spiking_objectives(features):
    """

    :param features: dict
    :return: dict
    """
    objectives = {}
    if features['amp'] is None: #Cell was firing spontaenously so no rheobase value
        for objective in objective_names:
            objectives[objective] = 1e9
    else:
        for objective in objective_names:
            objectives[objective] = 0.
        rheobase = features['amp']
        target_f_I_stability = experimental_f_I_slope * np.log(features['amp'] / rheobase)
        objectives['f_I stability'] = ((features['rate'] - target_f_I_stability) / (0.001 * target_f_I_stability)) ** 2. \
                                      + ((features['spike_count'] - 1.) / 0.002) ** 2.
        for target in ['v_th', 'ADP', 'AHP', 'stability', 'slow_depo', 'dend_amp', 'ais_delay']:
            # don't penalize AHP or slow_depo less than target
            if not ((target == 'AHP' and features[target] < target_val[target]) or
                        (target == 'slow_depo' and features[target] < target_val[target])):
                objectives[target] = ((target_val[target] - features[target]) / target_range[target]) ** 2.
        for i, this_adi in enumerate(features['adi']):
            if this_adi is not None and features['exp_adi'] is not None:
                objectives['adi'] += ((this_adi - features['exp_adi'][i]) / (0.01 * features['exp_adi'])) ** 2.
        target_f_I_adaptation = [experimental_f_I_slope * np.log((rheobase + i_inj_increment * (i + 1)) / rheobase)
                      for i in range(num_increments)]
        for i, this_rate in enumerate(features['f_I']):
            objectives['f_I adaptation'] += ((this_rate - target_f_I_adaptation[i]) / (0.01 *
                                                                                       target_f_I_adaptation[i])) ** 2.

    return objectives

@interactive
def spike_shape_features(x, plot=False):
    """
    :param local_x: array [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor, soma.sh_nas/x,
                    axon.gkbar factor, dend.gkabar factor]
    :param plot: bool
    :return: float
    """
    start_time = time.time()
    cell.reinit_mechanisms(reset_cable=True, from_file=True)
    update_na_ka_stability(x)
    # sim.cvode_state = True
    soma_vm = offset_vm('soma', v_active)
    result = {'v_rest': soma_vm}
    stim_dur = 150.
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=stim_dur)
    duration = equilibrate + stim_dur
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    d_amp = 0.01
    amp = max(0., i_th['soma'] - 0.02)
    while not spike:
        sim.modify_stim(0, amp=amp)
        sim.run(v_active)
        vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
        if np.any(vm[:int(equilibrate/dt)] > -30.):
            print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
            return None
        if np.any(vm[int(equilibrate/dt):int((equilibrate+50.)/dt)] > -30.):
            spike = True
        elif amp >= 0.4:
            print 'Process %i: Aborting - rheobase outside target range' % (os.getpid())
            return None
        else:
            amp += d_amp
            if sim.verbose:
                print 'increasing amp to %.3f' % amp
    sim.parameters['amp'] = amp
    i_th['soma'] = amp
    spike_times = cell.spike_detector.get_recordvec().to_python()
    peak, threshold, ADP, AHP = get_spike_shape(vm, spike_times)
    dend_vm = np.interp(t, sim.tvec, sim.get_rec('dend')['vec'])
    th_x = np.where(vm[int(equilibrate / dt):] >= threshold)[0][0] + int(equilibrate / dt)
    if len(spike_times) > 1:
        end = min(th_x + int(10. / dt), int((spike_times[1] - 5.)/dt))
    else:
        end = th_x + int(10. / dt)
    result['soma_peak'] = peak
    result['v_th'] = threshold
    result['ADP'] = ADP
    result['AHP'] = AHP
    result['rheobase'] = amp
    result['spike_count'] = len(spike_times)
    dend_peak = np.max(dend_vm[th_x:end])
    dend_pre = np.mean(dend_vm[th_x - int(0.2 / dt):th_x - int(0.1 / dt)])
    result['dend_amp'] = (dend_peak - dend_pre) / (peak - threshold)

    # calculate AIS delay
    ais_vm = np.interp(t, sim.tvec, sim.get_rec('ais')['vec'])
    ais_dvdt = np.gradient(ais_vm, dt)
    axon_vm = np.interp(t, sim.tvec, sim.get_rec('axon')['vec'])
    axon_dvdt = np.gradient(axon_vm, dt)
    left = th_x - int(2. / dt)
    right = th_x + int(5. / dt)
    ais_peak = np.max(ais_dvdt[left:right])
    ais_peak_t = np.where(ais_dvdt[left:right] == ais_peak)[0][0] * dt
    axon_peak = np.max(axon_dvdt[left:right])
    axon_peak_t = np.where(axon_dvdt[left:right] == axon_peak)[0][0] * dt
    if axon_peak_t >= ais_peak_t + dt:
        result['ais_delay'] = 0.
    else:
        result['ais_delay'] = ais_peak_t + dt - axon_peak_t
    print 'Process %i took %.1f s to find spike rheobase at amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    if plot:
        sim.plot()
    return result


@interactive
def sim_f_I(amp, x, extend_dur=False, plot=False):
    """

    :param amp: float
    :param local_x: array
    :param plot: bool
    :return: dict
    """
    cell.reinit_mechanisms(reset_cable=True, from_file=True)
    update_na_ka_stability(x)
    # sim.cvode_state = True
    soma_vm = offset_vm('soma', v_active)
    sim.parameters['amp'] = amp
    start_time = time.time()
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=stim_dur, amp=amp)
    if extend_dur:
        duration = equilibrate + stim_dur + 150. #extend duration of simulation to find rebound
    else:
        duration = equilibrate + stim_dur
    sim.tstop = duration
    sim.run(v_active)
    if plot:
        sim.plot()
    spike_times = np.subtract(cell.spike_detector.get_recordvec().to_python(), equilibrate)
    t = np.arange(0., duration, dt)
    stability = 0.
    result = {}
    result['spike_times'] = spike_times
    result['amp'] = amp
    rate = len(spike_times) / stim_dur * 1000.
    result['rate'] = rate
    vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
    v_rest = np.mean(vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
    v_before = np.max(vm[int((equilibrate - 50.) / dt):int((equilibrate - 1.) / dt)])
    v_after = np.max(vm[-int(50. / dt):-1])
    stability += abs(v_before - v_rest) + abs(v_after - v_rest)
    v_min_late = np.min(vm[int((equilibrate + stim_dur - 20.) / dt):int((equilibrate + stim_dur - 1.) / dt)])
    result['stability'] = stability
    result['v_min_late'] = v_min_late
    print 'Process %i took %.1f s to run simulation with I_inj amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    return result

@interactive
def offset_vm(description, vm_target=None):
    """

    :param description: str
    :param vm_target: float
    """
    if vm_target is None:
        vm_target = v_init
    sim.modify_stim(0, amp=0.)
    node = rec_nodes[description]
    loc = rec_locs[description]
    rec_dict = sim.get_rec(description)
    sim.modify_stim(1, node=node, loc=loc, amp=0.)
    rec = rec_dict['vec']
    offset = True
    sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    sim.modify_stim(1, amp=i_holding[description])
    sim.run(vm_target)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < vm_target - 0.5:
        i_holding[description] += 0.01
        while offset:
            if sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(1, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest < vm_target - 0.5:
                i_holding[description] += 0.01
            else:
                offset = False
    elif v_rest > vm_target + 0.5:
        i_holding[description] -= 0.01
        while offset:
            if sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(1, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                i_holding[description] -= 0.01
            else:
                offset = False
    sim.tstop = duration
    return v_rest

@interactive
def get_spike_shape(vm, spike_times):
    """

    :param vm: array
    :return: tuple of float: (v_peak, th_v, ADP, AHP)
    """
    start = int((equilibrate+1.)/dt)
    vm = vm[start:]
    dvdt = np.gradient(vm, dt)
    th_x = np.where(dvdt > th_dvdt)[0]
    if th_x.any():
        th_x = th_x[0] - int(1.6/dt)
    else:
        th_x = np.where(vm > -30.)[0][0] - int(2./dt)
    th_v = vm[th_x]
    v_before = np.mean(vm[th_x-int(0.1/dt):th_x])
    v_peak = np.max(vm[th_x:th_x+int(5./dt)])
    x_peak = np.where(vm[th_x:th_x+int(5./dt)] == v_peak)[0][0]
    if len(spike_times) > 1:
        end = max(th_x+x_peak, int((spike_times[1] - 5.) / dt) - start)
    else:
        end = len(vm)
    v_AHP = np.min(vm[th_x+x_peak:end])
    x_AHP = np.where(vm[th_x+x_peak:end] == v_AHP)[0][0]
    AHP = v_before - v_AHP
    # if spike waveform includes an ADP before an AHP, return the value of the ADP in order to increase error function
    rising_x = np.where(dvdt[th_x+x_peak:th_x+x_peak+x_AHP] > 0.)[0]
    if rising_x.any():
        v_ADP = np.max(vm[th_x+x_peak+rising_x[0]:th_x+x_peak+x_AHP])
        ADP = v_ADP - th_v
    else:
        ADP = 0.
    return v_peak, th_v, ADP, AHP

@interactive
def update_mech_dict(x, update_function):
    update_function(x)
    cell.export_mech_dict(cell.mech_file_path)

@interactive
def update_na_ka_stability(x):
    """

    :param x: array ['soma.gbar_nas', 'dend.gbar_nas', 'axon.gbar_nax', 'ais.gbar_nax', 'soma.gkabar', 'dend.gkabar',
                       'soma.gkdrbar', 'axon.gkbar', 'soma.sh_nas/x', 'ais.sha_nas', 'soma.gCa factor',
                       'soma.gCadepK factor', 'soma.gkmbar', 'ais.gkmbar']
    """
    if spines is False:
        cell.reinit_mechanisms(reset_cable=True)
    cell.modify_mech_param('soma', 'nas', 'gbar', x[param_indexes['soma.gbar_nas']])
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[param_indexes['soma.gkabar']])
    cell.modify_mech_param('soma', 'kap', 'gkabar', x[param_indexes['soma.gkabar']])
    slope = (x[param_indexes['dend.gkabar']] - x[param_indexes['soma.gkabar']]) / 300.
    cell.modify_mech_param('soma', 'nas', 'sh', x[param_indexes['soma.sh_nas/x']])
    for sec_type in ['apical']:
        cell.reinitialize_subset_mechanisms(sec_type, 'nas')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=x[param_indexes['soma.gkabar']]+slope*75., replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300.,
                               value=x[param_indexes['soma.gkabar']]+slope*300., replace=False)
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
        cell.modify_mech_param(sec_type, 'nas', 'sha', 5.)
        cell.modify_mech_param(sec_type, 'nas', 'gbar', x[param_indexes['dend.gbar_nas']])
    cell.set_terminal_branch_na_gradient()
    cell.reinitialize_subset_mechanisms('axon_hill', 'kap')
    cell.reinitialize_subset_mechanisms('axon_hill', 'kdr')
    cell.modify_mech_param('ais', 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param('ais', 'kap', 'gkabar', x[param_indexes['soma.gkabar']] * x[param_indexes['axon.gkbar factor']])
    cell.modify_mech_param('axon', 'kdr', 'gkdrbar', origin='ais')
    cell.modify_mech_param('axon', 'kap', 'gkabar', origin='ais')
    cell.modify_mech_param('axon_hill', 'nax', 'sh', x[param_indexes['axon.gkbar factor']])
    cell.modify_mech_param('axon_hill', 'nax', 'gbar', soma_na_gbar)
    cell.modify_mech_param('axon', 'nax', 'gbar', x[param_indexes['axon.gbar_nax']])
    for sec_type in ['ais', 'axon']:
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='axon_hill')
    cell.modify_mech_param('soma', 'Ca', 'gcamult', x[param_indexes['soma.gCa factor']])
    cell.modify_mech_param('soma', 'CadepK', 'gcakmult', x[param_indexes['soma.gCadepK factor']])
    cell.modify_mech_param('soma', 'km3', 'gkmbar', x[param_indexes['soma.gkmbar']])
    cell.modify_mech_param('ais', 'km3', 'gkmbar', x[param_indexes['ais.gkmbar']])
    cell.modify_mech_param('axon_hill', 'km3', 'gkmbar', origin='soma')
    cell.modify_mech_param('axon', 'km3', 'gkmbar', origin='ais')
    cell.modify_mech_param('ais', 'nax', 'sha', x[param_indexes['ais.sha_nas']])
    cell.modify_mech_param('ais', 'nax', 'gbar', x[param_indexes['ais.gbar_nax']])
    if spines is False:
        cell.correct_for_spines()
    cell.set_terminal_branch_na_gradient()


@interactive
def make_param_dict(x, param_names):
    """

    :param x:
    :param param_names:
    :return:
    """
    return {param_name: x[ind] for ind, param_name in enumerate(param_names)}


@interactive
def make_param_arr(x_dict, param_names):
    """

    :param x_dict:
    :param param_names:
    :return:
    """
    return [x_dict[param_name] for param_name in param_names]

@interactive
def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
        sim.export_to_file(f)

@interactive
def get_voltage_traces(x_val, group_sizes, combined_rec_filename=None, discard=True):
    """
    Run simulations on the engines with the given parameter values, have the engines export their results to .hdf5,
    and then read in and plot the results.

    :param x_val: dict
    :param group_size: int
    """
    x_arr = make_param_arr(x_val, param_names)
    if combined_rec_filename is None:
        combined_rec_filename = 'combined_sim_output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'_pid'+str(os.getpid())
    rec_file_list = []
    rec_file_list.extend(compute_features([x_arr], group_size=group_sizes, export=True))
    combined_rec_filename = combine_output_files(rec_file_list)
    if discard:
        for rec_filename in rec_file_list:
            os.remove(data_dir + rec_filename + '.hdf5')
    # plot_best_voltage_traces(combined_rec_filename)
    return combined_rec_filename

@interactive
def plot_best_voltage_traces(combined_rec_filename):
    with h5py.File(data_dir+combined_rec_filename+'.hdf5', 'r') as f:
        for trial in f.itervalues():
            amplitude = trial.attrs['amp']
            fig, axes = plt.subplots(1)
            for rec in trial['rec'].itervalues():
                axes.plot(trial['time'], rec, label=rec.attrs['description'])
            axes.legend(loc='best', frameon=False, framealpha=0.5)
            axes.set_xlabel('Time (ms)')
            axes.set_ylabel('Vm (mV)')
            axes.set_title('Optimize f_I and spike adaptation: I_inj amplitude %.2f' % amplitude)
            clean_axes(axes)
            fig.tight_layout()
        plt.show()
        plt.close()


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1,sys.argv)+1):])

