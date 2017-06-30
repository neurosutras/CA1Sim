__author__ = 'Grace Ng'
import click
from ipyparallel import interactive
from ipyparallel import Client
# from IPython.display import clear_output
from specify_cells3 import *
from moopgen import *
from plot_results import *

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
default_get_features = ('get_stability_features', 'get_adaptation_features')
default_organize_features = ('organize_stability_features', 'organize_adaptation_features')
default_get_objectives = 'get_spiking_objectives'

default_x0_dict = {'soma.gkabar': 2.108E-02, 'soma.gkdrbar': 4.299E-02, 'soma.sh_nas/x': 1.219E+00,
                   'axon.gkbar factor': 1.225E+00, 'dend.gkabar factor': 1.285E-01, 'soma.gCa factor': 3.364E-01,
                   'soma.gCadepK factor': 4.096E+00, 'soma.gkmbar': 4.286E-03} # Err: 4.6192E+05
default_bounds_dict = {'soma.gkabar': (0.01, 0.05), 'soma.gkdrbar': (0.01, 0.06), 'soma.sh_nas/x': (0.1, 6.),
                   'axon.gkbar factor': (1., 3.), 'dend.gkabar factor': (0.1, 5.), 'soma.gCa factor': (0.1, 5.),
                   'soma.gCadepK factor': (0.1, 5.), 'soma.gkmbar': (0.0005, 0.005)}
default_feature_names = ['v_rest', 'v_th', 'soma_peak', 'ADP', 'AHP', 'stability', 'ais_delay', 'slow_depo', 'dend_amp',
                         'spike_count', 'v_min_late', 'rate', 'adi', 'exp_adi', 'f_I']
default_objective_names = ['f_I stability', 'v_th', 'ADP', 'AHP', 'stability', 'slow_depo', 'dend_amp', 'ais_delay',
                           'f_I adaptation', 'adi']
default_target_val = {'v_rest': v_init, 'v_th': -48., 'soma_peak': 40., 'ADP': 0., 'AHP': 4., 'stability': 0.,
                      'ais_delay': 0., 'slow_depo': 20., 'dend_amp': 0.3}
default_target_range = {'v_rest': 0.25, 'v_th': .01, 'soma_peak': 2., 'ADP': 0.01, 'AHP': .005, 'stability': 1.,
                        'ais_delay': 0.0005, 'slow_depo': 0.5, 'dend_amp': 0.0002}
default_optimization_title = 'spiking'
# we should load defaults from a file if we're going to be running optimizations with many more parameters
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
@click.option("--organize-features", nargs=2, type=str, default=())
@click.option("--get-objectives", type=str, default=None)
@click.option("--group-size", type=int, default=2)
@click.option("--pop-size", type=int, default=100)
@click.option("--seed", type=int, default=None)
@click.option("--max-iter", type=int, default=30)
@click.option("--max-gens", type=int, default=None)
@click.option("--path-length", type=int, default=1)
@click.option("--adaptive-step-factor", type=float, default=0.9)
@click.option("--niter-success", type=int, default=None)
@click.option("--survival-rate", type=float, default=0.2)
@click.option("--disp", is_flag=True)
def main(cluster_id, spines, mech_file_path, neurotree_file_path, neurotree_index, param_file_path, param_gen,
         get_features, organize_features, get_objectives, group_size, pop_size, seed, max_iter, max_gens, path_length,
         adaptive_step_factor, niter_success, survival_rate, disp):
    """

    :param cluster_id: str
    :param spines: bool
    :param mech_file_path: str (path)
    :param neurotree_file_path: str (path)
    :param neurotree_index: int
    :param param_file_path: str (path)
    :param param_gen: str (must refer to callable in globals())
    :param get_features: tuple of two str (must refer to callable in globals())
    :param organize_features: tuple of two str (must refer to callable in globals())
    :param get_objectives: str (must refer to callable in globals())
    :param group_size: int
    :param pop_size: int
    :param seed: int
    :param max_iter: int
    :param max_gens: int
    :param path_length: int
    :param adaptive_step_factor: float in [0., 1.]
    :param niter_success: int
    :param survival_rate: float
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

    if param_file_path is None:
        param_file_path = default_param_file_path

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
        param_names = params_dict['x0'].keys()
        x0 = [params_dict['x0'][key] for key in param_names]
        bounds = [params_dict['bounds'][key] for key in param_names]
        feature_names = params_dict['feature_names']
        objective_names = params_dict['objective_names']
        target_val = params_dict['target_val']
        target_range = params_dict['target_range']
        optimization_title = params_dict['optimization_title']
    else:
        param_names = default_x0_dict.keys()
        x0 = [default_x0_dict[key] for key in param_names]
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
    if not organize_features:
        organize_features = default_organize_features
    if organize_features[0] not in globals() or not callable(globals()[organize_features[0]]):
        raise NameError('Multi-Objective Optimization: organize_features: %s has not been imported, or is not a callable '
                        'function.' % organize_features[0])
    if organize_features[1] not in globals() or not callable(globals()[organize_features[1]]):
        raise NameError('Multi-Objective Optimization: organize_features: %s has not been imported, or is not a callable '
                        'function.' % organize_features[1])
    if get_objectives is None:
        get_objectives = default_get_objectives
    if get_objectives not in globals() or not callable(globals()[get_objectives]):
        raise NameError('Multi-Objective Optimization: get_objectives: %s has not been imported, or is not a callable '
                        'function.' % get_objectives)

    globals()['group_size'] = group_size
    globals()['pop_size'] = pop_size

    if group_size > num_procs:
        group_size = num_procs
        print 'Multi-Objective Optimization: group_size adjusted to not exceed num_processes: %i' % num_procs
    un_utilized = num_procs % group_size
    blocks = pop_size / (num_procs / group_size)
    if blocks * (num_procs / group_size) < pop_size:
        blocks += 1

    print 'Multi-Objective Optimization: %s; Total processes: %i; Population size: %i; Group size: %i; ' \
          'Feature calculator: %s; Organize features: %s; Objective calculator: %s; Blocks / generation: %i' % \
          (param_gen, num_procs, pop_size, group_size, get_features, organize_features, get_objectives, blocks)
    if un_utilized > 0:
        print 'Multi-Objective Optimization: %i processes are unutilized' % un_utilized
    sys.stdout.flush()

    param_gen = globals()[param_gen]
    globals()['param_gen'] = param_gen
    get_features = (globals()[get_features[0]], globals()[get_features[1]])
    globals()['get_features'] = get_features
    organize_features = (globals()[organize_features[0]], globals()[organize_features[1]])
    globals()['organize_features'] = organize_features
    get_objectives = globals()[get_objectives]
    globals()['get_objectives'] = get_objectives

    c[:].execute('from parallel_optimize_spiking import *', block=True)
    c[:].map_sync(init_engine, [spines] * num_procs, [mech_file_path] * num_procs, [neurotree_file_path] * num_procs,
                  [neurotree_index] * num_procs, [param_file_path] * num_procs, [disp] * num_procs)

    global local_param_gen
    local_param_gen = param_gen(x0, param_names, feature_names, objective_names, pop_size, bounds=bounds, seed=seed,
                                max_iter=max_iter, max_gens=max_gens, path_length=path_length,
                                adaptive_step_factor=adaptive_step_factor, niter_success=niter_success,
                                survival_rate=survival_rate, disp=disp)

    for generation in local_param_gen():
        features, objectives = compute_features(generation, group_size=group_size, disp=disp)
        local_param_gen.update_population(features, objectives)
    local_param_gen.storage.save(data_dir+history_filename)
    local_param_gen.storage.plot()


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
    if param_file_path is None:
        param_file_path = default_param_file_path

    global x0
    global param_names
    global param_indexes

    if param_file_path is not None:
        params_dict = read_from_pkl(param_file_path)
        param_names = params_dict['x0'].keys()
        x0 = [params_dict['x0'][key] for key in param_names]
    else:
        param_names = default_x0_dict.keys()
        x0 = [default_x0_dict[key] for key in param_names]

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
def compute_features(generation, group_size=1, disp=False, voltage_traces=False):
    """

    :param generation: list of array
    :param group_size: int
    :param disp: bool
    :return: tuple of list of dict
    """
    pop_size = len(generation)
    usable_procs = num_procs - (num_procs % group_size)
    client_ranges = [range(start, start+group_size) for start in range(0, usable_procs, group_size)]
    extra_args = [{} for i in range(pop_size)] #later replaced by updated_extra_args to store rheobase values for
                                                #the second round of getting features
    updated_extra_args = {}
    results = []
    final_results = {}
    for ind in range(len(get_features)):
        feature_function = get_features[ind]
        generation_list = generation
        pop_ids = range(pop_size)
        while len(pop_ids) > 0 or len(results) > 0:
            num_groups = min(len(client_ranges), len(pop_ids))
            if num_groups > 0:
                results.extend(map(feature_function, [generation_list.pop(0) for i in range(num_groups)],
                                   [pop_ids.pop(0) for i in range(num_groups)],
                                   [client_ranges.pop(0) for i in range(num_groups)],
                                   [extra_args.pop(0) for i in range(num_groups)],
                                   [voltage_traces for i in range(num_groups)]))
            if np.any([this_result['async_result'].ready() for this_result in results]) or \
                    np.any([this_result['async_result'] is None for this_result in results]):
                for this_result in results:
                    if this_result['async_result'] is None:
                        client_ranges.append(this_result['client_range'])
                        for feature in feature_names:
                            final_results[this_result['pop_id']][feature] = None
                        updated_extra_args[this_result['pop_id']] = {}
                        results.remove(this_result)
                        if disp:
                            print 'Individual: %i, aborting - cell is spontaneously firing' % \
                                  (this_result['pop_id'])
                    elif this_result['async_result'].ready():
                        client_ranges.append(this_result['client_range'])
                        if disp:
                            flush_engine_buffer(this_result['async_result'])
                        organize_function = organize_features[ind]
                        final_results[this_result['pop_id']] = organize_function(this_result['async_result'].get(),
                                                                        this_result['extra_args'])
                        final_results[this_result['pop_id']].update(this_result['final_results'])
                        updated_extra_args[this_result['pop_id']] = this_result['extra_args']
                        results.remove(this_result)
                        if disp:
                            print 'Individual: %i, computing features took %.2f s' % \
                                  (this_result['pop_id'], this_result['async_result'].wall_time)
                    sys.stdout.flush()
            else:
                time.sleep(1.)
        extra_args = [updated_extra_args[pop_id] for pop_id in range(pop_size)]
    features = [final_results[pop_id] for pop_id in range(pop_size)]
    if voltage_traces is False:
        objectives = map(get_objectives, features)
        return features, objectives
    else:
        return features

@interactive
def get_stability_features(x, pop_id, client_range, extra_args, voltage_traces=False):
    """
    Distribute simulations across available engines for testing spike stability.
    :param x: array
    :param pop_id: int
    :param client_range: list of ints
    :param extra_args: dict
    :return: dict
    """
    results = c[client_range[0]].apply(spike_shape_features, x)
    if results is None:
        print 'Process %i: Aborting - Cell is spontaneously firing.'
        return {'final_results': {}, 'pop_id': pop_id, 'client_range': client_range, 'async_result': None}
    else:
        final_results = results.get()
        rheobase = final_results['amp']
        spike_stability_params = [[rheobase+0.05, 300.], [rheobase+0.5, 100.]]
        dv = c[client_range]
        result = dv.map_async(spike_stability_features, spike_stability_params, [x] * len(spike_stability_params))
        return {'final_results': final_results, 'pop_id': pop_id, 'client_range': client_range, 'async_result': result,
                'extra_args': {'rheobase': rheobase, 'v_th': final_results['v_th']}}

@interactive
def organize_stability_features(async_result, extra_args):
    """

    :param async_result: list of dict (each of two dictionaries contains keys for stability, v_min_late, rate, and amp
    for a particular stimulation)
    :return:
    """
    final_results = {}
    temp_dict = {}
    temp_dict['amp'] = []
    temp_dict['rate'] = []
    for i, this_dict in enumerate(async_result):
        temp_dict['amp'].append(this_dict['amp'])
        temp_dict['rate'].append(this_dict['rate'])
        if 'stability' not in final_results:
            final_results['stability'] = this_dict['stability']
        else:
            final_results['stability'] += this_dict['stability']
        if 'slow_depo' not in final_results:
            final_results['slow_depo'] = this_dict['v_min_late'] - extra_args['v_th']
        else:
            final_results['slow_depo'] += this_dict['v_min_late'] - extra_args['v_th']
    indexes = range(len(temp_dict['rate']))
    indexes.sort(key=temp_dict['amp'].__getitem__)
    temp_dict['amp'] = map(temp_dict['amp'].__getitem__, indexes)
    temp_dict['rate'] = map(temp_dict['rate'].__getitem__, indexes)
    final_results['rate'] = temp_dict['rate'][0]
    final_results['f_I amp'] = temp_dict['amp'][0] #This is used later to calculate the target firing rate for this
                                                  #value of current
    return final_results

@interactive
def get_adaptation_features(x, pop_id, client_range, extra_args, voltage_traces=False):
    """
    :param x: array
    :param pop_id: int
    :param client_range: list of ints
    :param extra_args: dict
    :return: dict
    """
    if not extra_args: #cell was spontaneously firing, so no rheobase was found
        return {'final_results': {}, 'pop_id': pop_id, 'client_range': client_range, 'async_result': None}
    dv = c[client_range]
    rheobase = extra_args['rheobase']

    # Calculate firing rates for a range of I_inj amplitudes using a stim duration of 500 ms
    result = dv.map_async(sim_f_I, [rheobase + i_inj_increment * (i + 1) for i in range(num_increments)])
    if voltage_traces is True:
        global rec_file_list
        dv.execute('export_sim_results()')
        rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir + filename + '.hdf5')]
    return {'final_results': {}, 'pop_id': pop_id, 'client_range': client_range, 'async_result': result,
            'extra_args': {'rheobase': rheobase}}

@interactive
def organize_adaptation_features(async_result, extra_args):
    """

    :param async_result: list of dict (each of two dictionaries contains keys for stability, v_min_late, rate, and amp
    for a particular stimulation)
    :return:
    """
    final_results = {}
    amps = []
    final_results['adi'] = []
    final_results['f_I'] = []
    for i, trial in enumerate(async_result):
        amps.append(trial['amp'])
        spike_times = trial['spike_times']
        if len(spike_times) < 3:
            adi = None
            exp_adi = None
        elif len(spike_times) > len(experimental_spike_times):
            adi = get_adaptation_index(spike_times[:len(experimental_spike_times)])
            exp_adi = experimental_adaptation_indexes[len(experimental_spike_times) - 3]
        else:
            adi = get_adaptation_index(spike_times)
            exp_adi = experimental_adaptation_indexes[len(spike_times) - 3]
        final_results['adi'].append(adi)
        final_results['exp_adi'].append(exp_adi)
        this_rate = len(spike_times) / stim_dur * 1000.
        final_results['f_I'].append(this_rate)
    indexes = range(len(final_results['f_I']))
    indexes.sort(key=amps.__getitem__)
    final_results['adi'] = map(final_results['adi'].__getitem__, indexes)
    final_results['exp_adi'] = map(final_results['exp_adi'].__getitem__, indexes)
    final_results['f_I'] = map(final_results['f_I'].__getitem__, indexes)
    return final_results

@interactive
def get_spiking_objectives(features):
    """

    :param features: dict
    :return: dict
    """
    objectives = {}
    if not features: #Cell was firing spontaenously so no rheobase value
        for objective in objective_names:
            objectives[objective] = 1e9
    else:
        for objective in objective_names:
            objectives[objective] = 0.
        rheobase = features['amp']
        target_f_I_stability = experimental_f_I_slope * np.log(features['f_I amp'] / rheobase)
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
        elif amp >= 0.4: #this was implemented in spike_adaptation_engine; is it also applicable here?
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
    result['v_th'] = threshold
    result['ADP'] = ADP
    result['AHP'] = AHP
    result['amp'] = amp
    result['spike_count'] = len(spike_times)

    return result

@interactive
def spike_stability_features(input_param, x, plot=False):
    """

    :param amp: float
    :param local_x: array
    :param plot: bool
    :return: dict
    """
    amp = input_param[0]
    stim_dur = input_param[1]
    sim.parameters['amp'] = amp
    start_time = time.time()
    update_na_ka_stability(x)
    # sim.cvode_state = True
    soma_vm = offset_vm('soma', v_active)
    # sim.cvode_state = False
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=stim_dur)
    duration = equilibrate + stim_dur + 100.
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    stability = 0.
    result = {}
    sim.modify_stim(0, amp=amp)
    sim.run(v_active)
    if plot:
        sim.plot()
    spike_times = np.subtract(cell.spike_detector.get_recordvec().to_python(), equilibrate)
    rate = len(spike_times) / stim_dur * 1000.
    result['rate'] = rate
    result['amp'] = amp
    vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
    v_rest = np.mean(vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
    v_before = np.max(vm[int((equilibrate - 50.) / dt):int((equilibrate - 1.) / dt)])
    v_after = np.max(vm[-int(50. / dt):-1])
    stability += abs(v_before - v_rest) + abs(v_after - v_rest)
    v_min_late = np.min(vm[int((equilibrate + stim_dur - 20.) / dt):int((equilibrate + stim_dur - 1.) / dt)])
    result['stability'] = stability
    result['v_min_late'] = v_min_late
    print 'Process %i took %.1f s to test spike stability with amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
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
def sim_f_I(amp, x, plot=False):
    """

    :param amp: float
    :param local_x: array
    :param plot: bool
    :return: dict
    """
    update_spike_adaptation(x)
    # sim.cvode_state = True
    soma_vm = offset_vm('soma', v_active)
    # print 'Process %i: Getting here - after offset_vm' % os.getpid()
    # sim.cvode_state = False
    sim.parameters['amp'] = amp
    start_time = time.time()
    sim.modify_stim(0, dur=stim_dur, amp=amp)
    duration = equilibrate + stim_dur
    sim.tstop = duration
    sim.run(v_active)
    if plot:
        sim.plot()
    spike_times = np.subtract(cell.spike_detector.get_recordvec().to_python(), equilibrate)
    result = {}
    result['spike_times'] = spike_times
    result['amp'] = amp
    print 'Process %i took %i s to run simulation with I_inj amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    return result

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

    :param x: array [soma.gkabar, soma.gkdrbar, soma.sh_nas/x, axon.gkbar factor, dend.gkabar factor,
            soma.gCa factor, soma.gCadepK factor, soma.gkmbar]
    """
    if spines is False:
        cell.reinit_mechanisms(reset_cable=True)
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[1])
    cell.modify_mech_param('soma', 'kap', 'gkabar', x[0])
    slope = (x[4] - 1.) * x[0] / 300.
    cell.modify_mech_param('soma', 'nas', 'gbar', soma_na_gbar)
    cell.modify_mech_param('soma', 'nas', 'sh', x[2])
    for sec_type in ['apical']:
        cell.reinitialize_subset_mechanisms(sec_type, 'nas')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=x[0]+slope*75., replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300., value=x[0]+slope*300.,
                               replace=False)
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
        cell.modify_mech_param(sec_type, 'nas', 'sha', 5.)
        cell.modify_mech_param(sec_type, 'nas', 'gbar', soma_na_gbar)
    cell.set_terminal_branch_na_gradient()
    cell.reinitialize_subset_mechanisms('axon_hill', 'kap')
    cell.reinitialize_subset_mechanisms('axon_hill', 'kdr')
    cell.modify_mech_param('ais', 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param('ais', 'kap', 'gkabar', x[0] * x[3])
    cell.modify_mech_param('axon', 'kdr', 'gkdrbar', origin='ais')
    cell.modify_mech_param('axon', 'kap', 'gkabar', origin='ais')
    cell.modify_mech_param('axon_hill', 'nax', 'sh', x[2])
    cell.modify_mech_param('axon_hill', 'nax', 'gbar', soma_na_gbar)
    cell.modify_mech_param('axon', 'nax', 'gbar', soma_na_gbar * 2.)
    for sec_type in ['ais', 'axon']:
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='axon_hill')
    cell.modify_mech_param('soma', 'Ca', 'gcamult', x[5])
    cell.modify_mech_param('soma', 'CadepK', 'gcakmult', x[6])
    cell.modify_mech_param('soma', 'km3', 'gkmbar', x[7])
    cell.modify_mech_param('ais', 'km3', 'gkmbar', x[7] * 3.)
    cell.modify_mech_param('axon_hill', 'km3', 'gkmbar', origin='soma')
    cell.modify_mech_param('axon', 'km3', 'gkmbar', origin='ais')
    if spines is False:
        cell.correct_for_spines()
    cell.set_terminal_branch_na_gradient()

@interactive
def update_spike_adaptation(x):
    """

    :param x: array [soma.gCa factor, soma.gCadepK factor, soma.gkmbar]
    """

    cell.modify_mech_param('soma', 'Ca', 'gcamult', x[0])
    cell.modify_mech_param('soma', 'CadepK', 'gcakmult', x[1])
    cell.modify_mech_param('soma', 'km3', 'gkmbar', x[2])
    for sec_type in ['axon_hill', 'ais', 'axon']:
        cell.reinitialize_subset_mechanisms(sec_type, 'km3')

@interactive
def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
        sim.export_to_file(f)

@interactive
def get_best_voltage_traces(storage, group_size, n=1, discard=True):
    """
    Run simulations on the engines with the last best set of parameters, have the engines export their results to .hdf5,
    and then read in and plot the results.
    :param storage: :class:PopulationStorage
    """
    best_inds = storage.get_best(n)
    if len(best_inds) == 0:
        raise Exception('Storage object is empty')
    for ind in best_inds:
        features = compute_features([ind.x], group_size=group_size, voltage_traces=True)
    plot_best_voltage_traces(discard)

@interactive
def plot_best_voltage_traces(discard=True):
    for rec_filename in rec_file_list:
        with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
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
    if discard:
        for rec_filename in rec_file_list:
            os.remove(data_dir + rec_filename + '.hdf5')

if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1,sys.argv)+1):])

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
num_increments = 8