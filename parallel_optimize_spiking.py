__author__ = 'Grace Ng'
import click
from ipyparallel import interactive
from ipyparallel import Client
# from IPython.display import clear_output
from specify_cells3 import *
from moopgen import *

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

experimental_f_I_slope = 53.  # Hz/ln(pA); rate = slope * ln(current - rheobase)
# GC experimental f-I data from Kowalski J...Pernia-Andrade AJ, Hippocampus, 2016

default_mech_file_path = data_dir + '042717 GC optimizing spike stability.pkl'
default_neurotree_file_path = morph_dir + '121516_DGC_trees.pkl'
default_param_gen = 'BGen'
default_get_features = 'get_spiking_features'
default_get_objectives = 'get_spiking_objectives'

default_x0_dict = {'soma.gkabar': 2.108E-02, 'soma.gkdrbar': 4.299E-02, 'soma.sh_nas/x': 1.219E+00,
                   'axon.gkbar factor': 1.225E+00, 'dend.gkabar factor': 1.285E-01, 'soma.gCa factor': 3.364E-01,
                   'soma.gCadepK factor': 4.096E+00, 'soma.gkmbar': 4.286E-03} # Err: 4.6192E+05
default_bounds_dict = {'soma.gkabar': (0.01, 0.05), 'soma.gkdrbar': (0.01, 0.06), 'soma.sh_nas/x': (0.1, 6.),
                   'axon.gkbar factor': (1., 3.), 'dend.gkabar factor': (0.1, 5.), 'soma.gCa factor': (0.1, 5.),
                   'soma.gCadepK factor': (0.1, 5.), 'soma.gkmbar': (0.0005, 0.005)}
default_feature_names = ['v_rest', 'v_th', 'soma_peak', 'ADP', 'AHP', 'stability', 'ais_delay', 'slow_depo', 'dend_amp',
                         'spike_count', 'spike_rate', 'v_min_late']
default_objective_names = ['na_ka_stability', 'f_I']
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
@click.option("--get-features", type=str, default=None)
@click.option("--get-objectives", type=str, default=None)
@click.option("--group-size", type=int, default=3)
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
         get_features, get_objectives, group_size, pop_size, seed, max_iter, max_gens, path_length,
         adaptive_step_factor, niter_success, survival_rate, disp):
    """

    :param cluster_id: str
    :param spines: bool
    :param mech_file_path: str (path)
    :param neurotree_file_path: str (path)
    :param neurotree_index: int
    :param param_file_path: str (path)
    :param param_gen: str (must refer to callable in globals())
    :param get_features: str (must refer to callable in globals())
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

    if get_features is None:
        get_features = default_get_features
    if get_features not in globals() or not callable(globals()[get_features]):
        raise NameError('Multi-Objective Optimization: get_features: %s has not been imported, or is not a callable '
                        'function.' % get_features)

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
          'Feature calculator: %s; Objective calculator: %s; Blocks / generation: %i' % \
          (param_gen, num_procs, pop_size, group_size, get_features, get_objectives, blocks)
    if un_utilized > 0:
        print 'Multi-Objective Optimization: %i processes are unutilized' % un_utilized

    param_gen = globals()[param_gen]
    globals()['param_gen'] = param_gen
    get_features = globals()[get_features]
    globals()['get_features'] = get_features
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
def compute_features(generation, group_size=1, disp=False):
    """

    :param generation: list of array
    :param group_size: int
    :param disp: bool
    :return: tuple of list of dict
    """
    pop_size = len(generation)
    pop_ids = range(pop_size)
    client_ranges = [range(start, start+group_size) for start in range(0, num_procs, group_size)]
    results = []
    final_results = {}
    while len(pop_ids) > 0 or len(results) > 0:
        num_groups = min(len(client_ranges), len(pop_ids))
        if num_groups > 0:
            results.extend(map(get_features, [generation.pop(0) for i in range(num_groups)],
                               [pop_ids.pop(0) for i in range(num_groups)],
                               [client_ranges.pop(0) for i in range(num_groups)]))
        if np.any([this_result['async_result'].ready() for this_result in results]):
            for this_result in results:
                if this_result['async_result'].ready():
                    client_ranges.append(this_result['client_range'])
                    if disp:
                        flush_engine_buffer(this_result['async_result'])
                    this_feature_dict = {key: value for result_dict in this_result['async_result'].get()
                                         for key, value in result_dict.iteritems()}
                    final_results[this_result['pop_id']] = this_feature_dict
                    results.remove(this_result)
                    if disp:
                        print 'Individual: %i, computing features took %.2f s' % \
                              (this_result['pop_id'], this_result['async_result'].wall_time)
        else:
            time.sleep(1.)
    features = [final_results[pop_id] for pop_id in range(pop_size)]
    objectives = map(get_objectives, features)
    return features, objectives


@interactive
def get_spiking_features(x, pop_id, client_range):
    """
    Distribute simulations across available engines for optimization of leak conductance density gradient.
    :param x: array
    :return: float
    """
    spike_stability_params = [[0.05, 300.], [0.5, 100.]]
    dv = c[client_range]
    result = dv.map_async(spike_shape_stability, spike_stability_params, [x] * len(spike_stability_params))
    return {'pop_id': pop_id, 'client_range': client_range, 'async_result': result}


@interactive
def get_pas_objectives(features):
    """

    :param features: dict
    :return: dict
    """
    objectives = {}

    return objectives

@interactive
def spike_shape_stability(input_params, x):
    start_time = time.time()
    result = spike_shape_features(x)


@interactive
def spike_shape_features(x):
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
            #Will this be problematic?


        if np.any(vm[int(equilibrate/dt):int((equilibrate+50.)/dt)] > -30.):
            spike = True
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
    result['v_th'] = threshold
    result['ADP'] = ADP
    result['AHP'] = AHP
    result['amp'] = amp
    result['spike_count'] = len(spike_times)
    return result


@interactive
def compute_spike_stability_features(input_param, local_x=None, plot=False):
    """

    :param amp: float
    :param local_x: array
    :param plot: bool
    :return: dict
    """
    amp = input_param[0]
    stim_dur = input_param[1]
    sim.parameters['amp'] = amp
    if local_x is None:
        local_x = x
    start_time = time.time()
    update_na_ka_stability(local_x)
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
def update_mech_dict(x):
    update_na_ka_stability(x)
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
def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
        sim.export_to_file(f)


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1,sys.argv)+1):])