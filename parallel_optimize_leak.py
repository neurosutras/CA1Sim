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

Hierarchical optimization:
I) optimize g_pas for target rinp at soma, trunk bifurcation, and tuft bifurcation [without h].
II) optimize ghbar_h for target rinp at soma, trunk bifurcation, and tuft bifurcation, while also optimizing for v_rest
offset between soma and tuft, and EPSP shape changes between proximal and distal synapses measured at the soma.
III) optimize gbar_nax/nas/sh/sha, gkabar_kap/d, gkdrbar for target na spike threshold, AHP amp, and vm stability

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass


script_filename = 'parallel_optimize_leak.py'

equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.02
# dt = 0.002
amp = 0.3
th_dvdt = 10.
v_init = -77.
v_active = -77.
i_holding = {'soma': 0., 'dend': 0., 'distal_dend': 0.}
soma_ek = -77.

default_mech_file_path = data_dir + '042717 GC optimizing spike stability.pkl'
default_neurotree_file_path = morph_dir + '121516_DGC_trees.pkl'
default_param_gen = 'BGen'
default_get_features = 'get_Rinp_features'
default_get_objectives = 'get_pas_objectives'

default_x0_dict = {'soma.g_pas': 1.050E-10, 'dend.g_pas slope': 1.058E-08, 'dend.g_pas tau': 3.886E+01}  # Error: 4.187E-09
default_bounds_dict = {'soma.g_pas': (1.0E-18, 1.0E-6), 'dend.g_pas slope': (1.0E-12, 1.0E-4),
                  'dend.g_pas tau': (25., 400.)}
default_feature_names = ['soma R_inp', 'dend R_inp', 'distal_dend R_inp']
default_objective_names = ['soma R_inp', 'dend R_inp', 'distal_dend R_inp']
default_target_val = {'soma R_inp': 295., 'dend R_inp': 375.}
default_target_range = {'soma R_inp': 0.5, 'dend R_inp': 1.}
default_optimization_title = 'leak'

# Load defaults from a file containing many parameters. Use write_to_file function in function_lib to generate this pkl file
# Note: click only recognizes this as a path if the string is copied as is into the command line; cannot type a variable
# name that stores this string
default_param_file_path = None #'data/leak_default_param_file.pkl'

@click.command()
@click.option("--cluster-id", type=str, default=None)
@click.option("--spines", is_flag=True)
@click.option("--mech-file-path", type=click.Path(exists=False, file_okay=True, dir_okay=False), default=None)
@click.option("--neurotree-file-path", type=click.Path(exists=False, file_okay=True, dir_okay=False), default=None)
@click.option("--neurotree-index", type=int, default=0)
@click.option("--param-file-path", type=click.Path(exists=False, file_okay=True, dir_okay=False), default=None)
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
    history_filename = '%s %s %s optimization history.hdf5' % \
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
    sys.stdout.flush()

    param_gen = globals()[param_gen]
    globals()['param_gen'] = param_gen
    get_features = globals()[get_features]
    globals()['get_features'] = get_features
    get_objectives = globals()[get_objectives]
    globals()['get_objectives'] = get_objectives

    c[:].execute('from parallel_optimize_leak import *', block=True)
    c[:].map_sync(init_engine, [spines] * num_procs, [mech_file_path] * num_procs, [neurotree_file_path] * num_procs,
                  [neurotree_index] * num_procs, [param_file_path] * num_procs, [disp] * num_procs)
    print 'Initiated engines.'
    sys.stdout.flush()
    global local_param_gen
    local_param_gen = param_gen(x0, param_names, feature_names, objective_names, pop_size, bounds=bounds, seed=seed,
                                max_iter=max_iter, max_gens=max_gens, path_length=path_length,
                                adaptive_step_factor=adaptive_step_factor, niter_success=niter_success,
                                survival_rate=survival_rate, disp=disp)

    for ind, generation in enumerate(local_param_gen()):
        features, objectives = compute_features(generation, group_size=group_size, disp=disp)
        local_param_gen.update_population(features, objectives)
        local_param_gen.storage.save_gen(data_dir+history_filename, ind)
    # local_param_gen.storage.plot()
    # get_best_voltage_traces(local_param_gen.storage, group_size)


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

    # get the most distal terminal branch > 300 um from the soma
    candidate_branches = []
    candidate_end_distances = []
    for branch in (branch for branch in cell.apical if cell.is_terminal(branch)):
        if cell.get_distance_to_node(cell.tree.root, branch, 0.) >= 300.:
            candidate_branches.append(branch)
            candidate_end_distances.append(cell.get_distance_to_node(cell.tree.root, branch, 1.))
    index = candidate_end_distances.index(max(candidate_end_distances))
    distal_dend = candidate_branches[index]
    distal_dend_loc = 1.

    global rec_locs
    rec_locs = {'soma': 0., 'dend': dend_loc, 'distal_dend': distal_dend_loc}
    global rec_nodes
    rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'distal_dend': distal_dend}
    global rec_filename
    rec_filename = 'sim_output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'_pid'+str(os.getpid())

    global sim
    # sim = QuickSim(duration, verbose=False)
    sim = QuickSim(duration, cvode=False, dt=dt, verbose=False)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    sim.parameters['spines'] = spines


@interactive
def compute_features(generation, group_size=1, disp=False, voltage_traces=False):
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
                               [client_ranges.pop(0) for i in range(num_groups)],
                               [voltage_traces for i in range(num_groups)]))
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
    if voltage_traces is False:
        objectives = map(get_objectives, features)
        return features, objectives
    else:
        return features


@interactive
def get_Rinp_features(x, pop_id, client_range, voltage_traces=False):
    """
    Distribute simulations across available engines for optimization of leak conductance density gradient.
    :param x: array (soma.g_pas, dend.g_pas slope, dend.g_pas tau, dend.g_pas xhalf)
    :return: float
    """
    sec_list = ['soma', 'dend', 'distal_dend']
    dv = c[client_range]
    result = dv.map_async(get_Rinp_for_section, sec_list, [x] * len(sec_list))
    if voltage_traces is True:
        global rec_file_list
        dv.execute('export_sim_results()')
        rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir + filename + '.hdf5')]
    return {'pop_id': pop_id, 'client_range': client_range, 'async_result': result}


@interactive
def get_pas_objectives(features):
    """

    :param features: dict
    :return: dict
    """
    objectives = {}
    for feature_name in ['soma R_inp', 'dend R_inp']:
        objective_name = feature_name
        objectives[objective_name] = ((target_val[objective_name] - features[feature_name]) /
                                                  target_range[objective_name]) ** 2.
    this_feature = features['distal_dend R_inp'] - features['dend R_inp']
    objective_name = 'distal_dend R_inp'
    if this_feature < 0.:
        objectives[objective_name] = (this_feature / target_range['dend R_inp']) ** 2.
    else:
        objectives[objective_name] = 0.
    return objectives


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
def update_mech_dict(x):
    update_pas_exp(x)
    cell.export_mech_dict(cell.mech_file_path)


@interactive
def update_pas_exp(x):
    """

    x0 = ['soma.g_pas': 2.28e-05, 'dend.g_pas slope': 1.58e-06, 'dend.g_pas tau': 58.4]
    :param x: array [soma.g_pas, dend.g_pas slope, dend.g_pas tau]
    """
    if spines is False:
        cell.reinit_mechanisms(reset_cable=True)
    cell.modify_mech_param('soma', 'pas', 'g', x[param_indexes['soma.g_pas']])
    cell.modify_mech_param('apical', 'pas', 'g', origin='soma', slope=x[param_indexes['dend.g_pas slope']],
                           tau=x[param_indexes['dend.g_pas tau']])
    for sec_type in ['axon_hill', 'axon', 'ais', 'apical', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')
    if spines is False:
        cell.correct_for_spines()


@interactive
def get_Rinp_for_section(section, x):
    """
    Inject a hyperpolarizing step current into the specified section, and return the steady-state input resistance.
    :param section: str
    :return: dict: {str: float}
    """
    start_time = time.time()
    sim.tstop = duration
    sim.parameters['section'] = section
    sim.parameters['target'] = 'Rinp'
    sim.parameters['optimization'] = 'pas'
    if spines:
        sim_description = 'with_spines'
    else:
        sim_description = 'no_spines'
    sim.parameters['description'] = sim_description
    amp = -0.05
    update_pas_exp(x)
    cell.zero_na()
    offset_vm(section)
    loc = rec_locs[section]
    node = rec_nodes[section]
    rec = sim.get_rec(section)
    sim.modify_stim(0, node=node, loc=loc, amp=amp, dur=stim_dur)
    sim.run(v_init)
    Rinp = get_Rinp(np.array(sim.tvec), np.array(rec['vec']), equilibrate, duration, amp)[2]
    result = {}
    result[section+' R_inp'] = Rinp
    print 'Process:', os.getpid(), 'calculated Rinp for %s in %.1f s, Rinp: %.1f' % (section, time.time() - start_time,
                                                                                    Rinp)
    return result


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
    # plot_best_voltage_traces(discard)

@interactive
def plot_best_voltage_traces(discard=True):
    for rec_filename in rec_file_list:
        with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
            for trial in f.itervalues():
                target = trial.attrs['target']
                section = trial.attrs['section']
                optimization = trial.attrs['optimization']
                fig, axes = plt.subplots(1)
                for rec in trial['rec'].itervalues():
                    axes.plot(trial['time'], rec, label=rec.attrs['description'])
                axes.legend(loc='best', frameon=False, framealpha=0.5)
                axes.set_xlabel('Time (ms)')
                axes.set_ylabel('Vm (mV)')
                axes.set_title('Optimize %s: %s (%s)' % (optimization, target, section))
                clean_axes(axes)
                fig.tight_layout()
                plt.show()
                plt.close()
    if discard:
        for rec_filename in rec_file_list:
            os.remove(data_dir + rec_filename + '.hdf5')


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1,sys.argv)+1):])