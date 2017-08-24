__author__ = 'Grace Ng'
from specify_cells3 import *
from plot_results import *

"""
Aims for spike initiation at initial segment by increasing nax density and decreasing activation V1/2 relative to soma,
axon_hill, and axon compartments. Extend linear kap gradient into basals and obliques, aim for 60% spike attenuation
at bifurcation of trunk and tuft.

Optimize g_pas for target rinp at soma, trunk bifurcation, and tuft bifurcation [without h].

Modify a YAML file to include parameters necessary for this script. Then, set up ipcluster and run parallel_optimize_main.py
Current YAML filepath: data/optimize_leak_defaults.yaml
"""

context = Context()


def setup_module_from_file(param_file_path='data/optimize_leak_defaults.yaml', output_dir='data', rec_file_path=None,
                           export_file_path=None, verbose=True):
    """

    :param param_file_path: str (.yaml file path)
    :param output_dir: str (dir path)
    :param rec_file_path: str (.hdf5 file path)
    :param export_file_path: str (.hdf5 file path)
    :param verbose: bool
    """
    params_dict = read_from_yaml(param_file_path)
    param_gen = params_dict['param_gen']
    param_names = params_dict['param_names']
    default_params = params_dict['default_params']
    x0 = params_dict['x0']
    bounds = [params_dict['bounds'][key] for key in param_names]
    feature_names = params_dict['feature_names']
    objective_names = params_dict['objective_names']
    target_val = params_dict['target_val']
    target_range = params_dict['target_range']
    optimization_title = params_dict['optimization_title']
    kwargs = params_dict['kwargs']  # Extra arguments
    kwargs['verbose'] = verbose
    update_params = params_dict['update_params']
    update_params_funcs = []
    for update_params_func_name in update_params:
        func = globals().get(update_params_func_name)
        if not callable (func):
            raise Exception('Multi-Objective Optimization: update_params: %s is not a callable function.'
                            % (update_params_func_name))
    if rec_file_path is None:
        rec_file_path = 'data/sim_output' + datetime.datetime.today().strftime('%m%d%Y%H%M') + \
                   '_pid' + str(os.getpid()) + '.hdf5'
    if export_file_path is None:
        export_file_path = 'data/%s_%s_%s_optimization_exported_traces.hdf5' % \
                           (datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)
    context.update(locals())
    config_engine(update_params_funcs, param_names, default_params, rec_file_path, export_file_path, **kwargs)


def config_controller(export_file_path, **kwargs):
    """

    :param export_file_path: str (path)
    """
    context.update(locals())
    set_constants()


def config_engine(update_params_funcs, param_names, default_params, rec_file_path, export_file_path, output_dur, disp,
                  mech_file_path, neurotree_file_path, neurotree_index, spines, **kwargs):
    """
    :param update_params_funcs: list of function references
    :param param_names: list of str
    :param default_params: dict
    :param rec_file_path: str
    :param export_file_path: str
    :param output_dur: str (dir path)
    :param disp: bool
    :param mech_file_path: str
    :param neurotree_file_path: str
    :param neurotree_index: int
    :param spines: bool
    """
    neurotree_dict = read_from_pkl(neurotree_file_path)[neurotree_index]
    param_indexes = {param_name: i for i, param_name in enumerate(param_names)}
    context.update(locals())
    context.update(kwargs)
    set_constants()
    setup_cell(**kwargs)


def set_constants():
    equilibrate = 250.  # time to steady-state
    stim_dur = 500.
    duration = equilibrate + stim_dur
    dt = 0.02
    amp = 0.3
    th_dvdt = 10.
    v_init = -77.
    v_active = -77.
    i_holding = {'soma': 0., 'dend': 0., 'distal_dend': 0.}
    soma_ek = -77.
    soma_na_gbar = 0.04
    context.update(locals())


def setup_cell(verbose=False, cvode=False, **kwargs):
    """

    :param verbose: bool
    :param cvode: bool
    """
    cell = DG_GC(neurotree_dict=context.neurotree_dict, mech_file_path=context.mech_file_path, full_spines=context.spines)
    context.cell = cell

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

    rec_locs = {'soma': 0., 'dend': dend_loc, 'distal_dend': distal_dend_loc}
    context.rec_locs = rec_locs
    rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'distal_dend': distal_dend}
    context.rec_nodes = rec_nodes

    equilibrate = context.equilibrate
    stim_dur = context.stim_dur
    duration = context.duration
    dt = context.dt

    sim = QuickSim(duration, cvode=cvode, dt=dt, verbose=verbose)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)
    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    sim.parameters['spines'] = context.spines
    context.sim = sim

    context.spike_output_vec = h.Vector()
    cell.spike_detector.record(context.spike_output_vec)


def update_mech_dict(x, update_function, mech_file_path):
    for update_func in context.update_params_funcs:
        update_func(x, context)
    context.cell.export_mech_dict(mech_file_path)


def get_Rinp_features(indiv, c, client_range, export=False):
    """
    Distribute simulations across available engines for testing spike stability.
    :param indiv: dict {'pop_id': pop_id, 'x': x arr, 'features': features dict}
    :param c: Client object
    :param client_range: list of ints
    :param export: False (for exporting voltage traces)
    :return: dict
    """
    dv = c[client_range]
    x = indiv['x']
    sec_list = ['soma', 'dend', 'distal_dend']
    result = dv.map_async(compute_Rinp_features, sec_list, [x] * len(sec_list), [export] * len(sec_list))
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result}


def get_objectives(features, objective_names, target_val, target_range):
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
    return features, objectives


def compute_Rinp_features(section, x, export=False):
    """
    Inject a hyperpolarizing step current into the specified section, and return the steady-state input resistance.
    :param section: str
    :param x: array
    :param export: bool
    :return: dict: {str: float}
    """
    start_time = time.time()
    context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    for update_func in context.update_params_funcs:
        update_func(x, context)
    context.cell.zero_na()

    duration = context.duration
    stim_dur = context.stim_dur
    equilibrate = context.equilibrate
    v_init = context.v_init
    context.sim.tstop = duration
    context.sim.parameters['section'] = section
    context.sim.parameters['target'] = 'Rinp'
    context.sim.parameters['optimization'] = 'pas'
    context.sim.parameters['description'] = 'Rinp ' + section
    amp = -0.05
    context.sim.parameters['amp'] = amp
    offset_vm(section)
    loc = context.rec_locs[section]
    node = context.rec_nodes[section]
    rec = context.sim.get_rec(section)
    context.sim.modify_stim(0, node=node, loc=loc, amp=amp, dur=stim_dur)
    context.sim.run(v_init)
    Rinp = get_Rinp(np.array(context.sim.tvec), np.array(rec['vec']), equilibrate, duration, amp)[2]
    result = {}
    result[section+' R_inp'] = Rinp
    print 'Process:', os.getpid(), 'calculated Rinp for %s in %.1f s, Rinp: %.1f' % (section, time.time() - start_time,
                                                                                    Rinp)
    if export:
        export_sim_results()
    return result


def offset_vm(description, vm_target=None):
    """

    :param description: str
    :param vm_target: float
    """
    if vm_target is None:
        vm_target = context.v_init
    context.sim.modify_stim(0, amp=0.)
    node = context.rec_nodes[description]
    loc = context.rec_locs[description]
    rec_dict = context.sim.get_rec(description)
    context.sim.modify_stim(1, node=node, loc=loc, amp=0.)
    rec = rec_dict['vec']
    offset = True

    equilibrate = context.equilibrate
    dt = context.dt
    duration = context.duration

    context.sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    context.sim.modify_stim(1, amp=context.i_holding[description])
    context.sim.run(vm_target)
    vm = np.interp(t, context.sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < vm_target - 0.5:
        context.i_holding[description] += 0.01
        while offset:
            if context.sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (context.i_holding[description], description)
            context.sim.modify_stim(1, amp=context.i_holding[description])
            context.sim.run(vm_target)
            vm = np.interp(t, context.sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest < vm_target - 0.5:
                context.i_holding[description] += 0.01
            else:
                offset = False
    elif v_rest > vm_target + 0.5:
        context.i_holding[description] -= 0.01
        while offset:
            if context.sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (context.i_holding[description], description)
            context.sim.modify_stim(1, amp=context.i_holding[description])
            context.sim.run(vm_target)
            vm = np.interp(t, context.sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                context.i_holding[description] -= 0.01
            else:
                offset = False
    context.sim.tstop = duration
    return v_rest


def update_pas_exp(x, local_context=None):
    """

    x0 = ['soma.g_pas': 2.28e-05, 'dend.g_pas slope': 1.58e-06, 'dend.g_pas tau': 58.4]
    :param x: array [soma.g_pas, dend.g_pas slope, dend.g_pas tau]
    """
    if local_context is None:
        local_context = context
    cell = local_context.cell
    param_indexes = local_context.param_indexes
    default_params = local_context.default_params
    if local_context.spines is False:
        cell.reinit_mechanisms(reset_cable=True)
    cell.modify_mech_param('soma', 'pas', 'g', find_param_value('soma.g_pas', x, param_indexes, default_params))
    cell.modify_mech_param('apical', 'pas', 'g', origin='soma', slope=find_param_value('dend.g_pas slope', x,
                                                                                       param_indexes, default_params),
                           tau=find_param_value('dend.g_pas tau', x, param_indexes, default_params))
    for sec_type in ['axon_hill', 'axon', 'ais', 'apical', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')
    if local_context.spines is False:
        cell.correct_for_spines()


def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(context.rec_file_path, 'a') as f:
        context.sim.export_to_file(f)