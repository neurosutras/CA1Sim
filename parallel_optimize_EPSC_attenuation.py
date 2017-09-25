__author__ = 'Aaron D. Milstein'
from specify_cells4 import *
from plot_results import *
from moopgen import *

"""
Varies properties of synaptic AMPA- and NMDA-type glutamate receptors at excitatory synapses to match target features
of spatiotemporal synaptic integration in dendrites. 

Optimization configurables specified in a YAML file.
Current YAML file_path: data/optimize_synaptic_defaults.yaml
"""

context = Context()


def setup_module_from_file(param_file_path='data/parallel_optimize_EPSC_attenuation_config.yaml', output_dir='data',
                           rec_file_path=None, export_file_path=None, verbose=True, disp=True):
    """
    :param param_file_path: str (.yaml file path)
    :param output_dir: str (dir path)
    :param rec_file_path: str (.hdf5 file path)
    :param export_file_path: str (.hdf5 file path)
    :param verbose: bool
    :param disp: bool
    """
    params_dict = read_from_yaml(param_file_path)
    param_gen = params_dict['param_gen']
    param_names = params_dict['param_names']
    default_params = params_dict['default_params']
    x0 = params_dict['x0']
    for param in default_params:
        params_dict['bounds'][param] = (default_params[param], default_params[param])
    bounds = [params_dict['bounds'][key] for key in param_names]
    if 'rel_bounds' in params_dict:
        rel_bounds = params_dict['rel_bounds']
    else:
        rel_bounds = None
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
        if not callable(func):
            raise Exception('Multi-Objective Optimization: update_params: %s is not a callable function.'
                            % update_params_func_name)
        update_params_funcs.append(func)
    if rec_file_path is None:
        rec_file_path = output_dir + '/sim_output' + datetime.datetime.today().strftime('%m%d%Y%H%M') + \
                        '_pid' + str(os.getpid()) + '.hdf5'
    if export_file_path is None:
        export_file_path = output_dir + '/%s_%s_%s_optimization_exported_traces.hdf5' % \
                           (datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)
    x0_array = param_dict_to_array(x0, param_names)
    context.update(locals())
    config_engine(update_params_funcs, param_names, default_params, rec_file_path, export_file_path, output_dir, disp,
                  **kwargs)
    update_submodule_params(x0_array)


def config_controller(export_file_path, **kwargs):
    """

    :param export_file_path: str
    """
    context.update(locals())
    set_constants()


def config_engine(update_params_funcs, param_names, default_params, rec_file_path, export_file_path, output_dur, disp,
                  mech_file_path, neuroH5_file_path, neuroH5_index, spines, **kwargs):
    """
    :param update_params_funcs: list of function references
    :param param_names: list of str
    :param default_params: dict
    :param rec_file_path: str
    :param export_file_path: str
    :param output_dur: str (dir path)
    :param disp: bool
    :param mech_file_path: str
    :param neuroH5_file_path: str
    :param neuroH5_index: int
    :param spines: bool
    """
    neuroH5_dict = read_from_pkl(neuroH5_file_path)[neuroH5_index]
    param_indexes = {param_name: i for i, param_name in enumerate(param_names)}
    context.update(locals())
    context.update(kwargs)
    set_constants()
    setup_cell(**kwargs)


def set_constants():
    """

    :return:
    """
    seed_offset = 8. * 2e6
    num_branches = 2
    branch_names = ['branch%i' % i for i in xrange(num_branches)]
    syn_conditions = ['control', 'TTX']
    ISI = {'long': 100., 'short': 10.}  # inter-stimulus interval for synaptic stim (ms)
    units_per_sim = 5
    equilibrate = 250.  # time to steady-state
    stim_dur = 150.
    sim_duration = {'long': equilibrate + units_per_sim * ISI['long'] + 50,
                    'short': equilibrate + units_per_sim * ISI['short'] + 150.,
                    'default': equilibrate + stim_dur}
    trace_baseline = 10.
    duration = max(sim_duration.values())
    th_dvdt = 10.
    dt = 0.02
    v_init = -77.
    v_active = -77.
    syn_types = ['EPSC']
    local_random = random.Random()
    i_holding = {'soma': 0.}
    i_th = {'soma': 0.1}

    target_iEPSP_amp = 1.
    context.update(locals())


def setup_cell(verbose=False, cvode=False, daspk=False, **kwargs):
    """

    :param verbose: bool
    :param cvode: bool
    :param daspk: bool
    """
    cell = DG_GC(neuroH5_dict=context.neuroH5_dict, mech_file_path=context.mech_file_path,
                 full_spines=context.spines)
    context.cell = cell
    context.local_random.seed(int(context.seed_offset + context.neuroH5_index))

    # Choose apical branches to measure attenuation of EPSPs generated by dendritic injection of EPSC-shaped currents.
    # Each branch must be > 85 um from the soma.
    candidate_branches = [apical for apical in cell.apical if
                          85. < cell.get_distance_to_node(cell.tree.root, apical) < 135]
    context.local_random.shuffle(candidate_branches)

    syn_list = []
    i = 0
    while len(syn_list) < 2 and i < len(candidate_branches):
        branch = candidate_branches[i]
        parents = [syn.branch.parent for syn in syn_list]
        if branch.parent not in parents:
            syn = Synapse(context.cell, branch, loc=0., syn_types=context.syn_types, stochastic=False)
            syn_list.append(syn)
        i += 1
    if len(syn_list) < context.num_branches:
        print 'parallel_optimize_EPSC_attenuation: cell with index %i: fewer than target number of branches meet ' \
              'criterion.' % context.neuroH5_index
    context.syn_list = syn_list

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
    rec_locs = {'soma': 0., 'dend': dend_loc, 'local_branch': 0.}
    context.rec_locs = rec_locs
    rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'local_branch': dend}
    context.rec_nodes = rec_nodes

    equilibrate = context.equilibrate
    duration = context.duration
    dt = context.dt
    stim_dur = context.stim_dur

    sim = QuickSim(duration, cvode=cvode, daspk=daspk, dt=dt, verbose=verbose)
    sim.parameters['duration'] = duration
    sim.parameters['equilibrate'] = equilibrate
    sim.parameters['spines'] = context.spines
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur, description='step')
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration, description='offset')
    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    context.sim = sim

    context.spike_output_vec = h.Vector()
    cell.spike_detector.record(context.spike_output_vec)


def update_submodule_params(x, local_context=None):
    """

    :param x: array
    :param local_context: :class:'Context'
    """
    if local_context is None:
        local_context = context
    local_context.cell.reinit_mechanisms(from_file=True)
    if not local_context.spines:
        local_context.cell.correct_g_pas_for_spines()
    for update_func in local_context.update_params_funcs:
        update_func(x, local_context)


def get_iEPSP_features_long_ISI(indiv, c=None, client_range=None, export=False):
    """
    Distribute simulations across available engines for measuring iEPSP amplitude.
    :param indiv: dict {'pop_id': pop_id, 'x': x arr, 'features': features dict}
    :param c: :class:'ipyparallel.Client'
    :param client_range: list of int
    :param export: bool
    :return: dict
    """
    if c is not None:
        if client_range is None:
            client_range = range(len(c))
        dv = c[client_range]
        map_func = dv.map_async
    else:
        map_func = map
    x = indiv['x']
    ISI_key = 'long'
    result = map_func(compute_iEPSP_amp_features, [x] * len(context.syn_list), range(len(context.syn_list)),
                      [ISI_key] * len(context.syn_list), [None] * len(context.syn_list),
                      [export] * len(context.syn_list))
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result,
            'filter_features': filter_iEPSP_features}


def get_iEPSP_features_short_ISI(indiv, c=None, client_range=None, export=False):
    """
    Distribute simulations across available engines for measuring iEPSP amplitude.
    :param indiv: dict {'pop_id': pop_id, 'x': x arr, 'features': features dict}
    :param c: :class:'ipyparallel.Client'
    :param client_range: list of int
    :param export: bool
    :return: dict
    """
    if c is not None:
        if client_range is None:
            client_range = range(len(c))
        dv = c[client_range]
        map_func = dv.map_async
    else:
        map_func = map
    x = indiv['x']
    ISI_key = 'short'
    if 'features' in indiv and 'imax' in indiv['features'] and len(indiv['features']['imax']) > 0:
        imax = [indiv['features']['imax'][syn_index] for syn_index in xrange(len(indiv['features']['imax']))]
    else:
        imax = [None for i in xrange(len(context.syn_list))]
    result = map_func(compute_iEPSP_amp_features, [x] * len(context.syn_list), range(len(context.syn_list)),
                      [ISI_key] * len(context.syn_list), imax, [export] * len(context.syn_list))
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result,
            'filter_features': filter_iEPSP_features}


def filter_iEPSP_features(computed_result_list, current_features, target_val, target_range, export=False):
    """

    :param computed_result_list: list of dict (each dict contains results from a single simulation)
    :param current_features: dict
    :param target_val: dict of float
    :param target_range: dict of float
    :param export: bool
    :return: dict
    """
    traces = {}
    iEPSP_amp = {}
    imax = {}
    ISI_key = computed_result_list[0]['traces'].iterkeys().next()
    for this_result_dict in computed_result_list:
        for syn_index in this_result_dict['traces'][ISI_key]:
            traces[syn_index] = this_result_dict['traces'][ISI_key][syn_index]
            iEPSP_amp[syn_index] = {}
            for rec_name in traces[syn_index]:
                iEPSP_amp[syn_index][rec_name] = np.max(traces[syn_index][rec_name])
            imax[syn_index] = this_result_dict['imax'][syn_index]
    new_features = {'iEPSP_amp_' + ISI_key: iEPSP_amp,
                    'iEPSP_traces_' + ISI_key: traces,
                    'imax': imax}
    if export:
        processed_export_file_path = context.export_file_path.replace('.hdf5', '_processed.hdf5')
        with h5py.File(processed_export_file_path, 'a') as f:
            if 'time' not in f:
                f.create_group('time')
            if 'traces' not in f:
                f.create_group('traces')
            if ISI_key not in f['time']:
                t = np.arange(-context.trace_baseline, context.sim_duration[ISI_key] - context.equilibrate, context.dt)
                f['time'].create_dataset(ISI_key, compression='gzip', compression_opts=9, data=t)
            if ISI_key not in f['traces']:
                f['traces'].create_group(ISI_key)
            for rec_name in traces.itervalues().next():
                this_mean_trace = np.mean([traces[syn_index][rec_name] for syn_index in traces], axis=0)
                f['traces'][ISI_key].create_dataset(rec_name, compression='gzip', compression_opts=9,
                                                    data=this_mean_trace)
    return new_features


def get_objectives(features, target_val, target_range):
    """

    :param features: dict
    :param target_val: dict of float
    :param target_range: dict of float
    :return: tuple of dict
    """
    objectives = {}
    objective_names = ['EPSC_attenuation_long_ISI', 'EPSC_attenuation_short_ISI', 'EPSC_amplification_soma',
                       'EPSC_amplification_dend']
    features['EPSC_attenuation_long_ISI'] = np.mean([features['iEPSP_amp_long'][i]['soma'] /
                                                     features['iEPSP_amp_long'][i]['local_branch'] for
                                                     i in xrange(context.num_branches)])
    features['EPSC_attenuation_short_ISI'] = np.mean([features['iEPSP_amp_short'][i]['soma'] /
                                                      features['iEPSP_amp_short'][i]['local_branch'] for
                                                      i in xrange(context.num_branches)])
    features['EPSC_amplification_soma'] = np.mean([features['iEPSP_amp_short'][i]['soma'] /
                                                   features['iEPSP_amp_long'][i]['soma'] for
                                                   i in xrange(context.num_branches)])
    features['EPSC_amplification_dend'] = np.mean([features['iEPSP_amp_short'][i]['local_branch'] /
                                                   features['iEPSP_amp_long'][i]['local_branch'] for
                                                   i in xrange(context.num_branches)])
    for objective_name in objective_names:
        objectives[objective_name] = ((target_val[objective_name] - features[objective_name]) /
                                      target_range[objective_name]) ** 2.
    return features, objectives


def iEPSP_amp_error(x, syn_index):
    """

    :param x: array
    :param syn_index: int
    :return: float
    """
    syn = context.syn_list[syn_index]
    syn.target(context.syn_types[0]).imax = x[0]
    syn.source.play(h.Vector([context.equilibrate]))
    context.sim.tstop = context.equilibrate + context.ISI['long']
    context.sim.run(context.v_init)
    syn.source.play(h.Vector())
    rec = context.sim.get_rec('soma')['vec']
    t = np.arange(0., context.duration, context.dt)
    vm = np.interp(t, context.sim.tvec, rec)
    baseline = np.mean(vm[int((context.equilibrate - 3.) / context.dt):int((context.equilibrate - 1.) / context.dt)])
    vm -= baseline
    iEPSP_amp = np.max(vm[int(context.equilibrate / context.dt):])
    Err = ((iEPSP_amp - context.target_iEPSP_amp) / 0.01) ** 2.
    if context.disp:
        print '%s.imax: %.4f, soma iEPSP amp: %.3f' % (context.syn_types[0], x[0], iEPSP_amp)
    return Err


def compute_iEPSP_amp_features(x, syn_index, ISI_key, imax=None, export=False, plot=False):
    """

    :param x: arr
    :param syn_index: int
    :param ISI_key: str
    :param imax: float
    :param export: bool
    :param plot: bool
    :return: dict
    """
    start_time = time.time()
    update_submodule_params(x, context)

    soma_vm = offset_vm('soma', context.v_init)

    if imax is None:
        imax = 0.05
        result = optimize.minimize(iEPSP_amp_error, [imax], args=(syn_index,), method='L-BFGS-B',
                                   bounds=[(imax / 10., imax * 10.)],
                                   options={'ftol': 1e-3, 'xtol': 1e-3, 'disp': True, 'maxfun': 10})
        imax = result.x[0]

    syn = context.syn_list[syn_index]
    syn.source.play(h.Vector([context.equilibrate + i * context.ISI[ISI_key] for i in range(context.units_per_sim)]))
    syn.target(context.syn_types[0]).imax = imax
    context.sim.parameters['input_loc'] = syn.branch.type
    context.sim.modify_rec(context.sim.get_rec_index('local_branch'), node=syn.branch)
    description = 'compute_iEPSP_amp_features: %s' % (context.branch_names[syn_index])
    duration = context.sim_duration[ISI_key]
    context.sim.tstop = duration
    context.sim.parameters['duration'] = duration
    context.sim.parameters['description'] = 'compute_iEPSP_amp_features'
    context.sim.run(context.v_init)
    dt = context.dt
    equilibrate = context.equilibrate
    interp_t = np.arange(0., duration, dt)
    trace_baseline = context.trace_baseline

    result = {'traces': {ISI_key: {syn_index: {}}},
              'imax': {syn_index: imax}}
    start = int(equilibrate / dt)
    end = int(duration / dt)
    trace_start = start - int(trace_baseline / dt)
    baseline_start, baseline_end = int(start - 3. / dt), int(start - 1. / dt)
    for rec in context.sim.rec_list:
        interp_vm = np.interp(interp_t, context.sim.tvec, rec['vec'])
        baseline = np.mean(interp_vm[baseline_start:baseline_end])
        corrected_vm = interp_vm[trace_start:end] - baseline
        result['traces'][ISI_key][syn_index][rec['description']] = np.array(corrected_vm)
    syn.source.play(h.Vector())

    if context.disp:
        print 'Process: %i: %s took %.3f s' % (os.getpid(), description, time.time() - start_time)
    if plot:
        context.sim.plot()
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
    step_stim_index = context.sim.get_stim_index('step')
    offset_stim_index = context.sim.get_stim_index('offset')
    context.sim.modify_stim(step_stim_index, amp=0.)
    node = context.rec_nodes[description]
    loc = context.rec_locs[description]
    rec_dict = context.sim.get_rec(description)
    context.sim.modify_stim(offset_stim_index, node=node, loc=loc, amp=0.)
    rec = rec_dict['vec']
    offset = True

    equilibrate = context.equilibrate
    dt = context.dt
    duration = context.duration

    context.sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    context.sim.modify_stim(offset_stim_index, amp=context.i_holding[description])
    context.sim.run(vm_target)
    vm = np.interp(t, context.sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < vm_target - 0.5:
        context.i_holding[description] += 0.01
        while offset:
            if context.sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (context.i_holding[description], description)
            context.sim.modify_stim(offset_stim_index, amp=context.i_holding[description])
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
            context.sim.modify_stim(offset_stim_index, amp=context.i_holding[description])
            context.sim.run(vm_target)
            vm = np.interp(t, context.sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                context.i_holding[description] -= 0.01
            else:
                offset = False
    context.sim.tstop = duration
    return v_rest


def get_spike_shape(vm, spike_times):
    """

    :param vm: array
    :param spike_times: array
    :return: tuple of float: (v_peak, th_v, ADP, AHP)
    """
    equilibrate = context.equilibrate
    dt = context.dt
    th_dvdt = context.th_dvdt

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
        end = max(th_x + x_peak + int(2./dt), int((spike_times[1] - 4.) / dt) - start)
    else:
        end = len(vm)
    v_AHP = np.min(vm[th_x+x_peak:end])
    x_AHP = np.where(vm[th_x+x_peak:end] == v_AHP)[0][0]
    AHP = v_before - v_AHP
    # if spike waveform includes an ADP before an AHP, return the value of the ADP in order to increase error function
    ADP = 0.
    rising_x = np.where(dvdt[th_x+x_peak+1:th_x+x_peak+x_AHP-1] > 0.)[0]
    if rising_x.any():
        v_ADP = np.max(vm[th_x+x_peak+1+rising_x[0]:th_x+x_peak+x_AHP])
        pre_ADP = np.mean(vm[th_x+x_peak+1+rising_x[0] - int(0.1/dt):th_x+x_peak+1+rising_x[0]])
        ADP += v_ADP - pre_ADP
    falling_x = np.where(dvdt[th_x + x_peak + x_AHP + 1:end] < 0.)[0]
    if falling_x.any():
        v_ADP = np.max(vm[th_x + x_peak + x_AHP + 1: th_x + x_peak + x_AHP + 1 + falling_x[0]])
        ADP += v_ADP - v_AHP
    return v_peak, th_v, ADP, AHP


def update_nap_params(x, local_context=None):
    """
    :param x: array ['soma.gbar_nas', 'dend.gbar_nas', 'dend.gbar_nas slope', 'dend.gbar_nas min', 'dend.gbar_nas bo',
                    'axon.gbar_nax', 'ais.gbar_nax', 'soma.gkabar', 'dend.gkabar', 'soma.gkdrbar', 'axon.gkabar',
                    'soma.sh_nas/x', 'soma.sha_nas/x', 'ais.sha_nax', 'soma.gCa factor', 'soma.gCadepK factor',
                    'soma.gkmbar', 'ais.gkmbar']
    """
    if local_context is None:
        local_context = context
    cell = local_context.cell
    param_indexes = local_context.param_indexes
    cell.modify_mech_param('soma', 'nas', 'sha', x[param_indexes['soma.sha_nas/x']])
    for sec_type in ['apical']:
        cell.modify_mech_param(sec_type, 'nas', 'sha', origin='soma')
    cell.modify_mech_param('axon_hill', 'nax', 'sha', x[param_indexes['soma.sha_nas/x']])
    cell.modify_mech_param('axon', 'nax', 'sha', origin='axon_hill')
    cell.modify_mech_param('ais', 'nax', 'sha', x[param_indexes['ais.sha_nax']] + x[param_indexes['soma.sha_nas/x']])


def update_spike_shape_params(x, local_context=None):
    """
    :param x: array ['soma.gbar_nas', 'dend.gbar_nas', 'dend.gbar_nas slope', 'dend.gbar_nas min', 'dend.gbar_nas bo',
                    'axon.gbar_nax', 'ais.gbar_nax', 'soma.gkabar', 'dend.gkabar', 'soma.gkdrbar', 'axon.gkabar',
                    'soma.sh_nas/x', 'ais.sha_nax', 'soma.gCa factor', 'soma.gCadepK factor', 'soma.gkmbar',
                    'ais.gkmbar']
    """
    if local_context is None:
        local_context = context
    cell = local_context.cell
    param_indexes = local_context.param_indexes
    cell.modify_mech_param('soma', 'nas', 'gbar', x[param_indexes['soma.gbar_nas']])
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[param_indexes['soma.gkdrbar']])
    cell.modify_mech_param('soma', 'kap', 'gkabar', x[param_indexes['soma.gkabar']])
    slope = (x[param_indexes['dend.gkabar']] - x[param_indexes['soma.gkabar']]) / 300.
    cell.modify_mech_param('soma', 'nas', 'sh', x[param_indexes['soma.sh_nas/x']])
    for sec_type in ['apical']:
        cell.reinitialize_subset_mechanisms(sec_type, 'nas')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=(x[param_indexes['soma.gkabar']] + slope * 75.), replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300.,
                               value=(x[param_indexes['soma.gkabar']] + slope * 300.), replace=False)
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
        cell.modify_mech_param(sec_type, 'nas', 'sha', 0.)  # 5.)
        cell.modify_mech_param(sec_type, 'nas', 'gbar',
                               x[param_indexes['dend.gbar_nas']])
        cell.modify_mech_param(sec_type, 'nas', 'gbar', origin='parent', slope=x[param_indexes['dend.gbar_nas slope']],
                               min=x[param_indexes['dend.gbar_nas min']],
                               custom={'method': 'custom_gradient_by_branch_order',
                                       'branch_order': x[param_indexes['dend.gbar_nas bo']]}, replace=False)
        cell.modify_mech_param(sec_type, 'nas', 'gbar', origin='parent',
                               slope=x[param_indexes['dend.gbar_nas slope']], min=x[param_indexes['dend.gbar_nas min']],
                               custom={'method': 'custom_gradient_by_terminal'}, replace=False)
    cell.reinitialize_subset_mechanisms('axon_hill', 'kap')
    cell.reinitialize_subset_mechanisms('axon_hill', 'kdr')
    cell.modify_mech_param('ais', 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param('ais', 'kap', 'gkabar', x[param_indexes['axon.gkabar']])
    cell.modify_mech_param('axon', 'kdr', 'gkdrbar', origin='ais')
    cell.modify_mech_param('axon', 'kap', 'gkabar', origin='ais')
    cell.modify_mech_param('axon_hill', 'nax', 'sh', x[param_indexes['soma.sh_nas/x']])
    cell.modify_mech_param('axon_hill', 'nax', 'gbar', x[param_indexes['soma.gbar_nas']])
    cell.modify_mech_param('axon', 'nax', 'gbar', x[param_indexes['axon.gbar_nax']])
    for sec_type in ['ais', 'axon']:
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='axon_hill')
    cell.modify_mech_param('soma', 'Ca', 'gcamult', x[param_indexes['soma.gCa factor']])
    cell.modify_mech_param('soma', 'CadepK', 'gcakmult', x[param_indexes['soma.gCadepK factor']])
    cell.modify_mech_param('soma', 'km3', 'gkmbar', x[param_indexes['soma.gkmbar']])
    cell.modify_mech_param('ais', 'km3', 'gkmbar', x[param_indexes['ais.gkmbar']])
    cell.modify_mech_param('axon_hill', 'km3', 'gkmbar', origin='soma')
    cell.modify_mech_param('axon', 'km3', 'gkmbar', origin='ais')
    cell.modify_mech_param('ais', 'nax', 'sha', x[param_indexes['ais.sha_nax']])
    cell.modify_mech_param('ais', 'nax', 'gbar', x[param_indexes['ais.gbar_nax']])


def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(context.rec_file_path, 'a') as f:
        context.sim.export_to_file(f)


def plot_exported_synaptic_features(processed_export_file_path):
    """

    :param processed_export_file_path: str
    :return:
    """
    from matplotlib import cm
    with h5py.File(processed_export_file_path, 'r') as f:
        clustered_branch_names = [key for key in f if key not in ['random', 'time']]
        syn_conditions = f[clustered_branch_names[0]].keys()
        rec_names = [key for key in
                     f[clustered_branch_names[0]][syn_conditions[0]]['compound_EPSP_traces'].itervalues().next()]
        colors = list(cm.Paired(np.linspace(0, 1, len(syn_conditions))))
        for branch in clustered_branch_names:
            for rec in rec_names:
                fig, axes = plt.subplots(2, 2, sharey=True)
                fig.suptitle('Branch: %s, Rec: %s' % (branch, rec))
                ordered_syn_conditions = ['expected', 'control'] + [syn_condition for syn_condition in syn_conditions if
                                                                    syn_condition not in ['expected', 'control']]
                for i, syn_condition in enumerate(ordered_syn_conditions):
                    col = i / 2
                    row = i % 2
                    for num_syns in f[branch][syn_condition]['compound_EPSP_traces']:
                        axes[col][row].plot(f['time']['compound_EPSP'],
                                     f[branch][syn_condition]['compound_EPSP_traces'][num_syns][rec], c='k')
                    axes[col][row].set_xlabel('Time (ms)')
                    axes[col][row].set_title(syn_condition)
                axes[0][0].set_ylabel('Compound EPSP amplitude (mV)')
                clean_axes(axes)
                fig.tight_layout()
                fig.subplots_adjust(top=0.85)

        fig, axes = plt.subplots(1, len(clustered_branch_names), sharey=True)
        if len(clustered_branch_names) == 1:
            axes = [axes]
        fig.suptitle('Rec: soma')
        diagonal = np.linspace(0., np.max(f[branch]['expected']['compound_EPSP_amp']), 10)
        for i, branch in enumerate(clustered_branch_names):
            for j, syn_condition in enumerate((syn_condition for syn_condition in syn_conditions if
                                               syn_condition != 'expected')):
                axes[i].plot(f[branch]['expected']['compound_EPSP_amp'], f[branch][syn_condition]['compound_EPSP_amp'],
                             c=colors[j], label=syn_condition)
                axes[i].set_title('Branch: %s' % branch)
                axes[i].set_xlabel('Expected EPSP amp (mV)')
            axes[i].plot(diagonal, diagonal, c='lightgrey', linestyle='--')
        axes[0].set_ylabel('Actual EPSP amp (mV)')
        axes[0].legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes)
        fig.tight_layout()
        fig.subplots_adjust(top=0.85)

        fig, axes = plt.subplots(1, len(rec_names))
        if len(rec_names) == 1:
            axes = [axes]
        for i, rec in enumerate(rec_names):
            for j, syn_condition in enumerate(f['random']):
                axes[i].plot(f['time']['unitary_EPSP'], f['random'][syn_condition]['mean_unitary_EPSP_traces'][rec],
                             c=colors[j], label=syn_condition)
                axes[i].set_xlabel('Time (ms)')
                axes[i].set_title(rec)
        axes[0].set_ylabel('Unitary EPSP amplitude (mV)')
        axes[0].legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes)
        fig.tight_layout()
    plt.show()
    plt.close()
