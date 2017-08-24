__author__ = 'Grace Ng'
from specify_cells3 import *
from plot_results import *

"""
Submodule used by parallel_optimize to tune spike shape, f-I curve, and spike adaptation.
Requires a YAML file to specify required configuration parameters. 
Requires use of an ipyparallel client.
"""

context = Context()


def setup_module_from_file(param_file_path='data/optimize_spiking_defaults.yaml', output_dir='data', rec_file_path=None,
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
        if not callable(func):
            raise Exception('Multi-Objective Optimization: update_params: %s is not a callable function.'
                            % (update_params_func_name))
    if rec_file_path is None:
        rec_file_path = output_dir + '/sim_output' + datetime.datetime.today().strftime('%m%d%Y%H%M') + \
                   '_pid' + str(os.getpid()) + '.hdf5'
    if export_file_path is None:
        export_file_path = output_dir + '%s_%s_%s_optimization_exported_traces.hdf5' % \
                           (datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)
    context.update(locals())
    config_engine(update_params_funcs, param_names, default_params, rec_file_path, export_file_path, **kwargs)


def config_controller(export_file_path, **kwargs):
    """

    :param export_file_path: str
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
    """

    """
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

    # GC experimental spike adaptation data from Brenner...Aldrich, Nat. Neurosci., 2005
    experimental_spike_times = [0., 8.57331572, 21.79656539, 39.24702774, 60.92470277, 83.34214003, 109.5640687,
                                137.1598415, 165.7067371, 199.8546896, 236.2219287, 274.3857332, 314.2404227,
                                355.2575958,
                                395.8520476, 436.7635403]
    experimental_adaptation_indexes = []
    for i in range(3, len(experimental_spike_times) + 1):
        experimental_adaptation_indexes.append(get_adaptation_index(experimental_spike_times[:i]))
    experimental_f_I_slope = 53.  # Hz/ln(pA); rate = slope * ln(current - rheobase)
    # GC experimental f-I data from Kowalski J...Pernia-Andrade AJ, Hippocampus, 2016
    i_inj_increment = 0.05
    num_increments = 10
    context.update(locals())


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


def setup_cell(verbose=False, cvode=False, **kwargs):
    """

    :param verbose: bool
    :param cvode: bool
    """
    cell = DG_GC(neurotree_dict=context.neurotree_dict, mech_file_path=context.mech_file_path,
                 full_spines=context.spines)
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
    axon_seg_locs = [seg.x for seg in cell.axon[2].sec]

    rec_locs = {'soma': 0., 'dend': dend_loc, 'ais': 1., 'axon': axon_seg_locs[0]}
    context.rec_locs = rec_locs
    rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'ais': cell.axon[1], 'axon': cell.axon[2]}
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


def update_mech_dict(x, mech_file_path):
    """

    :param x: array
    :param mech_file_path: str
    """
    for update_func in context.update_params_funcs:
        update_func(x, context)
    context.cell.export_mech_dict(mech_file_path)


def get_stability_features(indiv, c, client_range, export=False):
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
    result = dv.map_async(compute_stability_features, [x], [export])
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result}


def get_fI_features(indiv, c, client_range, export=False):
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
    num_incr = context.num_increments
    i_inj_increment = context.i_inj_increment
    result = dv.map_async(compute_fI_features, [rheobase + i_inj_increment * (i + 1) for i in range(num_incr)],
                          [x] * num_incr, [False] * (num_incr-1) + [True], [export] * num_incr)
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result,
            'filter_features': filter_fI_features}


def filter_fI_features(get_result, old_features, export=False):
    """

    :param get_result: list of dict (each dict has the results from a particular simulation)
    :param old_features: dict
    :param export: bool
    :return: dict
    """
    amps = []
    new_features = {}
    new_features['adi'] = []
    new_features['exp_adi'] = []
    new_features['f_I'] = []
    for i, this_dict in enumerate(get_result):
        amps.append(this_dict['amp'])
        if 'vm_stability' in this_dict.keys():
            new_features['vm_stability'] = this_dict['vm_stability']
        if 'rebound_firing' in this_dict.keys():
            new_features['rebound_firing'] = this_dict['rebound_firing']
        if 'v_min_late' in this_dict.keys():
            new_features['slow_depo'] = this_dict['v_min_late'] - old_features['v_th']
        spike_times = this_dict['spike_times']
        experimental_spike_times = context.experimental_spike_times
        experimental_adaptation_indexes = context.experimental_adaptation_indexes
        stim_dur = context.stim_dur
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
    adapt_ind = range(len(new_features['f_I']))
    adapt_ind.sort(key=amps.__getitem__)
    new_features['adi'] = map(new_features['adi'].__getitem__, adapt_ind)
    new_features['exp_adi'] = map(new_features['exp_adi'].__getitem__, adapt_ind)
    new_features['f_I'] = map(new_features['f_I'].__getitem__, adapt_ind)
    amps = map(amps.__getitem__, adapt_ind)
    if export:
        context.processed_export_file_path = context.export_file_path.replace('.hdf5', '_processed.hdf5')
        with h5py.File(context.processed_export_file_path, 'a') as f:
            group = f.create_group(str(len(f)))
            group.create_dataset('amps', compression='gzip', compression_opts=9, data=amps)
            group.create_dataset('adi', compression='gzip', compression_opts=9, data=new_features['adi'])
            group.create_dataset('exp_adi', compression='gzip', compression_opts=9, data=new_features['exp_adi'])
            group.create_dataset('f_I', compression='gzip', compression_opts=9, data=new_features['f_I'])
            num_increments = context.num_increments
            i_inj_increment = context.i_inj_increment
            rheobase = old_features['rheobase']
            exp_f_I = [context.experimental_f_I_slope * np.log((rheobase + i_inj_increment * (i + 1)) / rheobase)
                          for i in range(num_increments)]
            group.create_dataset('exp_f_I', compression='gzip', compression_opts=9, data=exp_f_I)
    return new_features


def get_objectives(features, objective_names, target_val, target_range):
    """

    :param features: dict
    :param objective_names: list of str
    :param target_val: dict
    :param target_range: dict
    :return: tuple of dict
    """
    if features is None:  # No rheobase value found
        objectives = None
    else:
        objectives = {}
        rheobase = features['rheobase']
        for target in ['v_th', 'ADP', 'AHP', 'spont_firing', 'rebound_firing', 'vm_stability', 'ais_delay',
                       'slow_depo', 'dend_amp', 'soma_peak', 'th_count']:
            # don't penalize AHP or slow_depo less than target
            if not ((target == 'AHP' and features[target] < target_val[target]) or
                        (target == 'slow_depo' and features[target] < target_val[target])):
                objectives[target] = ((target_val[target] - features[target]) / target_range[target]) ** 2.
            else:
                objectives[target] = 0.
        objectives['adi'] = 0.
        for i, this_adi in enumerate(features['adi']):
            if this_adi is not None and features['exp_adi'] is not None:
                objectives['adi'] += ((this_adi - features['exp_adi'][i]) / (0.01 * features['exp_adi'][i])) ** 2.
        features.pop('exp_adi')
        all_adi = []
        for adi in features['adi']:
            if adi is not None:
                all_adi.append(adi)
        features['adi'] = np.mean(all_adi)
        num_increments = context.num_increments
        i_inj_increment = context.i_inj_increment
        target_f_I = [context.experimental_f_I_slope * np.log((rheobase + i_inj_increment * (i + 1)) / rheobase)
                      for i in range(num_increments)]
        f_I_residuals = [(features['f_I'][i] - target_f_I[i]) for i in range(num_increments)]
        features['f_I_residuals'] = np.mean(np.abs(f_I_residuals))
        objectives['f_I_slope'] = 0.
        for i in range(num_increments):
            objectives['f_I_slope'] += (f_I_residuals[i] / (0.01 * target_f_I[i])) ** 2.
        I_inj = [np.log((rheobase + i_inj_increment * (i + 1)) / rheobase) for i in range(num_increments)]
        slope, intercept, r_value, p_value, std_err = stats.linregress(I_inj, features['f_I'])
        features['f_I_slope'] = slope
        features.pop('f_I')
    return features, objectives


def compute_stability_features(x, export=False, plot=False):
    """
    :param x: array
    :param export: bool
    :param plot: bool
    :return: float
    """
    start_time = time.time()
    result = {}
    context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    for update_func in context.update_params_funcs:
        update_func(x, context)

    v_active = context.v_active
    equilibrate = context.equilibrate
    dt = context.dt
    i_th = context.i_th

    soma_vm = offset_vm('soma', v_active)
    result['v_rest'] = soma_vm
    stim_dur = 150.
    context.sim.modify_stim(0, node=context.cell.tree.root, loc=0., dur=stim_dur)
    duration = equilibrate + stim_dur
    context.sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    d_amp = 0.01
    amp = max(0., i_th['soma'] - 0.02)
    while not spike:
        context.sim.modify_stim(0, amp=amp)
        context.sim.run(v_active)
        vm = np.interp(t, context.sim.tvec, context.sim.get_rec('soma')['vec'])
        if np.any(vm[:int(equilibrate/dt)] > -30.):
            if context.disp:
                print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
            return None
        if np.any(vm[int(equilibrate/dt):int((equilibrate+50.)/dt)] > -30.):
            spike = True
        elif amp >= 0.4:
            if context.disp:
                print 'Process %i: Aborting - rheobase outside target range' % (os.getpid())
            return None
        else:
            amp += d_amp
            if context.sim.verbose:
                print 'increasing amp to %.3f' % amp
    context.sim.parameters['amp'] = amp
    context.sim.parameters['description'] = 'spike shape'
    i_th['soma'] = amp
    spike_times = context.cell.spike_detector.get_recordvec().to_python()
    peak, threshold, ADP, AHP = get_spike_shape(vm, spike_times)
    result['v_th'] = threshold
    result['ADP'] = ADP
    result['AHP'] = AHP
    result['rheobase'] = amp
    result['spont_firing'] = len(np.where(spike_times < equilibrate)[0])
    result['th_count'] = len(spike_times)
    dend_vm = np.interp(t, context.sim.tvec, context.sim.get_rec('dend')['vec'])
    th_x = np.where(vm[int(equilibrate / dt):] >= threshold)[0][0] + int(equilibrate / dt)
    if len(spike_times) > 1:
        end = min(th_x + int(10. / dt), int((spike_times[1] - 5.)/dt))
    else:
        end = th_x + int(10. / dt)
    result['soma_peak'] = peak
    dend_peak = np.max(dend_vm[th_x:end])
    dend_pre = np.mean(dend_vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
    result['dend_amp'] = (dend_peak - dend_pre) / (peak - soma_vm)

    # calculate AIS delay
    soma_dvdt = np.gradient(vm, dt)
    ais_vm = np.interp(t, context.sim.tvec, context.sim.get_rec('ais')['vec'])
    ais_dvdt = np.gradient(ais_vm, dt)
    axon_vm = np.interp(t, context.sim.tvec, context.sim.get_rec('axon')['vec'])
    axon_dvdt = np.gradient(axon_vm, dt)
    left = th_x - int(2. / dt)
    right = th_x + int(5. / dt)
    soma_peak = np.max(soma_dvdt[left:right])
    soma_peak_t = np.where(soma_dvdt[left:right] == soma_peak)[0][0] * dt
    ais_peak = np.max(ais_dvdt[left:right])
    ais_peak_t = np.where(ais_dvdt[left:right] == ais_peak)[0][0] * dt
    axon_peak = np.max(axon_dvdt[left:right])
    axon_peak_t = np.where(axon_dvdt[left:right] == axon_peak)[0][0] * dt
    result['ais_delay'] = max(0., ais_peak_t + dt - soma_peak_t) + max(0., ais_peak_t + dt - axon_peak_t)
    if context.disp:
        print 'Process %i took %.1f s to find spike rheobase at amp: %.3f' % (os.getpid(), time.time() - start_time,
                                                                              amp)
    if plot:
        context.sim.plot()
    if export:
        export_sim_results()
    return result


def compute_fI_features(amp, x, extend_dur=False, export=False, plot=False):
    """

    :param amp: float
    :param x: array
    :param extend_dur: bool
    :param export: bool
    :param plot: bool
    :return: dict
    """
    context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    for update_func in context.update_params_funcs:
        update_func(x, context)
        sys.stdout.flush()

    soma_vm = offset_vm('soma', context.v_active)
    context.sim.parameters['amp'] = amp
    context.sim.parameters['description'] = 'f_I'
    start_time = time.time()

    stim_dur = context.stim_dur
    equilibrate = context.equilibrate
    v_active = context.v_active
    dt = context.dt

    context.sim.modify_stim(0, node=context.cell.tree.root, loc=0., dur=stim_dur, amp=amp)
    if extend_dur:
        duration = equilibrate + stim_dur + 100. #extend duration of simulation to find rebound
    else:
        duration = equilibrate + stim_dur
    context.sim.tstop = duration
    context.sim.run(v_active)
    if plot:
        context.sim.plot()
    spike_times = np.subtract(context.cell.spike_detector.get_recordvec().to_python(), equilibrate)
    t = np.arange(0., duration, dt)
    result = {}
    result['spike_times'] = spike_times
    result['amp'] = amp
    if extend_dur:
        vm = np.interp(t, context.sim.tvec, context.sim.get_rec('soma')['vec'])
        v_min_late = np.min(vm[int((equilibrate + stim_dur - 20.) / dt):int((equilibrate + stim_dur - 1.) / dt)])
        result['v_min_late'] = v_min_late
        v_rest = np.mean(vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
        v_after = np.max(vm[-int(50. / dt):-1])
        vm_stability = abs(v_after - v_rest)
        result['vm_stability'] = vm_stability
        result['rebound_firing'] = len(np.where(spike_times > stim_dur)[0])
    if context.disp:
        print 'Process %i took %.1f s to run simulation with I_inj amp: %.3f' % (os.getpid(), time.time() - start_time,
                                                                                 amp)
        sys.stdout.flush()
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


def update_na_ka_stability(x, local_context=None):
    """
    :param x: array ['soma.gbar_nas', 'dend.gbar_nas', 'dend.gbar_nas slope', 'dend.gbar_nas min', 'dend.gbar_nas bc',
                    'axon.gbar_nax', 'ais.gbar_nax', 'soma.gkabar', 'dend.gkabar', 'soma.gkdrbar', 'axon.gkabar',
                    'soma.sh_nas/x', 'ais.sha_nas', 'soma.gCa factor', 'soma.gCadepK factor', 'soma.gkmbar', 'ais.gkmbar']
    """
    if local_context is None:
        local_context = context
    cell = local_context.cell
    param_indexes = local_context.param_indexes
    default_params = local_context.default_params
    if local_context.spines is False:
        cell.correct_for_spines()
    cell.modify_mech_param('soma', 'nas', 'gbar', find_param_value('soma.gbar_nas', x, param_indexes, default_params))
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', find_param_value('soma.gkdrbar', x, param_indexes, default_params))
    cell.modify_mech_param('soma', 'kap', 'gkabar', find_param_value('soma.gkabar', x, param_indexes, default_params))
    slope = (find_param_value('dend.gkabar', x, param_indexes, default_params) -
             find_param_value('soma.gkabar', x, param_indexes, default_params)) / 300.
    cell.modify_mech_param('soma', 'nas', 'sh', find_param_value('soma.sh_nas/x', x, param_indexes, default_params))
    for sec_type in ['apical']:
        cell.reinitialize_subset_mechanisms(sec_type, 'nas')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=find_param_value('soma.gkabar', x, param_indexes, default_params)+slope*75.,
                               replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300.,
                               value=find_param_value('soma.gkabar', x, param_indexes, default_params)+slope*300.,
                               replace=False)
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
        cell.modify_mech_param(sec_type, 'nas', 'sha', 0.)  # 5.)
        cell.modify_mech_param(sec_type, 'nas', 'gbar',
                               find_param_value('dend.gbar_nas', x, param_indexes, default_params))
        cell.modify_mech_param(sec_type, 'nas', 'gbar', origin='parent',
                               slope=find_param_value('dend.gbar_nas slope', x, param_indexes, default_params),
                               min=find_param_value('dend.gbar_nas min', x, param_indexes, default_params),
                               custom={'method': 'custom_gradient_by_branch_ord', 'branch_order':
                                   find_param_value('dend.gbar_nas bo', x, param_indexes, default_params)},
                               replace=False)
        cell.modify_mech_param(sec_type, 'nas', 'gbar', origin='parent',
                               slope=find_param_value('dend.gbar_nas slope', x, param_indexes, default_params),
                               min=find_param_value('dend.gbar_nas min', x, param_indexes, default_params),
                               custom={'method': 'custom_gradient_by_terminal'}, replace=False)
    cell.reinitialize_subset_mechanisms('axon_hill', 'kap')
    cell.reinitialize_subset_mechanisms('axon_hill', 'kdr')
    cell.modify_mech_param('ais', 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param('ais', 'kap', 'gkabar', find_param_value('axon.gkabar', x, param_indexes, default_params))
    cell.modify_mech_param('axon', 'kdr', 'gkdrbar', origin='ais')
    cell.modify_mech_param('axon', 'kap', 'gkabar', origin='ais')
    cell.modify_mech_param('axon_hill', 'nax', 'sh',
                           find_param_value('soma.sh_nas/x', x, param_indexes, default_params))
    cell.modify_mech_param('axon_hill', 'nax', 'gbar', local_context.soma_na_gbar)
    cell.modify_mech_param('axon', 'nax', 'gbar', find_param_value('axon.gbar_nax', x, param_indexes, default_params))
    for sec_type in ['ais', 'axon']:
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='axon_hill')
    cell.modify_mech_param('soma', 'Ca', 'gcamult',
                           find_param_value('soma.gCa factor', x, param_indexes, default_params))
    cell.modify_mech_param('soma', 'CadepK', 'gcakmult', find_param_value('soma.gCadepK factor', x, param_indexes,
                                                                          default_params))
    cell.modify_mech_param('soma', 'km3', 'gkmbar', find_param_value('soma.gkmbar', x, param_indexes, default_params))
    cell.modify_mech_param('ais', 'km3', 'gkmbar', find_param_value('ais.gkmbar', x, param_indexes, default_params))
    cell.modify_mech_param('axon_hill', 'km3', 'gkmbar', origin='soma')
    cell.modify_mech_param('axon', 'km3', 'gkmbar', origin='ais')
    cell.modify_mech_param('ais', 'nax', 'sha', find_param_value('ais.sha_nas', x, param_indexes, default_params))
    cell.modify_mech_param('ais', 'nax', 'gbar', find_param_value('ais.gbar_nax', x, param_indexes, default_params))


def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(context.rec_file_path, 'a') as f:
        context.sim.export_to_file(f)
