__author__ = 'Grace Ng'
from specify_cells3 import *
from plot_results import *

"""
Iterates through every spine and activates AMPA_KIN synapses. Allows measurement of EPSP attenuation
and kinetics.

Import this script into parallel_optimize_main. Then, set up ipcluster and run parallel_optimize_main.py
"""
# param_file_path = 'data/optimize_EPSP_attenuation_defaults.yaml'

context = Context()

def config_engine(param_names, mech_file_path, neurotree_dict, spines, rec_filepath):
    """

    :param param_names: list of str
    :param mech_file_path: str
    :param neurotree_dict: dict
    :param spines: bool
    :param rec_filepath: str
    :return:
    """
    param_indexes = {param_name: i for i, param_name in enumerate(param_names)}
    prev_job_type = 'None'
    context.update(locals())

def init_EPSP_attenuation_engine():
    """

    :return:
    """
    equilibrate = 250.  # time to steady-state
    duration = 450.
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

def setup_cell():
    """

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
    axon_seg_locs = [seg.x for seg in cell.axon[2].sec]

    rec_locs = {'soma': 0., 'dend': dend_loc, 'ais': 1., 'axon': axon_seg_locs[0]}
    context.rec_locs = rec_locs
    rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'ais': cell.axon[1], 'axon': cell.axon[2]}
    context.rec_nodes = rec_nodes

    equilibrate = context.equilibrate
    stim_dur = context.stim_dur
    duration = context.duration
    dt = context.dt

    sim = QuickSim(duration, cvode=False, dt=dt, verbose=False)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)
    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    sim.parameters['spines'] = context.spines
    context.sim = sim

    context.spike_output_vec = h.Vector()
    cell.spike_detector.record(context.spike_output_vec)

def update_mech_dict(x, update_function, mech_file_path):
    update_function(x)
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
    init_spiking_engine() # So that controller core has access to constansts such as num_increments
    num_incr = context.num_increments
    i_inj_increment = context.i_inj_increment
    result = dv.map_async(compute_fI_features, [rheobase + i_inj_increment * (i + 1) for i in range(num_incr)],
                          [x] * num_incr, [False] * (num_incr-1) + [True], [export] * num_incr)
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result,
            'filter_features': filter_fI_features}

def filter_fI_features(get_result, old_features):
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
    adapt_ind.sort(key=temp_dict['amp'].__getitem__)
    new_features['adi'] = map(new_features['adi'].__getitem__, adapt_ind)
    new_features['exp_adi'] = map(new_features['exp_adi'].__getitem__, adapt_ind)
    new_features['f_I'] = map(new_features['f_I'].__getitem__, adapt_ind)
    return new_features

def get_objectives(features, objective_names, target_val, target_range):
    """

    :param features: dict
    :return: dict
    """
    if features is None: #No rheobase value found
        objectives = None
    else:
        objectives = {}
        for objective in objective_names:
            objectives[objective] = 0.
        rheobase = features['rheobase']
        for target in ['v_rest', 'v_th', 'ADP', 'AHP', 'spont_firing', 'rebound_firing', 'vm_stability', 'ais_delay',
                       'slow_depo', 'dend_amp', 'soma_peak', 'th_count', 'na_gradient']:
            # don't penalize AHP or slow_depo less than target
            if not ((target == 'AHP' and features[target] < target_val[target]) or
                        (target == 'slow_depo' and features[target] < target_val[target])):
                objectives[target] = ((target_val[target] - features[target]) / target_range[target]) ** 2.
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
        features['f_I_residuals'] = np.mean(f_I_residuals)
        for i in range(num_increments):
            objectives['f_I_slope'] += (f_I_residuals[i] / (0.01 * target_f_I[i])) ** 2.
        I_inj = [np.log((rheobase + i_inj_increment * (i + 1)) / rheobase) for i in range(num_increments)]
        slope, intercept, r_value, p_value, std_err = stats.linregress(I_inj, features['f_I'])
        features['f_I_slope'] = slope
        features.pop('f_I')
    return features, objectives

def compute_stability_features(x, export=False, plot=False):
    """
    :param local_x: array [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor, soma.sh_nas/x,
                    axon.gkbar factor, dend.gkabar factor]
    :param plot: bool
    :return: float
    """
    na_gradient = max(0., x[context.param_indexes['dend.gbar_nas']] - x[context.param_indexes['soma.gbar_nas']]) + \
                  max(0., x[context.param_indexes['soma.gbar_nas']] - x[context.param_indexes['ais.gbar_nax']]) + \
                  max(0., x[context.param_indexes['soma.gbar_nas']] - x[context.param_indexes['axon.gbar_nax']]) + \
                  max(0., x[context.param_indexes['axon.gbar_nax']] - x[context.param_indexes['ais.gbar_nax']])
    result = {'na_gradient': na_gradient}
    start_time = time.time()
    if context.prev_job_type == 'spiking':
        context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    else:
        init_spiking_engine()
        setup_cell()
    update_na_ka_stability(x)
    # sim.cvode_state = True

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
            print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
            return None
        if np.any(vm[int(equilibrate/dt):int((equilibrate+50.)/dt)] > -30.):
            spike = True
        elif amp >= 0.4:
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
    dend_pre = np.mean(dend_vm[th_x - int(0.2 / dt):th_x - int(0.1 / dt)])
    result['dend_amp'] = (dend_peak - dend_pre) / (peak - threshold)

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
    print 'Process %i took %.1f s to find spike rheobase at amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    if plot:
        context.sim.plot()
    if export:
        export_sim_results()
    context.prev_job_type = 'spiking'
    return result

def compute_fI_features(amp, x, extend_dur=False, export=False, plot=False):
    """

    :param amp: float
    :param local_x: array
    :param plot: bool
    :return: dict
    """
    if context.prev_job_type == 'spiking':
        context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    else:
        init_spiking_engine()
        setup_cell()
    update_na_ka_stability(x)
    # sim.cvode_state = True
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
    print 'starting sim at %.1f A' %amp
    sys.stdout.flush()
    context.sim.run(v_active)
    if plot:
        context.sim.plot()
    spike_times = np.subtract(context.cell.spike_detector.get_recordvec().to_python(), equilibrate)
    t = np.arange(0., duration, dt)
    result = {}
    result['spike_times'] = spike_times
    result['amp'] = amp
    rate = len(spike_times) / stim_dur * 1000.
    result['rate'] = rate
    if extend_dur:
        vm = np.interp(t, context.sim.tvec, context.sim.get_rec('soma')['vec'])
        v_min_late = np.min(vm[int((equilibrate + stim_dur - 20.) / dt):int((equilibrate + stim_dur - 1.) / dt)])
        result['v_min_late'] = v_min_late
        v_rest = np.mean(vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
        v_after = np.max(vm[-int(50. / dt):-1])
        vm_stability = abs(v_after - v_rest)
        result['vm_stability'] = vm_stability
        result['rebound_firing'] = len(np.where(spike_times > stim_dur)[0])
    print 'Process %i took %.1f s to run simulation with I_inj amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    sys.stdout.flush()
    if export:
        export_sim_results()
    context.prev_job_type = 'spiking'
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

def update_na_ka_stability(x):
    """

    :param x: array ['soma.gbar_nas', 'dend.gbar_nas', 'axon.gbar_nax', 'ais.gbar_nax', 'soma.gkabar', 'dend.gkabar',
                       'soma.gkdrbar', 'axon.gkbar', 'soma.sh_nas/x', 'ais.sha_nas', 'soma.gCa factor',
                       'soma.gCadepK factor', 'soma.gkmbar', 'ais.gkmbar']
    """
    param_indexes = context.param_indexes
    context.cell.modify_mech_param('soma', 'nas', 'gbar', x[param_indexes['soma.gbar_nas']])
    context.cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[param_indexes['soma.gkdrbar']])
    context.cell.modify_mech_param('soma', 'kap', 'gkabar', x[param_indexes['soma.gkabar']])
    slope = (x[param_indexes['dend.gkabar']] - x[param_indexes['soma.gkabar']]) / 300.
    context.cell.modify_mech_param('soma', 'nas', 'sh', x[param_indexes['soma.sh_nas/x']])
    for sec_type in ['apical']:
        context.cell.reinitialize_subset_mechanisms(sec_type, 'nas')
        context.cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        context.cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        context.cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        context.cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=x[param_indexes['soma.gkabar']]+slope*75., replace=False)
        context.cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300.,
                               value=x[param_indexes['soma.gkabar']]+slope*300., replace=False)
        context.cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
        context.cell.modify_mech_param(sec_type, 'nas', 'sha', 5.)
        context.cell.modify_mech_param(sec_type, 'nas', 'gbar', x[param_indexes['dend.gbar_nas']])
    context.cell.set_terminal_branch_na_gradient()
    context.cell.reinitialize_subset_mechanisms('axon_hill', 'kap')
    context.cell.reinitialize_subset_mechanisms('axon_hill', 'kdr')
    context.cell.modify_mech_param('ais', 'kdr', 'gkdrbar', origin='soma')
    context.cell.modify_mech_param('ais', 'kap', 'gkabar', x[param_indexes['axon.gkbar']])
    context.cell.modify_mech_param('axon', 'kdr', 'gkdrbar', origin='ais')
    context.cell.modify_mech_param('axon', 'kap', 'gkabar', origin='ais')
    context.cell.modify_mech_param('axon_hill', 'nax', 'sh', x[param_indexes['soma.sh_nas/x']])
    context.cell.modify_mech_param('axon_hill', 'nax', 'gbar', context.soma_na_gbar)
    context.cell.modify_mech_param('axon', 'nax', 'gbar', x[param_indexes['axon.gbar_nax']])
    for sec_type in ['ais', 'axon']:
        context.cell.modify_mech_param(sec_type, 'nax', 'sh', origin='axon_hill')
    context.cell.modify_mech_param('soma', 'Ca', 'gcamult', x[param_indexes['soma.gCa factor']])
    context.cell.modify_mech_param('soma', 'CadepK', 'gcakmult', x[param_indexes['soma.gCadepK factor']])
    context.cell.modify_mech_param('soma', 'km3', 'gkmbar', x[param_indexes['soma.gkmbar']])
    context.cell.modify_mech_param('ais', 'km3', 'gkmbar', x[param_indexes['ais.gkmbar']])
    context.cell.modify_mech_param('axon_hill', 'km3', 'gkmbar', origin='soma')
    context.cell.modify_mech_param('axon', 'km3', 'gkmbar', origin='ais')
    context.cell.modify_mech_param('ais', 'nax', 'sha', x[param_indexes['ais.sha_nas']])
    context.cell.modify_mech_param('ais', 'nax', 'gbar', x[param_indexes['ais.gbar_nax']])
    if context.spines is False:
        context.cell.correct_for_spines()
    context.cell.set_terminal_branch_na_gradient()

def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(context.rec_filepath, 'a') as f:
        context.sim.export_to_file(f)