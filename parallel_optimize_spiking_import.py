__author__ = 'Grace Ng'
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
default_group_sizes = (1, 10)

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
default_feature_names = ['v_rest', 'v_th', 'ADP', 'AHP', 'spont_firing', 'rebound_firing', 'vm_stability', 'ais_delay',
                         'slow_depo', 'dend_amp', 'rheobase', 'soma_peak', 'adi', 'f_I_slope', 'f_I_residuals']
default_objective_names = ['v_rest', 'v_th', 'ADP', 'AHP', 'spont_firing', 'rebound_firing', 'vm_stability', 'ais_delay',
                           'slow_depo', 'dend_amp', 'soma_peak', 'adi', 'f_I_slope']
default_target_val = {'v_rest': v_init, 'v_th': -48., 'ADP': 0., 'AHP': 4., 'spont_firing': 0, 'rebound_firing': 0,
                      'vm_stability': 0., 'ais_delay': 0., 'slow_depo': 20., 'dend_amp': 0.3, 'soma_peak': 40.}
default_target_range = {'v_rest': 0.25, 'v_th': .01, 'ADP': 0.01, 'AHP': .005, 'spont_firing': 1, 'rebound_firing': 1,
                      'vm_stability': 1., 'ais_delay': 0.0005, 'slow_depo': 0.5, 'dend_amp': 0.0002, 'soma_peak': 2.}
default_optimization_title = 'spiking'

# param_file_path = 'data/optimize_spiking_defaults.yaml'


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
def get_stability_features(indiv, c, client_range, export=False):
    """
    Distribute simulations across available engines for testing spike stability.
    :param indiv: dict {'pop_id': pop_id, 'x': x arr, 'features': features dict}
    :param client_range: list of ints
    :param export: False (for exporting voltage traces)
    :return: dict
    """
    x = indiv['x']
    dv = c[client_range]
    result = dv.map_async(spike_shape_features, [x], [export])
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result}

@interactive
def process_stability_features(get_result, old_features):
    """

    :param get_result: list with one dict
    :return: dict
    """
    return get_result[0]

@interactive
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
    sys.stdout.flush()
    # Calculate firing rates for a range of I_inj amplitudes using a stim duration of 500 ms
    result = dv.map_async(sim_f_I, [rheobase + i_inj_increment * (i + 1) for i in range(num_increments)],
                          [x] * num_increments, [False] * (num_increments-1) + [True], [export] * num_increments)
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result,
            'process_features': process_fI_features}

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
        if 'vm_stability' in this_dict.keys():
            new_features['vm_stability'] = this_dict['vm_stability']
        if 'rebound_firing' in this_dict.keys():
            new_features['rebound_firing'] = this_dict['rebound_firing']
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

    adapt_ind = range(len(new_features['f_I']))
    adapt_ind.sort(key=temp_dict['amp'].__getitem__)
    new_features['adi'] = map(new_features['adi'].__getitem__, adapt_ind)
    new_features['exp_adi'] = map(new_features['exp_adi'].__getitem__, adapt_ind)
    new_features['f_I'] = map(new_features['f_I'].__getitem__, adapt_ind)
    return new_features

@interactive
def get_spiking_objectives(features, objective_names, target_val, target_range):
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
                       'slow_depo', 'dend_amp', 'soma_peak']:
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
        target_f_I = [experimental_f_I_slope * np.log((rheobase + i_inj_increment * (i + 1)) / rheobase)
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

@interactive
def spike_shape_features(x, export=False, plot=False):
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
    sim.parameters['description'] = 'spike shape'
    i_th['soma'] = amp
    spike_times = cell.spike_detector.get_recordvec().to_python()
    peak, threshold, ADP, AHP = get_spike_shape(vm, spike_times)
    result['v_th'] = threshold
    result['ADP'] = ADP
    result['AHP'] = AHP
    result['rheobase'] = amp
    result['spont_firing'] = len(np.where(spike_times < equilibrate))
    dend_vm = np.interp(t, sim.tvec, sim.get_rec('dend')['vec'])
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
    if export:
        export_sim_results()
    return result


@interactive
def sim_f_I(amp, x, extend_dur=False, export=False, plot=False):
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
    sim.parameters['description'] = 'f_I'
    start_time = time.time()
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=stim_dur, amp=amp)
    if extend_dur:
        duration = equilibrate + stim_dur + 100. #extend duration of simulation to find rebound
    else:
        duration = equilibrate + stim_dur
    sim.tstop = duration
    print 'starting sim at %.1f A' %amp
    sys.stdout.flush()
    sim.run(v_active)
    if plot:
        sim.plot()
    spike_times = np.subtract(cell.spike_detector.get_recordvec().to_python(), equilibrate)
    t = np.arange(0., duration, dt)
    result = {}
    result['spike_times'] = spike_times
    result['amp'] = amp
    rate = len(spike_times) / stim_dur * 1000.
    result['rate'] = rate
    vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
    v_min_late = np.min(vm[int((equilibrate + stim_dur - 20.) / dt):int((equilibrate + stim_dur - 1.) / dt)])
    result['v_min_late'] = v_min_late
    if extend_dur:
        v_rest = np.mean(vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
        v_after = np.max(vm[-int(50. / dt):-1])
        vm_stability = abs(v_after - v_rest)
        result['vm_stability'] = vm_stability
        result['rebound_firing'] = len(np.where(spike_times > stim_dur))
    print 'Process %i took %.1f s to run simulation with I_inj amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    sys.stdout.flush()
    if export:
        export_sim_results()
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
def update_na_ka_stability(x):
    """

    :param x: array ['soma.gbar_nas', 'dend.gbar_nas', 'axon.gbar_nax', 'ais.gbar_nax', 'soma.gkabar', 'dend.gkabar',
                       'soma.gkdrbar', 'axon.gkbar', 'soma.sh_nas/x', 'ais.sha_nas', 'soma.gCa factor',
                       'soma.gCadepK factor', 'soma.gkmbar', 'ais.gkmbar']
    """
    if spines is False:
        cell.reinit_mechanisms(reset_cable=True)
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
    cell.modify_mech_param('ais', 'kap', 'gkabar', x[param_indexes['axon.gkbar']])
    cell.modify_mech_param('axon', 'kdr', 'gkdrbar', origin='ais')
    cell.modify_mech_param('axon', 'kap', 'gkabar', origin='ais')
    cell.modify_mech_param('axon_hill', 'nax', 'sh', x[param_indexes['soma.sh_nas/x']])
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
def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f)

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