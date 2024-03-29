__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
import os
"""

Extend linear kap gradient into basals and obliques, see how it does. Aim for 60% spike attenuation ~200 uM along trunk.

Aims for spike initiation at initial segment by increasing nax density and decreasing activation V1/2 relative to soma,
axon_hill, and axon compartments.

Hierarchical optimization:
I) optimize g_pas for target rinp at soma, trunk 50 uM from SLM border, and tuft 50 uM into SLM [without h]
II) optimize ghbar_h for target rinp trunk 50 uM from SLM border, and tuft 50 uM into SLM
# III) optimize e_pas and vhalfl_h for target v_rest
IV) optimize gbar_nax/nas/sh, gkabar_kap/d, gkdrbar for target na spike threshold and AHP amp
"""

morph_filename = 'EB2-late-bifurcation.swc'

#mech_filename = '080615 rebalanced na_ka ampa nmda - EB2'
mech_filename = '100215 template dendritic excitability'


def check_bounds(x, param_name):
    """
    For optimize_polish, based on simplex algorithm, check that the current set of parameters are within the bounds.
    :param x: array
    :param param_name: str
    :return: bool
    """
    for i in range(len(x)):
        if ((xmin[param_name][i] is not None and x[i] < xmin[param_name][i]) or
            (xmax[param_name][i] is not None and x[i] > xmax[param_name][i])):
            return False
    return True


def update_pas(x):
    """

    x0 = [2.28e-05, 1.58e-06, 58.4]
    :param x: array (soma.g_pas, trunk.g_pas slope, trunk.g_pas tau)
    """
    e_k = -90.  # cell.mech_dict['soma']['ions']['ek']['value']
    distance = cell.get_distance_to_node(cell.tree.root, distal_trunk, 1.)
    slope = (e_k - v_init) / (np.exp(distance / x[2]) - 1.)
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.modify_mech_param('trunk', 'pas', 'g', origin='soma', slope=x[1], tau=x[2])
    cell.modify_mech_param('trunk', 'pas', 'e', origin='soma', slope=slope, tau=x[2])
    for sec_type in ['apical', 'tuft']:
        cell.modify_mech_param(sec_type, 'pas', 'e', origin='trunk')
    for sec_type in ['basal', 'axon_hill', 'axon', 'ais', 'trunk', 'apical', 'tuft', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')


def update_h(x):
    """

    x0 = [1.07e-09, 0.0204, 28.7, 265.5]
    :param x: array (soma.ghbar_h, trunk.ghbar_h slope, trunk.ghbar_h tau, trunk.ghbar_h xhalf)
    """
    cell.modify_mech_param('soma', 'h', 'ghbar', x[0])
    cell.modify_mech_param('trunk', 'h', 'ghbar', origin='soma', slope=x[1], tau=x[2], xhalf=x[3])
    cell.modify_mech_param('basal', 'h', 'ghbar', origin='soma')
    for sec_type in ['apical', 'tuft']:
        cell.modify_mech_param(sec_type, 'h', 'ghbar', origin='trunk')


def update_v_rest(x):
    """

    :param x: array (soma.ghbar_h, trunk.ghbar_h slope, trunk.ghbar_h tau, trunk.ghbar_h xhalf)
    :return:
    """
    e_pas = cell.mech_dict['soma']['pas']['e']['value']
    e_k = cell.mech_dict['soma']['ions']['ek']['value']
    distance = cell.get_distance_to_node(cell.tree.root, distal_trunk, 1.)
    slope = (e_k - e_pas) / ((1. / (1. + np.exp((x[3] - distance) / x[2]))) - (1. / (1. + np.exp(x[3] / x[2]))))
    cell.modify_mech_param('trunk', 'pas', 'e', origin='soma', slope=slope, tau=x[2], xhalf=x[3])
    for sec_type in ['apical', 'tuft']:
        cell.modify_mech_param(sec_type, 'pas', 'e', origin='trunk')


def update_na_ka(x):
    """

    :param x: array [soma.sh_nas, soma.gkabar, dend.gkabar factor, soma.gkdrbar, ais.gbar_nax factor, axon.kap factor,
                    axon.kdr factor]
    """
    cell.modify_mech_param('soma', 'nas', 'gbar', soma_na_gbar)
    cell.modify_mech_param('soma', 'nas', 'sh', x[0])
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[3])
    cell.modify_mech_param('soma', 'kap', 'gkabar', x[1])
    slope = (x[2]-1.)*x[1]/300.
    for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=x[1]+slope*75., replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300., value=x[1]*x[2], replace=False)
        cell.modify_mech_param(sec_type, 'nas', 'sh', origin='soma')
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param('trunk', 'nas', 'gbar', origin='soma')
    cell.modify_mech_param('basal', 'nas', 'gbar', origin='soma', slope=-0.04/200., min=0.)
    cell.modify_mech_param('apical', 'nas', 'gbar', origin='trunk', slope=-0.04/200., min=0.)
    cell.modify_mech_param('tuft', 'nas', 'gbar', origin='trunk')
    cell.modify_mech_param('axon_hill', 'nax', 'gbar', 2. * soma_na_gbar)
    cell.modify_mech_param('axon_hill', 'nax', 'sh', x[0])
    cell.modify_mech_param('axon_hill', 'kap', 'gkabar', x[1] * x[5])
    cell.modify_mech_param('axon_hill', 'kdr', 'gkdrbar', x[3] * x[6])
    cell.modify_mech_param('ais', 'nax', 'gbar', 2. * soma_na_gbar * x[4])
    cell.modify_mech_param('axon', 'nax', 'gbar', origin='axon_hill')
    for sec_type in ['ais', 'axon']:
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='axon_hill')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='axon_hill')
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='axon_hill')


def zero_h():
    """

    """
    cell.modify_mech_param('soma', 'h', 'ghbar', 0.)
    cell.mech_dict['trunk']['h']['ghbar']['value'] = 0.
    cell.mech_dict['trunk']['h']['ghbar']['slope'] = 0.
    for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
        cell.reinitialize_subset_mechanisms(sec_type, 'h')


def zero_na():
    """

    """
    for sec_type in ['axon_hill', 'ais']:
        cell.modify_mech_param(sec_type, 'nax', 'gbar', 0.)
    cell.reinitialize_subset_mechanisms('axon', 'nax')
    cell.modify_mech_param('soma', 'nas', 'gbar', 0.)
    for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
        cell.reinitialize_subset_mechanisms(sec_type, 'nas')


def get_spike_shape(vm):
    """

    :param t: array
    :param vm: array
    :return: tuple of float: (v_peak, th_v, ADP, AHP)
    """
    vm = vm[int((equilibrate+0.4)/dt):]
    dvdt = np.gradient(vm, dt)
    th_x = np.where(dvdt > th_dvdt)[0]
    if th_x.any():
        th_x = th_x[0] - int(1.6/dt)
    else:
        th_x = np.where(vm > -30.)[0][0] - int(1.8/dt)
    th_v = vm[th_x]
    v_before = np.mean(vm[th_x-int(0.1/dt):th_x])
    v_peak = np.max(vm[th_x:th_x+int(5./dt)])
    x_peak = np.where(vm[th_x:th_x+int(5./dt)]==v_peak)[0][0]
    v_AHP = np.min(vm[th_x+x_peak:th_x+int(12./dt)])
    x_AHP = np.where(vm[th_x+x_peak:th_x+int(12./dt)]==v_AHP)[0][0]
    AHP = v_before - v_AHP
    # if spike waveform includes an ADP before an AHP, return the value of the ADP in order to increase error function
    rising_x = np.where(dvdt[th_x+x_peak:th_x+x_peak+x_AHP] > 0.)[0]
    if rising_x.any():
        v_ADP = np.max(vm[th_x+x_peak+rising_x[0]:th_x+x_peak+x_AHP])
        ADP = v_ADP - v_AHP
    else:
        ADP = 0.
    return v_peak, th_v, ADP, AHP


def offset_vm(sec_type):
    """

    :param sec_type: str
    """
    sim.modify_stim(0, amp=0.)
    if sec_type == 'soma':
        sim.modify_stim(1, node=cell.tree.root, loc=0., amp=i_holding[sec_type])
        rec = sim.rec_list[0]['vec']
    elif sec_type == 'trunk':
        sim.modify_stim(1, node=trunk, loc=1., amp=i_holding[sec_type])
        rec = sim.rec_list[1]['vec']
    elif sec_type == 'tuft':
        sim.modify_stim(1, node=tuft, loc=0., amp=i_holding[sec_type])
        rec = sim.rec_list[2]['vec']
    offset = True
    sim.tstop = equilibrate
    sim.run(v_init)
    t = np.arange(0, equilibrate, dt)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < v_init - 0.5:
        i_holding[sec_type] += 0.025
        while offset:
            if sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (i_holding[sec_type], sec_type)
            sim.modify_stim(1, amp=i_holding[sec_type])
            sim.run(v_init)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest < v_init - 0.5:
                i_holding[sec_type] += 0.025
            else:
                offset = False
    elif v_rest > v_init + 0.5:
        i_holding[sec_type] -= 0.025
        while offset:
            if sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (i_holding[sec_type], sec_type)
            sim.modify_stim(1, amp=i_holding[sec_type])
            sim.run(v_init)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > v_init + 0.5:
                i_holding[sec_type] -= 0.025
            else:
                offset = False
    sim.tstop = duration
    return initial_v_rest


def pas_error(x, plot=0):
    """
    :param x: array (soma.g_pas, trunk.g_pas slope, trunk.g_pas tau)
    :param plot: int
    :return: float
    """
    if np.any(x < 0.):
        return 1e9
    start_time = time.time()
    sim.tstop = duration
    result = {}
    amp = -0.15
    zero_h()
    zero_na()
    update_pas(x)
    offset_vm('soma')
    sim.modify_stim(0, node=cell.tree.root, loc=0., amp=amp, dur=stim_dur)
    sim.run(v_init)
    result['soma'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[2]
    offset_vm('trunk')
    sim.modify_stim(0, node=trunk, loc=1., amp=amp)
    sim.run(v_init)
    result['trunk'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[1]['vec']), equilibrate, duration, amp)[2]
    offset_vm('tuft')
    sim.modify_stim(0, node=tuft, loc=0., amp=amp)
    sim.run(v_init)
    result['tuft'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[2]['vec']), equilibrate, duration, amp)[2]
    sim.modify_stim(1, amp=0.)
    Err = 0.
    for target in result:
        Err += ((target_val['pas'][target] - result[target])/target_range['pas'][target])**2.
    print('Simulation took %.3f s' % (time.time()-start_time))
    print 'Process %i: [soma g_pas, trunk slope, trunk tau]: [%.2E, %.2E, %.1f], soma R_inp: %.1f, trunk R_inp: %.1f,' \
          ' tuft R_inp: %.1f' % (os.getpid(), x[0], x[1], x[2], result['soma'], result['trunk'], result['tuft'])
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if plot:
        sim.plot()
    else:
        return Err


def h_error(x, plot=0):
    """

    :param x: array (soma.ghbar_h, trunk.ghbar_h slope, trunk.ghbar_h tau, trunk.ghbar_h xhalf)
    :param plot: int
    :return: float
    """
    if np.any(x < 1.e-10) or x[2] < 25.:  # or x[2] + x[3] > 300.:
        return 1e9
    start_time = time.time()
    sim.tstop = duration
    result = {}
    amp = -0.15
    zero_na()
    update_h(x)
    #update_v_rest(x)
    soma_vm = offset_vm('soma')
    sim.modify_stim(0, node=cell.tree.root, loc=0., amp=amp, dur=stim_dur)
    sim.run(v_init)
    result['soma'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[2]
    offset_vm('trunk')
    sim.modify_stim(0, node=trunk, loc=1., amp=amp)
    sim.run(v_init)
    result['trunk'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[1]['vec']), equilibrate, duration, amp)[2]
    tuft_vm = offset_vm('tuft')
    sim.modify_stim(0, node=tuft, loc=0., amp=amp)
    sim.run(v_init)
    result['tuft'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[2]['vec']), equilibrate, duration, amp)[2]
    # result['v_rest_offset'] = tuft_vm - soma_vm
    v_rest_offset = tuft_vm - soma_vm
    sim.modify_stim(1, amp=0.)
    Err = 0.
    for target in result:
        Err += ((target_val['h'][target] - result[target])/target_range['h'][target])**2.
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [soma.ghbar, trunk ghbar slope, trunk ghbar tau, trunk ghbar xhalf]: [%.2E, %.2E, %.1f, %.1f],' \
          ' soma R_inp: %.1f, trunk R_inp: %.1f, tuft R_inp: %.1f, v_rest offset: %.1f' % (os.getpid(), x[0], x[1],
                                x[2], x[3], result['soma'], result['trunk'], result['tuft'], v_rest_offset)
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if plot:
        sim.plot()
    else:
        return Err


def v_rest_error(plot=0):
    """

    :param plot: int
    :return: float
    """
    start_time = time.time()
    sim.tstop = equilibrate
    result = {}
    sim.modify_stim(0, amp=0.)
    sim.modify_stim(1, amp=0.)
    sim.run(v_init)
    t = np.arange(0, equilibrate, dt)
    vm = np.interp(t, sim.tvec, sim.rec_list[0]['vec'])
    result['soma'] = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    vm = np.interp(t, sim.tvec, sim.rec_list[1]['vec'])
    trunk_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    vm = np.interp(t, sim.tvec, sim.rec_list[2]['vec'])
    tuft_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    result['tuft_offset'] = tuft_rest - result['soma']
    sim.tstop = duration
    Err = 0.
    for target in result:
        Err += ((target_val['v_rest'][target] - result[target])/target_range['v_rest'][target])**2.
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: soma v_rest: %.1f, trunk v_rest: %.1f, tuft v_rest: %.1f, tuft offset: %.1f' % (os.getpid(),
                                                        result['soma'], trunk_rest, tuft_rest, result['tuft_offset'])
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if plot:
        sim.plot()
    else:
        return Err


def na_ka_error(x, plot=0):
    """

    :param x: array [soma.sh_nas, soma.gkabar, dend.gkabar factor, soma.gkdrbar, ais.gbar_nax factor, axon.kap factor,
                    axon.kdr factor]
    :param plot: int
    :return: float
    """
    if not check_bounds(x, 'na_ka'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9
    start_time = time.time()
    update_na_ka(x)
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
    duration = equilibrate + 200.
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    amp = 0.05
    while not spike:
        sim.modify_stim(0, amp=amp)
        sim.run(v_init)
        vm = np.interp(t, sim.tvec, sim.rec_list[0]['vec'])
        if amp == 0.05 and np.any(vm[:int(equilibrate/dt)] > -30.):
            print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
            return 1e9
        if np.any(vm[int(equilibrate/dt):int((equilibrate+50.)/dt)] > -30.):
            spike = True
        else:
            amp += 0.05
            if sim.verbose:
                print 'increasing amp to %.3f' % amp
    peak, threshold, ADP, AHP = get_spike_shape(vm)
    result = {'th_v': threshold, 'soma_peak': peak, 'ADP': ADP, 'AHP': AHP}
    trunk_vm = np.interp(t, sim.tvec, sim.rec_list[1]['vec'])
    th_x = np.where(vm[int(equilibrate/dt):] >= threshold)[0][0] + int(equilibrate/dt)
    trunk_peak = np.max(trunk_vm[th_x:th_x+int(10./dt)])
    trunk_pre = np.mean(trunk_vm[th_x-int(0.2/dt):th_x-int(0.1/dt)])
    result['trunk_amp'] = (trunk_peak - trunk_pre) / (peak - threshold)
    ais_vm = np.interp(t, sim.tvec, sim.rec_list[3]['vec'])
    ais_dvdt = np.gradient(ais_vm, [dt])
    axon_vm = np.interp(t, sim.tvec, sim.rec_list[4]['vec'])
    axon_dvdt = np.gradient(axon_vm, [dt])
    left = th_x - int(2./dt)
    right = th_x + int(2./dt)
    ais_peak = np.max(ais_dvdt[left:right])
    ais_peak_t = np.where(ais_dvdt[left:right] == ais_peak)[0][0] * dt
    axon_peak = np.max(axon_dvdt[left:right])
    axon_peak_t = np.where(axon_dvdt[left:right] == axon_peak)[0][0] * dt
    if axon_peak_t > ais_peak_t:
        result['ais_delay'] = 0.
    else:
        result['ais_delay'] = ais_peak_t - axon_peak_t
    if plot:
        sim.plot()
    slow_depo = 0.
    stability = 0.
    for new_amp in [0.5, 0.75, 1.0]:
        print 'increasing amp to %.3f' % new_amp
        sim.modify_stim(0, amp=new_amp)
        sim.run(v_init)
        if plot:
            sim.plot()
        vm = np.interp(t, sim.tvec, sim.rec_list[0]['vec'])
        v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
        v_before = np.max(vm[int((equilibrate - 50.)/dt):int((equilibrate - 1.)/dt)])
        v_after = np.max(vm[-int(50./dt):-1])
        stability += abs((v_before - v_rest) + (v_after - v_rest))
        if stability > 15.:
            stability += 100.
            break
        v_min_late = np.min(vm[int((equilibrate + 80.)/dt):int((equilibrate + 90.)/dt)])
        this_slow_depo = v_min_late - threshold
        if this_slow_depo > 15.:
            slow_depo += 100.
            break
        elif this_slow_depo > 0.:
            slow_depo += this_slow_depo
    result['stability'] = stability
    result['slow_depo'] = slow_depo
    Err = 0.
    for target in result:
        # don't penalize AHP or slow_depo less than target
        if not ((target == 'AHP' and result[target] < target_val['na_ka'][target]) or
                (target == 'slow_depo' and result[target] < target_val['na_ka'][target])):
            Err += ((target_val['na_ka'][target] - result[target])/target_range['na_ka'][target])**2.
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [soma.sh_nas, soma.gkabar, dend.gkabar factor, soma.gkdrbar, ais.gbar_nax factor, ' \
          'axon.kap factor, axon.kdr factor]: [%.1f, %.3E, %.1f, %.3E, %.1f, %.1f, %.1f], amp: %.3f, threshold: %.1f, ' \
          'soma_peak: %.1f, trunk_amp: %.2f, ADP: %.1f, AHP: %.1f, stability: %.1f, ais_delay: %.2f, slow_depo: %.1f, ' \
          'v_rest: %.2f' % (os.getpid(), x[0], x[1], x[2], x[3], x[4], x[5], x[6], amp, threshold, peak,
            result['trunk_amp'], ADP, AHP, result['stability'], result['ais_delay'], result['slow_depo'], v_rest)
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if not plot:
        return Err


def optimize_polish((param_name, x), maxfev=None):
    """

    :param param_name: str
    :param x: array
    :param xmin: array
    :param xmax: array
    :param maxfev: int
    :return: array
    """
    if maxfev is None:
        maxfev = 200
    error_functions = {}
    error_functions['pas'] = pas_error
    error_functions['h'] = h_error
    error_functions['v_rest'] = v_rest_error
    error_functions['na_ka'] = na_ka_error

    result = optimize.minimize(error_functions[param_name], x, method='Nelder-Mead', options={'ftol': 1e-3,
                                                    'xtol': 1e-3, 'disp': True, 'maxiter': maxfev})
    formatted_x = '['+', '.join(['%.2E' % xi for xi in result.x])+']'
    print 'Process: %i optimized %s for %i iterations with Error: %.2f and x: %s' % (os.getpid(), param_name,
                                                                            result.nit, result.fun, formatted_x)
    return {param_name: {'x': result.x, 'Err': result.fun}}


def optimize_explore((param_name, x, xmin, xmax), maxfev=None):
    """

    :param param_name: str
    :param x: array
    :param xmin: array
    :param xmax: array
    :param maxfev: int
    :return: array
    """
    if maxfev is None:
        maxfev = 200
    error_functions = {}
    error_functions['pas'] = pas_error
    error_functions['h'] = h_error
    error_functions['v_rest'] = v_rest_error
    error_functions['na_ka'] = na_ka_error

    take_step = Normalized_Step(x, xmin, xmax)
    minimizer_kwargs = dict(method=null_minimizer)
    result = optimize.basinhopping(error_functions[param_name], x, niter=maxfev, niter_success=maxfev/2,
                                       disp=True, interval=20, minimizer_kwargs=minimizer_kwargs, take_step=take_step)
    formatted_x = '['+', '.join(['%.2E' % xi for xi in result.x])+']'
    print 'Process: %i optimized %s for %i iterations with Error: %.2f and x: %s' % (os.getpid(), param_name,
                                                                            result.nit, result.fun, formatted_x)
    return {param_name: {'x': result.x, 'Err': result.fun}}


equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.01
amp = 0.3
th_dvdt = 20.
v_init = -67.
spines = True  # False
i_holding = {'soma': 0., 'trunk': 0., 'tuft': 0.}

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=spines)

# look for a trunk bifurcation
trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
if trunk_bifurcation:
    trunk_branches = [branch for branch in trunk_bifurcation[0].children if branch.type == 'trunk']
    # get where the thickest trunk branch gives rise to the tuft
    trunk = max(trunk_branches, key=lambda node: node.sec(0.).diam)
    trunk = (node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type
                                                                            for child in node.children)).next()
else:
    trunk_bifurcation = [node for node in cell.trunk if 'tuft' in (child.type for child in node.children)]
    trunk = trunk_bifurcation[0]
tuft = (child for child in trunk.children if child.type == 'tuft').next()
distal_trunk = trunk

# look for a trunk location ~50 uM from SLM
distance = cell.get_distance_to_node(trunk, tuft, 0.)
while distance < 25.:
    trunk = trunk.parent
    distance = cell.get_distance_to_node(trunk, tuft, 0.)
# look for a tuft location ~50 uM from SR
distance = cell.get_distance_to_node(distal_trunk, tuft, 0.)
while distance < 25.:
    tuft_branches = [branch for branch in tuft.children if branch.type == 'tuft' and not cell.is_terminal(branch)]
    # get the thickest non-terminal tuft branch
    tuft = max(tuft_branches, key=lambda node: node.sec(0.).diam)
    distance = cell.get_distance_to_node(distal_trunk, tuft, 0.)

"""
cell.modify_mech_param('soma', 'ions', 'ek', -77.)
cell.reinit_mechanisms()
"""

cell.modify_mech_param('soma', 'nas', 'ar', 1.)
for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
    cell.modify_mech_param(sec_type, 'nas', 'ar', 0.8)
cell.modify_mech_param('soma', 'pas', 'e', v_init)
for sec_type in ['spine_neck', 'spine_head']:
    cell.modify_mech_param(sec_type, 'pas', 'e', origin='soma')

soma_na_gbar = 0.04
cell.modify_mech_param('soma', 'nas', 'gbar', soma_na_gbar)
cell.modify_mech_param('axon_hill', 'nax', 'gbar', 2. * soma_na_gbar)
cell.modify_mech_param('axon', 'nax', 'gbar', origin='axon_hill')
cell.modify_mech_param('ais', 'nax', 'gbar', 10. * soma_na_gbar)
cell.modify_mech_param('ais', 'nax', 'sha', -5.)
for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
    cell.modify_mech_param(sec_type, 'nas', 'sha', 5.)

cell.reinit_mechanisms()

sim = QuickSim(duration)  # , verbose=False)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)
sim.append_rec(cell, cell.tree.root, loc=0.)
sim.append_rec(cell, trunk, loc=1.)
sim.append_rec(cell, tuft, loc=0.)

axon_seg_locs = np.arange(30., 500., 30.) / 500.

sim.append_rec(cell, cell.axon[1], loc=1., description='ais')
sim.append_rec(cell, cell.axon[2], loc=axon_seg_locs[0], description='axon')

#the target values and acceptable ranges
target_val = {}
target_range = {}
target_val['pas'] = {'soma': 110., 'trunk': 97., 'tuft': 121.}
target_range['pas'] = {'soma': 1., 'trunk': 1., 'tuft': 1.}
target_val['h'] = {'soma': 66., 'trunk': 39., 'tuft': 42.}  # , 'v_rest_offset': 0.}
target_range['h'] = {'soma': 0.5, 'trunk': 1., 'tuft': 1.}  # , 'v_rest_offset': 1.}
target_val['v_rest'] = {'soma': v_init, 'tuft_offset': 0.}
target_range['v_rest'] = {'soma': 0.5, 'tuft_offset': 0.1}
target_val['na_ka'] = {'th_v': -52., 'soma_peak': 40., 'trunk_amp': 0.6, 'ADP': 0., 'AHP': 4., 'stability': 0.,
                       'ais_delay': 0., 'slow_depo': 20.}
target_range['na_ka'] = {'th_v': .2, 'soma_peak': 2., 'trunk_amp': 0.01, 'ADP': 0.01, 'AHP': .2, 'stability': 1.,
                         'ais_delay': 0.001, 'slow_depo': 1.}

x0 = {}
x1 = {}
xmin = {}
xmax = {}

#x0['pas'] = [1.56E-05, 3.52E-06, 71.2]
#x0['pas'] = [1.43E-05, 1.46E-06, 55.0]
# 100615:
#x0['pas'] = [1.097e-05, 2.892e-08, 30.2]
#x0['pas'] = [1.65E-05, 1.57E-08, 27.6]
#x0['pas'] = [1.717e-05, 2.901e-07, 41.0]
x0['pas'] = [1.23E-05, 3.79E-07, 42.8]
xmin['pas'] = [1.0E-09, 1.0E-09, 25.]
xmax['pas'] = [None, None, 400.]

#x0['h'] = [1.01E-08, 1.04E-02, 50.1, 361.5]
#x0['h'] = [3.098e-08, 7.467e-03, 63.9, 290.2]
# 100615:
#x0['h'] = [2.158e-09, 4.5e-02, 50.1, 359.7]
x0['h'] = [2.24E-09, 4.51E-02, 48.2, 365.0]
xmin['h'] = [1.0E-10, 1.e-5, 25., 50.]
xmax['h'] = [1.0E-7, 5.e-2, 400., 400.]

# (soma.sh_nas, soma.gkabar, dend.gkabar factor, soma.gkdrbar, ais.gbar_nax factor, axon.kap factor, axon.kdr factor)
# x0['na_ka'] = [6.4, 3.549E-02, 2., 1.181E-02, 10., 2.5]
x0['na_ka'] = [6.2, 2.605E-02, 6.8, 5.927E-03, 13.9, 1.2, 1.2]
# threshold: -51.9, soma_peak: 36.2, trunk_amp: 0.59, ADP: 0.0, AHP: 3.6, stability: 0.1, ais_delay: 0.00,
# slow_depo: 100.0, v_rest: -66.83
xmin['na_ka'] = [1., 0.005, 2., 0.005, 1., 1., 1.]
xmax['na_ka'] = [10., 0.05, 10., 0.05, 15., 10., 10.]

x1 = dict(x0)
log = []
niter = 200

update_pas(x1['pas'])
update_h(x1['h'])
update_na_ka(x1['na_ka'])

result = optimize_explore(('na_ka', x1['na_ka'], xmin['na_ka'], xmax['na_ka']), 600)
log.append(result)
x1['na_ka'] = result['na_ka']['x']

"""
result = optimize_polish(('h', x1['h']), niter)
log.append(result)
x1['h'] = result['h']['x']
update_h(x1['h'])

result = optimize_explore(('na_ka', x1['na_ka'], xmin['na_ka'], xmax['na_ka']), 400)
log.append(result)
x1['na_ka'] = result['na_ka']['x']
update_na_ka(x1['na_ka'])
"""

"""
result = optimize_explore(('pas', x1['pas'], xmin['pas'], xmax['pas']), niter)
log.append(result)
x1['pas'] = result['pas']['x']

update_pas(x1['pas'])
update_h(x1['h'])
update_na_ka(x1['na_ka'])
result = optimize_explore(('h', x1['h'], xmin['h'], xmax['h']), niter)
log.append(result)
x1['h'] = result['h']['x']
update_h(x1['h'])
#result = optimize_explore(('na_ka', x1['na_ka'], xmin['na_ka'], xmax['na_ka']), niter)
result = optimize_polish(('na_ka', x1['na_ka']), niter)
log.append(result)
x1['na_ka'] = result['na_ka']['x']
update_na_ka(x1['na_ka'])
"""
"""
result = optimize_polish(('pas', x1['pas']), niter)
log.append(result)
x1['pas'] = result['pas']['x']
update_pas(x1['pas'])
update_h(x1['h'])
update_na_ka(x1['na_ka'])

result = optimize_explore(('h', x1['h'], xmin['h'], xmax['h']), niter)
log.append(result)
x1['h'] = result['h']['x']
update_h(x1['h'])
update_v_rest(x1['h'])
result = optimize_explore(('na_ka', x1['na_ka'], xmin['na_ka'], xmax['na_ka']), niter)
log.append(result)
x1['na_ka'] = result['na_ka']['x']
update_na_ka(x1['na_ka'])
"""