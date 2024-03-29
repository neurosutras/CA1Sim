__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
import os
"""

Aims for spike initiation at initial segment by increasing nax density and decreasing activation V1/2 relative to soma,
axon_hill, and axon compartments.

Extends linear kap gradient into basals and obliques, aims for 60% spike attenuation at trunk bifurcation.

Hierarchical optimization (no spines):
I) optimize g_pas for target rinp
II) optimize gbar_nax, gkabar_kap/d, gkdrbar_kdr for target na spike shape
"""

#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '101415 simple_axon_model'


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
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.modify_mech_param('trunk', 'pas', 'g', origin='soma', slope=x[1], tau=x[2])
    for sec_type in ['basal', 'axon_hill', 'axon', 'ais', 'trunk', 'apical', 'tuft']:
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


def update_na_ka(x):
    """

    :param x: array [soma.sh_nax, soma.gkabar, soma.gkdrbar, trunk.ka factor]
    """
    cell.modify_mech_param('soma', 'nax', 'gbar', soma_na_gbar)
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[2])
    cell.modify_mech_param('soma', 'kap', 'gkabar', x[1])
    cell.modify_mech_param('soma', 'nax', 'sh', x[0])
    slope = (x[3] - 1.) * x[1] / 300.
    for sec_type in ['basal', 'apical', 'trunk', 'tuft']:
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=x[1]+slope*75., replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300., value=x[1]+slope*300.,
                               replace=False)
        cell.modify_mech_param(sec_type, 'nax', 'gbar', origin='soma')
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='soma')
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param('axon_hill', 'nax', 'gbar', soma_na_gbar * 2.)
    cell.modify_mech_param('ais', 'nax', 'gbar', soma_na_gbar * 2 * 5.)
    cell.modify_mech_param('axon', 'nax', 'gbar', origin='axon_hill')
    for sec_type in ['axon_hill', 'ais', 'axon']:
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma')
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='soma')


def update_na_ka_soma(x):
    """

    :param x: array [soma.sh_nax, soma.gkabar, soma.gkdrbar]
    """
    cell.modify_mech_param('soma', 'nax', 'gbar', soma_na_gbar)
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[2])
    cell.modify_mech_param('soma', 'kap', 'gkabar', x[1])
    cell.modify_mech_param('soma', 'nax', 'sh', x[0])
    for sec_type in ['basal', 'apical', 'trunk', 'tuft']:
        cell.modify_mech_param(sec_type, 'nax', 'gbar', origin='soma')
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='soma')
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param('axon_hill', 'nax', 'gbar', soma_na_gbar * 2.)
    cell.modify_mech_param('ais', 'nax', 'gbar', soma_na_gbar * 2 * 5.)
    cell.modify_mech_param('axon', 'nax', 'gbar', origin='axon_hill')
    for sec_type in ['axon_hill', 'ais', 'axon']:
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma')
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='soma')


def update_na_ka_dend(x):
    """

    :param x: array [soma.sh_nax, trunk.ka factor]
    """
    cell.modify_mech_param('soma', 'nax', 'sh', x[0])
    soma_gkabar = cell.mech_dict['soma']['kap']['gkabar']['value']
    slope = (x[1] - 1.) * soma_gkabar / 300.
    for sec_type in ['basal', 'apical', 'trunk', 'tuft']:
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=soma_gkabar+slope*75., replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300., value=soma_gkabar+slope*300.,
                               replace=False)
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='soma')
    for sec_type in ['axon_hill', 'ais', 'axon']:
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='soma')


def offset_vm(sec_type):
    """

    :param sec_type: str
    """
    sim.modify_stim(0, amp=0.)
    i_holding = {'soma': 0., 'trunk': 0., 'tuft': 0.}
    if sec_type == 'soma':
        sim.modify_stim(1, node=cell.tree.root, loc=0., amp=i_holding[sec_type])
        rec = sim.rec_list[0]['vec']
    elif sec_type == 'trunk':
        sim.modify_stim(1, node=trunk, loc=1., amp=i_holding[sec_type])
        rec = sim.rec_list[1]['vec']
    elif sec_type == 'tuft':
        sim.modify_stim(1, node=tuft, loc=1., amp=i_holding[sec_type])
        rec = sim.rec_list[2]['vec']
    offset = True
    sim.run(v_init)
    t = np.arange(0, duration, dt)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    if v_rest < v_init - 0.5:
        i_holding[sec_type] += 0.005
        while offset:
            if sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (i_holding[sec_type], sec_type)
            sim.modify_stim(1, amp=i_holding[sec_type])
            sim.run(v_init)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest < v_init:
                i_holding[sec_type] += 0.025
            else:
                offset = False
    else:
        i_holding[sec_type] -= 0.025
        while offset:
            if sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (i_holding[sec_type], sec_type)
            sim.modify_stim(1, amp=i_holding[sec_type])
            sim.run(v_init)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > v_init:
                i_holding[sec_type] -= 0.025
            else:
                offset = False
    return i_holding[sec_type]


def zero_na():
    """

    """
    for sec_type in ['axon_hill', 'ais', 'soma']:
        cell.modify_mech_param(sec_type, 'nax', 'gbar', 0.)
    for sec_type in ['axon', 'basal', 'trunk', 'apical', 'tuft']:
        cell.reinitialize_subset_mechanisms(sec_type, 'nax')


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
    x_peak = np.where(vm[th_x:th_x+int(5./dt)] == v_peak)[0][0]
    v_AHP = np.min(vm[th_x+x_peak:th_x+int(12./dt)])
    x_AHP = np.where(vm[th_x+x_peak:th_x+int(12./dt)] == v_AHP)[0][0]
    AHP = v_before - v_AHP
    # if spike waveform includes an ADP before an AHP, return the value of the ADP in order to increase error function
    rising_x = np.where(dvdt[th_x+x_peak:th_x+x_peak+x_AHP] > 0.)[0]
    if rising_x.any():
        v_ADP = np.max(vm[th_x+x_peak+rising_x[0]:th_x+x_peak+x_AHP])
        ADP = v_ADP - v_AHP
    else:
        ADP = 0.
    return v_peak, th_v, ADP, AHP


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
    zero_na()
    update_pas(x)
    #offset_vm('soma')
    sim.modify_stim(1, amp=0.)
    sim.modify_stim(0, node=cell.tree.root, loc=0., amp=amp, dur=stim_dur)
    sim.run(v_init)
    result['soma'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[2]
    #offset_vm('trunk')
    sim.modify_stim(0, node=trunk, loc=1., amp=amp)
    sim.run(v_init)
    result['trunk'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[1]['vec']), equilibrate, duration, amp)[2]
    #offset_vm('tuft')
    sim.modify_stim(0, node=tuft, loc=1., amp=amp)
    sim.run(v_init)
    result['tuft'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[2]['vec']), equilibrate, duration, amp)[2]
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
    if np.any(x < 1.e-8) or x[2] < 50.:  # or x[2] + x[3] > 300.:
        return 1e9
    start_time = time.time()
    sim.tstop = duration
    result = {}
    amp = -0.15
    zero_na()
    update_h(x)
    #offset_vm('soma')
    sim.modify_stim(1, amp=0.)
    sim.modify_stim(0, node=cell.tree.root, loc=0., amp=amp, dur=stim_dur)
    sim.run(v_init)
    result['soma'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[2]
    #offset_vm('trunk')
    sim.modify_stim(0, node=trunk, loc=1., amp=amp)
    sim.run(v_init)
    result['trunk'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[1]['vec']), equilibrate, duration, amp)[2]
    #offset_vm('tuft')
    sim.modify_stim(0, node=tuft, loc=1., amp=amp)
    sim.run(v_init)
    result['tuft'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[2]['vec']), equilibrate, duration, amp)[2]
    Err = 0.
    for target in result:
        Err += ((target_val['h'][target] - result[target])/target_range['h'][target])**2.
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [soma.ghbar, trunk ghbar slope, trunk ghbar tau, trunk ghbar xhalf]: [%.2E, %.2E, %.1f, %.1f],' \
          ' soma R_inp: %.1f, trunk R_inp: %.1f, tuft R_inp: %.1f' % (os.getpid(), x[0], x[1], x[2], x[3],
                                                                     result['soma'], result['trunk'], result['tuft'])
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
    sim.tstop = duration
    result = {}
    sim.modify_stim(0, amp=0.)
    sim.modify_stim(1, amp=0.)
    sim.run(v_init)
    t = np.arange(0, duration, dt)
    vm = np.interp(t, sim.tvec, sim.rec_list[0]['vec'])
    result['soma'] = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    vm = np.interp(t, sim.tvec, sim.rec_list[1]['vec'])
    trunk_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    vm = np.interp(t, sim.tvec, sim.rec_list[2]['vec'])
    tuft_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    result['tuft_offset'] = tuft_rest - result['soma']
    Err = 0.
    for target in result:
        Err += ((target_val['v_rest'][target] - result[target])/target_range['v_rest'][target])**2.
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [soma.e_pas, tuft.e_pas offset]: [%.2f, %.2f], ' \
          'soma v_rest: %.1f, trunk v_rest: %.1f, tuft v_rest: %.1f' % (os.getpid(), cell.soma[0].sec.e_pas,
                                                        result['tuft_offset'], result['soma'], trunk_rest, tuft_rest)
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if plot:
        sim.plot()
    else:
        return Err


def na_ka_error(x, plot=0):
    """

    :param x: array [soma.sh_nax, soma.gkabar, soma.gkdrbar, trunk.ka factor]
    :param plot: int
    :return: float
    """
    if not check_bounds(x, 'na_ka'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9
    start_time = time.time()
    update_na_ka(x)
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
    sim.modify_rec(1, node=distal_trunk, loc=1.)
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
        if new_amp != 1. and this_slow_depo > 15.:
            slow_depo += this_slow_depo + 100.
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
    print 'Process %i: [soma.sh_nax, soma.gkabar, soma.gkdrbar, trunk.ka factor]: [%.1f, %.3E, %.3E, %.1f], ' \
          'amp: %.3f, threshold: %.1f, soma_peak: %.1f, trunk_amp: %.3f, ADP: %.1f, AHP: %.1f, stability: %.1f, ' \
          'ais_delay: %.2f, slow_depo: %.1f, v_rest: %.2f' % (os.getpid(), x[0], x[1], x[2], x[3], amp, threshold,
            peak, result['trunk_amp'], ADP, AHP, result['stability'], result['ais_delay'], result['slow_depo'], v_rest)
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    sim.modify_rec(1, node=trunk, loc=1.)
    if not plot:
        return Err


def na_ka_error_soma(x, plot=0):
    """

    :param x: array [soma.sh_nax, soma.gkabar, soma.gkdrbar]
    :param plot: int
    :return: float
    """
    """
    if not check_bounds(x, 'na_ka'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9
    """
    start_time = time.time()
    update_na_ka_soma(x)
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
    th_x = np.where(vm[int(equilibrate/dt):] >= threshold)[0][0] + int(equilibrate/dt)
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
        if new_amp != 1. and this_slow_depo > 15.:
            slow_depo += this_slow_depo + 100.
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
    print 'Process %i: [soma.sh_nax, soma.gkabar, soma.gkdrbar]: [%.1f, %.3E, %.3E], amp: %.3f, threshold: %.1f, ' \
          'soma_peak: %.1f, ADP: %.1f, AHP: %.1f, stability: %.1f, ais_delay: %.2f, slow_depo: %.1f, ' \
          'v_rest: %.2f' % (os.getpid(), x[0], x[1], x[2], amp, threshold, peak, ADP, AHP,
            result['stability'], result['ais_delay'], result['slow_depo'], v_rest)
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if not plot:
        return Err


def na_ka_error_dend(x, plot=0):
    """

    :param x: array [soma.sh_nax, trunk.ka factor]
    :param plot: int
    :return: float
    """
    """
    if not check_bounds(x, 'na_ka'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9
    """
    start_time = time.time()
    update_na_ka_dend(x)
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
    sim.modify_rec(1, node=distal_trunk, loc=1.)
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
    result = {'th_v': threshold}
    trunk_vm = np.interp(t, sim.tvec, sim.rec_list[1]['vec'])
    th_x = np.where(vm[int(equilibrate/dt):] >= threshold)[0][0] + int(equilibrate/dt)
    trunk_peak = np.max(trunk_vm[th_x:th_x+int(10./dt)])
    trunk_pre = np.mean(trunk_vm[th_x-int(0.2/dt):th_x-int(0.1/dt)])
    result['trunk_amp'] = (trunk_peak - trunk_pre) / (peak - threshold)
    if plot:
        sim.plot()
    Err = 0.
    for target in result:
        Err += ((target_val['na_ka'][target] - result[target])/target_range['na_ka'][target])**2.
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [soma.sh_nax, trunk.ka factor]: [%.1f, %.1f], amp: %.3f, threshold: %.1f, trunk_amp: %.3f' % \
          (os.getpid(), x[0], x[1], amp, threshold, result['trunk_amp'])
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    sim.modify_rec(1, node=trunk, loc=1.)
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
    error_functions['na_ka_soma'] = na_ka_error_soma
    error_functions['na_ka_dend'] = na_ka_error_dend

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
        maxfev = 400
    error_functions = {}
    error_functions['pas'] = pas_error
    error_functions['h'] = h_error
    error_functions['v_rest'] = v_rest_error
    error_functions['na_ka'] = na_ka_error
    error_functions['na_ka_soma'] = na_ka_error_soma
    error_functions['na_ka_dend'] = na_ka_error_dend

    take_step = Normalized_Step(x, xmin, xmax)
    minimizer_kwargs = dict(method=null_minimizer)
    result = optimize.basinhopping(error_functions[param_name], x, niter=maxfev, niter_success=maxfev/2, disp=True,
                                   interval=20, minimizer_kwargs=minimizer_kwargs, take_step=take_step)
    formatted_x = '['+', '.join(['%.2E' % xi for xi in result.x])+']'
    print 'Process: %i optimized %s for %i iterations with Error: %.2f and x: %s' % (os.getpid(), param_name,
                                                                            result.nit, result.fun, formatted_x)
    return {param_name: {'x': result.x, 'Err': result.fun}}


equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.01
#i_holding = {'soma': 0., 'trunk': 0., 'tuft': 0.}
amp = 0.3
th_dvdt = 20.
v_init = -65.
spines = True  # False

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=False)

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
trunk = trunk_bifurcation[0]

#cell.axon[1].set_diam_bounds(2, 0.5)
#cell.axon[2].sec.diam = 0.5
#cell.reinit_mechanisms(reset_cable=1)

soma_na_gbar = 0.04
cell.modify_mech_param('soma', 'pas', 'e', v_init)
for sec_type in ['basal', 'axon_hill', 'ais', 'axon', 'trunk', 'apical', 'tuft']:
    cell.modify_mech_param(sec_type, 'pas', 'e', origin='soma')
cell.modify_mech_param('ais', 'nax', 'sha', -5.)
for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
    cell.modify_mech_param(sec_type, 'nax', 'sha', 5.)
cell.modify_mech_param('soma', 'km2', 'gkmbar', 0.001)

cell.modify_mech_param('axon_hill', 'km2', 'gkmbar', 0.002)
for sec_type in ['ais', 'axon']:
    cell.modify_mech_param(sec_type, 'km2', 'gkmbar', origin='axon_hill')
"""
for sec_type in ['axon_hill', 'ais', 'axon']:
    cell.modify_mech_param(sec_type, 'km2', 'gkmbar', origin='soma')
"""

sim = QuickSim(duration)  # , verbose=False)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)
sim.append_rec(cell, cell.tree.root, loc=0.)
sim.append_rec(cell, trunk, loc=1.)
sim.append_rec(cell, tuft, loc=1.)

axon_seg_locs = [seg.x for seg in cell.axon[2].sec]

sim.append_rec(cell, cell.axon[1], loc=1., description='ais')
sim.append_rec(cell, cell.axon[2], loc=axon_seg_locs[0], description='axon')

#the target values and acceptable ranges
target_val = {}
target_range = {}
#target_val['pas'] = {'soma': 110., 'trunk': 97., 'tuft': 121.}
target_range['pas'] = {'soma': 1., 'trunk': 1., 'tuft': 1.}
#target_val['h'] = {'soma': 66., 'trunk': 39., 'tuft': 42.}
target_val['pas'] = {'soma': 66., 'trunk': 39., 'tuft': 42.}
target_val['na_ka'] = {'th_v': -52., 'soma_peak': 40., 'trunk_amp': 0.6, 'ADP': 0., 'AHP': 4., 'stability': 0.,
                       'ais_delay': 0., 'slow_depo': 30.}
target_range['na_ka'] = {'th_v': .2, 'soma_peak': 2., 'trunk_amp': 0.01, 'ADP': 0.01, 'AHP': .2, 'stability': 1.,
                         'ais_delay': 0.001, 'slow_depo': 1.}

x0 = {}
x1 = {}
xmin = {}
xmax = {}


#x0['pas'] = [5.222e-09, 7.428e-05, 7.53e+01]
#x0['pas'] = [3.813E-07, 1.101E-04, 84.8]
x0['pas'] = [4.55E-07, 6.89E-05, 73.2]
xmin['pas'] = [1.0E-09, 1.0E-09, 25.]
xmax['pas'] = [None, None, 400.]

#[soma.sh_nax, soma.gkabar, soma.gkdrbar, trunk.ka factor]
#x0['na_ka'] = [1.0, 7.000E-03, 9.800E-03, 3.]
x0['na_ka'] = [2.7, 4.367E-02, 1.282E-02, 1.5]
xmin['na_ka'] = [-1., 0.001, 0.001, 1.]
xmax['na_ka'] = [10., 0.05, 0.05, 15.]


x1 = dict(x0)
log = []

update_pas(x1['pas'])
update_na_ka(x1['na_ka'])

#update_na_ka_soma([3.1, 4.433E-02, 1.215E-02])
#update_na_ka(x1['na_ka'])
#result = optimize_explore(('na_ka', x1['na_ka'], xmin['na_ka'], xmax['na_ka']), 600)
#result = optimize_polish(('na_ka', x1['na_ka']), 400)
#log.append(result)
#x1['na_ka'] = result['na_ka']['x']
#update_na_ka(x1['na_ka'])

#result = optimize_polish(('pas', x1['pas']))
#log.append(result)
#x1['pas'] = result['pas']['x']
