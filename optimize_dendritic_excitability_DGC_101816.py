__author__ = 'milsteina'
from specify_cells2 import *
from plot_results import *
import scipy.optimize as optimize
import os
import random
"""

Aims for spike initiation at initial segment by increasing nax density and decreasing activation V1/2 relative to soma,
axon_hill, and axon compartments. Extend linear kap gradient into basals and obliques, aim for 60% spike attenuation
at bifurcation of trunk and tuft.

Hierarchical optimization:
I) optimize g_pas for target rinp at soma, trunk bifurcation, and tuft bifurcation [without h].
II) optimize ghbar_h for target rinp at soma, trunk bifurcation, and tuft bifurcation, while also optimizing for v_rest
offset between soma and tuft, and EPSP shape changes between proximal and distal synapses measured at the soma.
III) optimize gbar_nax/nas/sh/sha, gkabar_kap/d, gkdrbar for target na spike threshold, AHP amp, and vm stability

Goal is to acheive R_inp = 75, and lower f-I gain.
"""

# morph_filename = 'EB2-late-bifurcation.swc'
# morph_filename = 'gulyas_pc2b_labeled.swc'
morph_filename = 'DGC_356550_edited.swc'
# mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'


class best_param_tracker(object):
    """

    """
    def __init__(self):
        self.x = []
        self.Err = 1e9

    def report(self):
        formatted_x = '[' + ', '.join(['%.2E' % xi for xi in self.x]) + ']'
        print 'Error: %.2f and x: %s' % (self.Err, formatted_x)


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
    soma_e_pas = cell.mech_dict['soma']['pas']['e']['value']
    distance = cell.get_distance_to_node(cell.tree.root, distal_trunk, 1.)
    slope = (current_tuft_e_pas - soma_e_pas) / (np.exp(distance / x[2]) - 1.)
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.modify_mech_param('trunk', 'pas', 'g', origin='soma', slope=x[1], tau=x[2])
    cell.modify_mech_param('trunk', 'pas', 'e', origin='soma', slope=slope, tau=x[2])
    for sec_type in ['apical', 'tuft']:
        cell.modify_mech_param(sec_type, 'pas', 'e', origin='trunk')
    for sec_type in ['basal', 'axon_hill', 'axon', 'ais', 'trunk', 'apical', 'tuft', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')


def update_h(x):
    """

    :param x: array [soma.ghbar_h, trunk.ghbar_h slope, trunk.ghbar_h tau, trunk.ghbar_h xhalf]
    """
    cell.modify_mech_param('soma', 'h', 'ghbar', x[0])
    cell.modify_mech_param('trunk', 'h', 'ghbar', origin='soma', slope=x[1], tau=x[2], xhalf=x[3])
    cell.modify_mech_param('basal', 'h', 'ghbar', origin='soma')
    for sec_type in ['apical', 'tuft']:
        cell.modify_mech_param(sec_type, 'h', 'ghbar', origin='trunk')


def update_v_rest(x):
    """

    :param x: array [soma.e_pas, tuft.e_pas]
    """
    distance = cell.get_distance_to_node(cell.tree.root, distal_trunk, 1.)
    tau = cell.mech_dict['trunk']['pas']['g']['tau']
    slope = (x[1] - x[0]) / (np.exp(distance / tau) - 1.)
    cell.modify_mech_param('soma', 'pas', 'e', x[0])
    cell.modify_mech_param('trunk', 'pas', 'e', origin='soma', slope=slope, tau=tau)
    for sec_type in ['apical', 'tuft']:
        cell.modify_mech_param(sec_type, 'pas', 'e', origin='trunk')


def update_sh_na(x):
    """

    :param x: array [ais.sh_nas]
    """
    cell.modify_mech_param('soma', 'nas', 'sh', x[0])
    for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
        cell.modify_mech_param(sec_type, 'nas', 'sh', origin='soma')
    cell.modify_mech_param('axon_hill', 'nax', 'sh', x[0])
    for sec_type in ['ais', 'axon']:
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='axon_hill')


def update_ais_delay(x):
    """

    :param x: array [ais.sha_nas, ais.gbar_nax factor]
    """
    cell.modify_mech_param('ais', 'nax', 'sha', x[0])
    axon_gbar_nax = cell.axon[2].sec.gbar_nax
    cell.modify_mech_param('ais', 'nax', 'gbar', axon_gbar_nax * x[1])


def update_na_ka_stability(x):
    """

    :param x: array [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor]
    """
    global original_ka_trunk_factor
    global original_soma_kap
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[1])
    cell.modify_mech_param('soma', 'kap', 'gkabar', x[0])
    new_ka_trunk_factor = original_soma_kap * original_ka_trunk_factor / x[0]
    new_ais_nax_factor = original_axon_nax_factor * original_ais_nax_factor / x[3]
    slope = (new_ka_trunk_factor - 1.) * x[0] / 300.
    cell.modify_mech_param('soma', 'nas', 'gbar', soma_na_gbar)
    for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
        cell.reinitialize_subset_mechanisms(sec_type, 'nas')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=x[0]+slope*75., replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300., value=x[0]+slope*300.,
                               replace=False)
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
    cell.set_terminal_branch_na_gradient()
    cell.modify_mech_param('axon_hill', 'nax', 'gbar', soma_na_gbar)
    cell.modify_mech_param('ais', 'nax', 'gbar', soma_na_gbar * x[3] * new_ais_nax_factor)
    cell.modify_mech_param('axon', 'nax', 'gbar', soma_na_gbar * x[3])
    cell.modify_mech_param('axon_hill', 'kap', 'gkabar', origin='soma')
    cell.modify_mech_param('axon_hill', 'kdr', 'gkdrbar', origin='soma')
    for sec_type in ['ais', 'axon']:
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', x[0] * x[2])


def print_new_factors():
    """
    Optimization of na_ka_stability change the values of current_ka_trunk_factor and
    current_ais_nax_factor necessary to keep max_ka_trunk and max_ais_nax constant.
    """
    new_ka_trunk_factor = original_soma_kap * original_ka_trunk_factor / \
                          cell.mech_dict['soma']['kap']['gkabar']['value']
    new_ais_nax_factor = cell.mech_dict['ais']['nax']['gbar']['value'] / \
                         cell.mech_dict['axon']['nax']['gbar']['value']
    print 'new_ka_trunk_factor: %.2E, new_ais_nax_factor: %.2E' % (new_ka_trunk_factor, new_ais_nax_factor)


def update_na_ka_dend(x):
    """

    :param x: array [soma.sh_nas, trunk.ka factor]
    """
    cell.modify_mech_param('soma', 'nas', 'sh', x[0])
    soma_gkabar = cell.mech_dict['soma']['kap']['gkabar']['value']
    slope = (x[1] - 1.) * soma_gkabar / 300.
    for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=soma_gkabar+slope*75., replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300., value=soma_gkabar+slope*300.,
                               replace=False)
        cell.modify_mech_param(sec_type, 'nas', 'sh', origin='soma')
    cell.modify_mech_param('axon_hill', 'nax', 'sh', x[0])
    for sec_type in ['ais', 'axon']:
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='axon_hill')


def offset_vm(sec_type, vm_target=None):
    """

    :param sec_type: str
    """
    if vm_target is None:
        vm_target = v_init
    sim.modify_stim(0, amp=0.)
    if sec_type == 'soma':
        sim.modify_stim(1, node=cell.tree.root, loc=0., amp=0.)
        rec = sim.rec_list[0]['vec']
    elif sec_type == 'trunk':
        sim.modify_stim(1, node=trunk, loc=1., amp=0.)
        rec = sim.rec_list[1]['vec']
    elif sec_type == 'tuft':
        sim.modify_stim(1, node=tuft, loc=1., amp=0.)
        rec = sim.rec_list[2]['vec']
    offset = True
    sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    sim.modify_stim(1, amp=i_holding[sec_type])
    sim.run(vm_target)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < vm_target - 0.5:
        i_holding[sec_type] += 0.01
        while offset:
            if sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (i_holding[sec_type], sec_type)
            sim.modify_stim(1, amp=i_holding[sec_type])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest < vm_target - 0.5:
                i_holding[sec_type] += 0.01
            else:
                offset = False
    elif v_rest > vm_target + 0.5:
        i_holding[sec_type] -= 0.01
        while offset:
            if sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (i_holding[sec_type], sec_type)
            sim.modify_stim(1, amp=i_holding[sec_type])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                i_holding[sec_type] -= 0.01
            else:
                offset = False
    sim.tstop = duration
    return initial_v_rest


def get_spike_shape(vm):
    """

    :param t: array
    :param vm: array
    :return: tuple of float: (v_peak, th_v, ADP, AHP)
    """
    vm = vm[int((equilibrate+0.4)/dt):]
    dvdt = np.gradient(vm, [dt])
    th_x = np.where(dvdt > th_dvdt)[0]
    if th_x.any():
        th_x = th_x[0] - int(1.6/dt)
    else:
        th_x = np.where(vm > -30.)[0][0] - int(2./dt)
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


def get_EPSP_dur(t, vm, equilibrate, duration):
    """

    :param t:
    :param vm:
    :param equilibrate:
    :param duration:
    :return:
    """
    dt = 0.001
    interp_t = np.arange(0., duration, dt)
    interp_vm = np.interp(interp_t, t, vm)
    left = int((equilibrate-3.)/dt)
    right = int((equilibrate-1.)/dt)
    baseline = np.mean(interp_vm[left:right])
    start = int(equilibrate/dt)
    end = int(duration/dt)
    interp_t = interp_t[start:end]
    interp_vm = interp_vm[start:end] - baseline
    amp = np.max(interp_vm)
    t_peak = np.where(interp_vm == amp)[0][0]
    interp_vm /= amp
    interp_t -= interp_t[0]
    rise_50 = np.where(interp_vm[:t_peak] >= 0.5)[0][0]
    decay_50 = np.where(interp_vm[t_peak:] <= 0.5)[0][0]
    decay_dur = interp_t[decay_50] + interp_t[t_peak] - interp_t[rise_50]
    return decay_dur


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
    cell.zero_h()
    cell.zero_na()
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
    sim.modify_stim(0, node=tuft, loc=1., amp=amp)
    sim.run(v_init)
    result['tuft'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[2]['vec']), equilibrate, duration, amp)[2]
    sim.modify_stim(1, amp=0.)
    Err = 0.
    for target in result:
        Err += ((target_val['pas'][target] - result[target])/target_range['pas'][target])**2.
    print('Simulation took %.3f s' % (time.time()-start_time))
    print 'Process %i: [soma g_pas, trunk slope, trunk tau]: [%.2E, %.2E, %.1f], soma R_inp: %.1f, ' \
          'trunk R_inp: %.1f, tuft R_inp: %.1f' % (os.getpid(), x[0], x[1], x[2], result['soma'],
                                                   result['trunk'], result['tuft'])
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if plot:
        sim.plot()
    if Err < best.Err:
        best.Err = float(Err)
        best.x = list(x)
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
    cell.zero_na()
    update_h(x)
    soma_vm = offset_vm('soma')
    result['v_rest'] = soma_vm
    sim.modify_stim(0, node=cell.tree.root, loc=0., amp=amp, dur=stim_dur)
    sim.run(v_init)
    result['soma'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[2]
    offset_vm('trunk')
    sim.modify_stim(0, node=trunk, loc=1., amp=amp)
    sim.run(v_init)
    result['trunk'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[1]['vec']), equilibrate, duration, amp)[2]
    sim.modify_stim(0, amp=0.)
    syn = trunk.spines[-1].synapses[0]
    syn.source.play(h.Vector([equilibrate]))
    sim.run(v_init)
    trunk_EPSP_dur = get_EPSP_dur(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration)
    syn.source.play(h.Vector())
    tuft_vm = offset_vm('tuft')
    sim.modify_stim(0, node=tuft, loc=1., amp=amp)
    sim.run(v_init)
    tuft_rinp = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[2]['vec']), equilibrate, duration, amp)[2]
    result['v_rest_offset'] = tuft_vm - soma_vm
    if plot:
        sim.plot()
    sim.modify_stim(0, amp=0.)
    syn = tuft.spines[-1].synapses[0]
    syn.source.play(h.Vector([equilibrate]))
    sim.run(v_init)
    tuft_EPSP_dur = get_EPSP_dur(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration)
    syn.source.play(h.Vector())
    result['delta_EPSP_dur'] = max(0., trunk_EPSP_dur - tuft_EPSP_dur)
    sim.modify_stim(1, amp=0.)
    Err = 0.
    for target in result:
        Err += ((target_val['h'][target] - result[target])/target_range['h'][target])**2.
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [soma.ghbar, trunk ghbar slope, trunk ghbar tau, trunk ghbar xhalf]: [%.2E, %.2E, %.1f, %.1f],' \
          ' soma R_inp: %.1f, trunk R_inp: %.1f, tuft R_inp: %.1f, v_rest: %.1f, v_rest offset: %.1f, ' \
          'EPSP_dur delta: %.2f' % (os.getpid(), x[0], x[1], x[2], x[3], result['soma'], result['trunk'], tuft_rinp,
                                    soma_vm, result['v_rest_offset'], result['delta_EPSP_dur'])
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if Err < best.Err:
        best.Err = float(Err)
        best.x = list(x)
    return Err


def v_rest_check(plot=0):
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
    t = np.arange(0., equilibrate, dt)
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


def v_rest_error(x, plot=0):
    """

    :param x: array (soma.e_pas, tuft.e_pas)
    :param plot: int
    :return: float
    """
    if x[1] > x[0] or not check_bounds(x, 'v_rest'):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    start_time = time.time()
    update_v_rest(x)
    sim.tstop = equilibrate
    result = {}
    sim.modify_stim(0, amp=0.)
    sim.modify_stim(1, amp=0.)
    sim.run(v_init)
    t = np.arange(0., equilibrate, dt)
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
    print 'Process %i: [soma.e_pas, tuft.e_pas]: [%.1f, %.1f], soma v_rest: %.1f, trunk v_rest: %.1f, tuft v_rest: ' \
          '%.1f, tuft offset: %.1f' % (os.getpid(), x[0], x[1], result['soma'], trunk_rest, tuft_rest,
                                       result['tuft_offset'])
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if plot:
        sim.plot()
    return Err


def na_ka_threshold_error(x, plot=0):
    """

    :param x: array [ais.sha_nas, ais.gbar_nax factor, soma.sh_nas, trunk.ka factor]
    :param plot: int
    :return: float
    """
    if not check_bounds(x[:2], 'ais_delay'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9
    if not check_bounds(x[2:], 'na_ka_dend'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9

    start_time = time.time()
    update_ais_delay(x[:2])
    update_na_ka_dend(x[2:])

    offset_vm('soma', v_active)
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
    sim.modify_rec(1, node=distal_trunk, loc=1.)
    duration = equilibrate + 200.
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    amp = 0.05
    while not spike:
        sim.modify_stim(0, amp=amp)
        sim.run(v_active)
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
    result = {}
    result['soma_peak'] = peak
    result['th_v'] = threshold
    trunk_vm = np.interp(t, sim.tvec, sim.rec_list[1]['vec'])
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
    if axon_peak_t > ais_peak_t + dt:
        result['ais_delay'] = 0.
    else:
        result['ais_delay'] = ais_peak_t + dt - axon_peak_t
    trunk_peak = np.max(trunk_vm[th_x:th_x + int(10. / dt)])
    trunk_pre = np.mean(trunk_vm[th_x - int(0.2 / dt):th_x - int(0.1 / dt)])
    result['trunk_amp'] = (trunk_peak - trunk_pre) / (peak - threshold)
    if plot:
        sim.plot()
    Err = 0.
    for target in result:
        Err += ((target_val['na_ka'][target] - result[target])/target_range['na_ka'][target])**2.
    # attempt to find the minimal combination that produces the desired delay
    Err += abs(x[0]) + x[1]
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [ais.sha_nas, ais.gbar_nas factor, sh_nax/s, trunk.ka factor]: [%.2f, %.2f, %.2f, %.2f], ' \
          'ais_delay: %.3E, soma_peak: %.1f, threshold: %.1f, trunk_amp: %.1f' % (os.getpid(), x[0], x[1], x[2], x[3],
                                                                                  result['ais_delay'], peak, threshold,
                                                                                  result['trunk_amp'])
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    sim.modify_rec(1, node=trunk, loc=1.)
    if Err < best.Err:
        best.Err = float(Err)
        best.x = list(x)
    return Err


def na_ka_stability_error(x, plot=0):
    """

    :param x: array [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor]
    :param plot: int
    :return: float
    """
    if not check_bounds(x, 'na_ka_stability'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9
    start_time = time.time()
    update_na_ka_stability(x)
    soma_vm = offset_vm('soma', v_active)
    result = {'v_rest': soma_vm}
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
    duration = equilibrate + 200.
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    amp = 0.05
    while not spike:
        sim.modify_stim(0, amp=amp)
        sim.run(v_active)
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
    if plot:
        sim.plot()
    result['ADP'] = ADP
    result['AHP'] = AHP
    slow_depo = 0.
    stability = 0.
    for new_amp in [0.25, 0.5, 0.75]:
        print 'increasing amp to %.3f' % (new_amp + amp)
        sim.modify_stim(0, amp=(new_amp + amp))
        sim.run(v_active)
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
        v_min_late = np.min(vm[int((equilibrate + 80.)/dt):int((equilibrate + 99.)/dt)])
        this_slow_depo = v_min_late - threshold
        if new_amp != 0.75 and v_min_late >= -35.:
            slow_depo += this_slow_depo + (0.75 - new_amp) * 1000.
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
    print 'Process %i: [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor]: ' \
          '[%.4f, %.4f, %.2f, %.2f], amp: %.3f, v_rest: %.1f, threshold: %.1f, ADP: %.1f, AHP: %.1f, ' \
          'stability: %.2f, slow_depo: %.2f' % (os.getpid(), x[0], x[1], x[2], x[3], amp, soma_vm,
                                                threshold, ADP, AHP, result['stability'], result['slow_depo'])
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if Err < best.Err:
        best.Err = float(Err)
        best.x = list(x)
    return Err


def ais_delay_error(x, plot=0):
    """
    Part of a set of error functions designed to optimize one parameter at a time by changing as few variables as
    possible. This function just changes the V1/2 of activation of the axon initial segment and the relative density of
    the nax sodium channel mechanism in the axon initial segment to accomplish initiation of the action potential in the
    initial segment.
    :param x: array [ais.sha_nas, ais.gbar_nax factor]
    :param plot: int
    :return: float
    """
    if x[0] > 0. or x[0] < -5. or x[1] < 1. or x[1] > 5.:
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9
    start_time = time.time()
    update_ais_delay(x)
    offset_vm('soma', v_active)
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
    duration = equilibrate + 200.
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    amp = 0.05
    while not spike:
        sim.modify_stim(0, amp=amp)
        sim.run(v_active)
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
    result = {}
    result['soma_peak'] = peak
    result['th_v'] = threshold
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
    if axon_peak_t > ais_peak_t + dt:
        result['ais_delay'] = 0.
    else:
        result['ais_delay'] = ais_peak_t + dt - axon_peak_t
    if plot:
        sim.plot()
    Err = 0.
    for target in result:
        Err += ((target_val['na_ka'][target] - result[target])/target_range['na_ka'][target])**2.
    # attempt to find the minimal combination that produces the desired delay
    Err += abs(x[0]) + x[1]
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [ais.sha_nas, ais.gbar_nas factor]: [%.2f, %.2f], ais_delay: %.3E, soma_peak: %.1f, ' \
          'threshold: %.1f' % (os.getpid(), x[0], x[1], result['ais_delay'], peak, threshold)
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if Err < best.Err:
        best.Err = float(Err)
        best.x = list(x)
    return Err


def na_ka_dend_error(x, plot=0):
    """

    :param x: array [soma.sh_nas, trunk.ka factor]
    :param plot: int
    :return: float
    """
    if not check_bounds(x, 'na_ka_dend'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9
    start_time = time.time()
    update_na_ka_dend(x)
    offset_vm('soma', v_active)
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
    sim.modify_rec(1, node=distal_trunk, loc=1.)
    duration = equilibrate + 200.
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    amp = 0.05
    while not spike:
        sim.modify_stim(0, amp=amp)
        sim.run(v_active)
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
    print 'Process %i: [sh_nax/s, trunk.ka factor]: [%.2f, %.2f], threshold: %.1f, trunk_amp: %.1f' % \
          (os.getpid(), x[0], x[1], threshold, result['trunk_amp'])
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    sim.modify_rec(1, node=trunk, loc=1.)
    if Err < best.Err:
        best.Err = float(Err)
        best.x = list(x)
    return Err


def optimize_ais_delay(x):
    """
    Simplex does not perform well for finding a good parameter set for ais initiation of APs. This is an ad hoc method
    to find a ballpark to feed into simplex.
    :param x: array
    :return: array
    """
    perturb = (-0.1, 0.2)
    ais_delay_error(x)
    while best.Err > 0.:
        index = random.randint(0, 1)
        x[index] += perturb[index]
        start_time = time.time()
        update_ais_delay(x)
        offset_vm('soma', v_active)
        sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
        duration = equilibrate + 200.
        sim.tstop = duration
        t = np.arange(0., duration, dt)
        spike = False
        amp = 0.05
        while not spike:
            sim.modify_stim(0, amp=amp)
            sim.run(v_active)
            vm = np.interp(t, sim.tvec, sim.rec_list[0]['vec'])
            if amp == 0.05 and np.any(vm[:int(equilibrate / dt)] > -30.):
                print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
                return 1e9
            if np.any(vm[int(equilibrate / dt):int((equilibrate + 50.) / dt)] > -30.):
                spike = True
            else:
                amp += 0.05
                if sim.verbose:
                    print 'increasing amp to %.3f' % amp
        peak, threshold, ADP, AHP = get_spike_shape(vm)
        result = {}
        result['soma_peak'] = peak
        result['th_v'] = threshold
        th_x = np.where(vm[int(equilibrate / dt):] >= threshold)[0][0] + int(equilibrate / dt)
        ais_vm = np.interp(t, sim.tvec, sim.rec_list[3]['vec'])
        ais_dvdt = np.gradient(ais_vm, [dt])
        axon_vm = np.interp(t, sim.tvec, sim.rec_list[4]['vec'])
        axon_dvdt = np.gradient(axon_vm, [dt])
        left = th_x - int(2. / dt)
        right = th_x + int(2. / dt)
        ais_peak = np.max(ais_dvdt[left:right])
        ais_peak_t = np.where(ais_dvdt[left:right] == ais_peak)[0][0] * dt
        axon_peak = np.max(axon_dvdt[left:right])
        axon_peak_t = np.where(axon_dvdt[left:right] == axon_peak)[0][0] * dt
        if axon_peak_t > ais_peak_t + dt:
            result['ais_delay'] = 0.
        else:
            result['ais_delay'] = ais_peak_t + dt - axon_peak_t
        Err = 0.
        target = 'ais_delay'
        Err += ((target_val['na_ka'][target] - result[target]) / target_range['na_ka'][target]) ** 2.
        # attempt to find the minimal combination that produces the desired delay
        print 'Simulation took %i s' % (time.time() - start_time)
        print 'Process %i: [ais.sha_nas, ais.gbar_nas factor]: [%.2f, %.2f], ais_delay: %.3E, soma_peak: %.1f, ' \
              'threshold: %.1f' % (os.getpid(), x[0], x[1], result['ais_delay'], peak, threshold)
        print 'Process %i: Error: %.4E' % (os.getpid(), Err)
        if Err < best.Err:
            best.Err = float(Err)
            best.x = list(x)
        if Err == 0.:
            return x


def optimize_na_ka_dend(x, maxfev=100):
    """
    Simplex does not perform well for finding a good parameter set for the voltage threshold of APs. This is an ad hoc
    method to find a ballpark to feed into simplex.
    :param x: array
    :return: array
    """
    perturb = (0.2, -0.05)
    target_keys = ('th_v', 'trunk_amp')
    iter = 0
    result = {}
    # best_result = None
    while iter < maxfev:
        start_time = time.time()
        update_na_ka_dend(x)
        offset_vm('soma', v_active)
        sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
        sim.modify_rec(1, node=distal_trunk, loc=1.)
        duration = equilibrate + 200.
        sim.tstop = duration
        t = np.arange(0., duration, dt)
        spike = False
        amp = 0.05
        while not spike:
            sim.modify_stim(0, amp=amp)
            sim.run(v_active)
            vm = np.interp(t, sim.tvec, sim.rec_list[0]['vec'])
            if amp == 0.05 and np.any(vm[:int(equilibrate / dt)] > -30.):
                print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
                return 1e9
            if np.any(vm[int(equilibrate / dt):int((equilibrate + 50.) / dt)] > -30.):
                spike = True
            else:
                amp += 0.05
                if sim.verbose:
                    print 'increasing amp to %.3f' % amp
        peak, threshold, ADP, AHP = get_spike_shape(vm)
        result['th_v'] = threshold
        trunk_vm = np.interp(t, sim.tvec, sim.rec_list[1]['vec'])
        th_x = np.where(vm[int(equilibrate / dt):] >= threshold)[0][0] + int(equilibrate / dt)
        trunk_peak = np.max(trunk_vm[th_x:th_x + int(10. / dt)])
        trunk_pre = np.mean(trunk_vm[th_x - int(0.2 / dt):th_x - int(0.1 / dt)])
        result['trunk_amp'] = (trunk_peak - trunk_pre) / (peak - threshold)
        Err = 0.
        for target in result:
            Err += ((target_val['na_ka'][target] - result[target]) / target_range['na_ka'][target]) ** 2.
        print 'Simulation took %i s' % (time.time() - start_time)
        print 'Process %i: [sh_nax/s, trunk.ka factor]: [%.2f, %.2f], threshold: %.1f, trunk_amp: %.1f' % \
              (os.getpid(), x[0], x[1], threshold, result['trunk_amp'])
        print 'Process %i: Error: %.4E' % (os.getpid(), Err)
        sim.modify_rec(1, node=trunk, loc=1.)
        if Err < best.Err:
            best.Err = float(Err)
            best.x = list(x)
            # best_result = dict(result)
        if Err < 10.:
            return x
        if iter == 0:
            index = random.randint(0, 1)
            direction = random.choice([-1, 1])
            x[index] += direction * perturb[index]
        else:
            perturbed = False
            for index in range(2):
                if result[target_keys[index]] < target_val['na_ka'][target_keys[index]]:
                    x[index] += perturb[index]
                    perturbed = True
                elif result[target_keys[index]] > target_val['na_ka'][target_keys[index]]:
                    x[index] -= perturb[index]
                    perturbed = True
            if not perturbed:
                index = random.randint(0, 1)
                direction = random.choice([-1, 1])
                x[index] += direction * perturb[index]
        iter += 1


def optimize_polish(param_name, x=None, maxfev=None):
    """

    :param param_name: str
    :param x: array
    :param maxfev: int
    :return: array
    """
    if maxfev is None:
        maxfev = 200
    if x is None:
        if param_name == 'na_ka_threshold':
            x = x1['ais_delay'] + x1['na_ka_dend']
        else:
            x = x1[param_name]
    result = optimize.minimize(error_functions[param_name], x, method='Nelder-Mead', options={'ftol': 1e-3,
                                                    'xtol': 1e-3, 'disp': True, 'maxiter': maxfev})
    formatted_x = '['+', '.join(['%.2E' % xi for xi in result.x])+']'
    print 'Process: %i optimized %s for %i iterations with Error: %.2f and x: %s' % (os.getpid(), param_name,
                                                                            result.nit, result.fun, formatted_x)
    return {param_name: {'x': result.x, 'Err': result.fun}}


def optimize_explore(param_name, x=None, local_xmin=None, local_xmax=None, maxfev=None):
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
    if x is None:
        if param_name == 'na_ka_threshold':
            x = x1['ais_delay'] + x1['na_ka_dend']
        else:
            x = x1[param_name]
    if local_xmin is None:
        if param_name == 'na_ka_threshold':
            local_xmin = xmin['ais_delay'] + xmin['na_ka_dend']
        else:
            local_xmin = xmin[param_name]
    if local_xmax is None:
        if param_name == 'na_ka_threshold':
            local_xmax = xmax['ais_delay'] + xmax['na_ka_dend']
        else:
            local_xmax = xmax[param_name]
    take_step = Normalized_Step(x, local_xmin, local_xmax)
    minimizer_kwargs = dict(method=null_minimizer)
    result = optimize.basinhopping(error_functions[param_name], x, niter=maxfev, niter_success=maxfev/2,
                                       disp=True, interval=20, minimizer_kwargs=minimizer_kwargs, take_step=take_step)
    formatted_x = '['+', '.join(['%.2E' % xi for xi in result.x])+']'
    print 'Process: %i optimized %s for %i iterations with Error: %.2f and x: %s' % (os.getpid(), param_name,
                                                                            result.nit, result.fun, formatted_x)
    if param_name == 'na_ka_stability':
        print_new_factors()
    return {param_name: {'x': result.x, 'Err': result.fun}}


equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.01
amp = 0.3
th_dvdt = 10.
v_init = -67.
v_active = -61.
spines = False  # True

cell = DG_GC(morph_filename, full_spines=spines)

soma_ek = -77.
#cell.modify_mech_param('soma', 'ions', 'ek', soma_ek)

soma_na_gbar = 0.04

"""
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
tuft = distal_trunk  # testing optimizing pas and h at tuft bifurcation instead of in thin tuft

for branch in [trunk, tuft]:
    syn = Synapse(cell, branch.spines[-1], ['AMPA_KIN'], stochastic=0)
cell.init_synaptic_mechanisms()

cell.modify_mech_param('ais', 'nax', 'sha', -2.8)
for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
    cell.modify_mech_param(sec_type, 'nas', 'sha', 5.)

cell.modify_mech_param('soma', 'km2', 'gkmbar', 0.0015)
cell.modify_mech_param('axon_hill', 'km2', 'gkmbar', origin='soma')
cell.modify_mech_param('ais', 'km2', 'gkmbar', 0.0075)
cell.modify_mech_param('axon', 'km2', 'gkmbar', origin='ais')

cell.reinit_mechanisms()
cell.set_terminal_branch_na_gradient()
"""

sim = QuickSim(duration)  # , verbose=False)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)
sim.append_rec(cell, cell.tree.root, loc=0.)
#sim.append_rec(cell, trunk, loc=1.)
#sim.append_rec(cell, tuft, loc=1.)  # loc=0.)

axon_seg_locs = [seg.x for seg in cell.axon[2].sec]

sim.append_rec(cell, cell.axon[1], loc=1., description='ais')
sim.append_rec(cell, cell.axon[2], loc=axon_seg_locs[0], description='axon')

i_holding = {'soma': 0.09, 'primary': 0.09, 'secondary': 0.09}
compartment_objects = {'soma': cell.tree.root}


#the target values and acceptable ranges
target_val = {}
target_range = {}
target_val['pas'] = {'soma': 295., 'dend': 375.}
target_range['pas'] = {'soma': 10., 'dend': 10.}
target_val['h'] = {'v_rest': v_init, 'soma': 75., 'trunk': 39., 'tuft': 42., 'v_rest_offset': 0., 'delta_EPSP_dur': 0.}
target_range['h'] = {'v_rest': 0.25, 'soma': 0.5, 'trunk': 1., 'tuft': 2., 'v_rest_offset': 0.1, 'delta_EPSP_dur': 0.05}
target_val['v_rest'] = {'soma': v_init, 'tuft_offset': 0.}
target_range['v_rest'] = {'soma': 0.25, 'tuft_offset': 0.1}
target_val['na_ka'] = {'v_rest': v_init, 'th_v': -51., 'soma_peak': 40., 'trunk_amp': 0.5, 'ADP': 0., 'AHP': 4.,
                       'stability': 0., 'ais_delay': 0., 'slow_depo': 25.}
target_range['na_ka'] = {'v_rest': 0.25, 'th_v': .2, 'soma_peak': 2., 'trunk_amp': 0.01, 'ADP': 0.01, 'AHP': .2,
                         'stability': 1., 'ais_delay': 0.001, 'slow_depo': 1.}

error_functions = {'pas': pas_error, 'h': h_error, 'v_rest': v_rest_error, 'ais_delay': ais_delay_error,
                   'na_ka_stability': na_ka_stability_error, 'na_ka_dend': na_ka_dend_error,
                   'na_ka_threshold': na_ka_threshold_error}

x0 = {}
x1 = {}
xmin = {}
xmax = {}

x0['pas'] = [6.27E-06, 2.05E-09, 35.4]  # Error: 17.67
xmin['pas'] = [1.0E-6, 1.0E-06, 25.]
xmax['pas'] = [2.0E-5, 2.0E-05, 200.]

x0['h'] = [8.45E-09, 4.27E-04, 2.51E+01, 8.24E+01]  # Error: 121.50
xmin['h'] = [1.0E-9, 1.e-4, 25., 50.]
xmax['h'] = [1.0E-7, 1.e-2, 200., 500.]

# [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor]
# remember to run print_new_factors() and update x0['na_ka_dend'][1] and x0['ais_delay'][1] accordingly
x0['na_ka_stability'] = [0.0585, 0.0371, 2.28, 4.82]  # Error: 88.40
xmin['na_ka_stability'] = [0.01, 0.01, 1., 1.5]
xmax['na_ka_stability'] = [0.06, 0.06, 5., 5.]

# [soma.sh_nas, trunk.ka factor]
x0['na_ka_dend'] = [1.97, 1.21]  # Error: 13.57
xmin['na_ka_dend'] = [0.1, 1.1]
xmax['na_ka_dend'] = [4., 5.]

# [ais.sha_nas, ais.gbar_nax factor]
x0['ais_delay'] = [-1.20, 1.57]  # Error: 29.16
xmin['ais_delay'] = [-5., 1.1]
xmax['ais_delay'] = [-1., 5.]

# [soma.e_pas, tuft.e_pas]
x0['v_rest'] = [-63., -77.]
xmin['v_rest'] = [v_init, soma_ek]
xmax['v_rest'] = [-63., -63.]

x1 = dict(x0)

"""
current_tuft_e_pas = x1['v_rest'][1]
update_v_rest(x1['v_rest'])
update_pas(x1['pas'])
update_h(x1['h'])

original_soma_kap = x1['na_ka_stability'][0]
original_axon_nax_factor = x1['na_ka_stability'][3]
original_ais_nax_factor = x1['ais_delay'][1]
original_ka_trunk_factor = x1['na_ka_dend'][1]

update_na_ka_stability(x1['na_ka_stability'])
update_na_ka_dend(x1['na_ka_dend'])
update_ais_delay(x1['ais_delay'])
"""

best = best_param_tracker()