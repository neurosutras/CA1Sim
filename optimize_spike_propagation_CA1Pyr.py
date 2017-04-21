__author__ = 'Aaron D. Milstein'
from specify_cells2 import *
from plot_results import *
import sys
import os
import random

morph_filename = 'EB2-late-bifurcation.swc'

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    mech_filename = '041817 CA1Pyr optimizing spike stability'


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


def get_spike_shape(vm):
    """

    :param vm: array
    :return: tuple of float: (v_peak, th_v, ADP, AHP)
    """
    vm = vm[int((equilibrate+1.)/dt):]
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
    # end = min(th_x + int(50. / dt), len(vm))
    end = len(vm)
    v_AHP = np.min(vm[th_x + x_peak:end])
    x_AHP = np.where(vm[th_x + x_peak:end] == v_AHP)[0][0]
    AHP = v_before - v_AHP
    # if spike waveform includes an ADP before an AHP, return the value of the ADP in order to increase error function
    rising_x = np.where(dvdt[th_x+x_peak:th_x+x_peak+x_AHP] > 0.)[0]
    if rising_x.any():
        v_ADP = np.max(vm[th_x+x_peak+rising_x[0]:th_x+x_peak+x_AHP])
        ADP = v_ADP - th_v
    else:
        ADP = 0.
    return v_peak, th_v, ADP, AHP


def update_ais_delay(x):
    """

    :param x: array [ais.sha_nas, ais.gbar_nax]
    """
    cell.modify_mech_param('ais', 'nax', 'sha', x[0])
    cell.modify_mech_param('ais', 'nax', 'gbar', x[1])


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
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    print 'Trying x: %s: %s' % (str(xlabels['ais_delay']), formatted_x)
    hist.x_values.append(x)
    if not check_bounds.within_bounds(x, 'ais_delay'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        Err = 1e9
        hist.error_values.append(Err)
        return Err
    start_time = time.time()
    update_ais_delay(x)
    offset_vm('soma', v_active)
    stim_dur = 150.
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=stim_dur)
    duration = equilibrate + stim_dur
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    d_amp = 0.01
    amp = i_th['soma'] - 0.02
    while not spike:
        sim.modify_stim(0, amp=amp)
        sim.run(v_active)
        vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
        if np.any(vm[:int(equilibrate/dt)] > -30.):
            print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
            return 1e9
        if np.any(vm[int(equilibrate/dt):int((equilibrate+50.)/dt)] > -30.):
            spike = True
        else:
            amp += d_amp
            if sim.verbose:
                print 'increasing amp to %.3f' % amp
    i_th['soma'] = amp
    peak, threshold, ADP, AHP = get_spike_shape(vm)
    result = {}
    result['soma_peak'] = peak
    result['th_v'] = threshold
    th_x = np.where(vm[int(equilibrate/dt):] >= threshold)[0][0] + int(equilibrate/dt)
    ais_vm = np.interp(t, sim.tvec, sim.get_rec('ais')['vec'])
    ais_dvdt = np.gradient(ais_vm, dt)
    axon_vm = np.interp(t, sim.tvec, sim.get_rec('axon')['vec'])
    axon_dvdt = np.gradient(axon_vm, dt)
    left = th_x - int(2./dt)
    right = th_x + int(5./dt)
    ais_peak = np.max(ais_dvdt[left:right])
    ais_peak_t = np.where(ais_dvdt[left:right] == ais_peak)[0][0] * dt
    axon_peak = np.max(axon_dvdt[left:right])
    axon_peak_t = np.where(axon_dvdt[left:right] == axon_peak)[0][0] * dt
    if axon_peak_t >= ais_peak_t + dt:
        result['ais_delay'] = 0.
    else:
        result['ais_delay'] = ais_peak_t + dt - axon_peak_t
    if plot:
        sim.plot()
    Err = 0.
    for target in result:
        Err += ((target_val['na_ka'][target] - result[target])/target_range['na_ka'][target])**2.
    # attempt to find the minimal combination that produces the desired delay
    for target in result:
        if target not in hist.features:
            hist.features[target] = []
        hist.features[target].append(result[target])
    for i, x_i in enumerate(x):
        Err += ((x_i - xmin['ais_delay'][i])/(0.05*(abs(xmin['ais_delay'][i]) - abs(xmax['ais_delay'][i]))))**2.
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [ais.sha_nas, ais.gbar_nas]: [%.3E, %.3E], ais_delay: %.3E, soma_peak: %.1f, ' \
          'threshold: %.1f' % (os.getpid(), x[0], x[1], result['ais_delay'], peak, threshold)
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    hist.error_values.append(Err)
    sys.stdout.flush()
    return Err


def optimize_ais_delay(x):
    """
    Simplex does not perform well for finding a good parameter set for ais initiation of APs. This is an ad hoc method
    to find a ballpark to feed into simplex.
    :param x: array
    :return: array
    """
    perturb = (-0.1, 0.01)
    ais_delay_error(x)
    min_Err = 1e9
    max_iter = 80
    iter = 0
    best_x = list(x)
    while min_Err > 0. and iter < max_iter:
        index = random.randint(0, 1)
        candidate_x = list(x)
        candidate_x[index] += perturb[index]
        if not check_bounds.within_bounds(candidate_x, 'ais_delay'):
            print 'Process %i: Iteration %i - Parameters outside optimization bounds.' % (os.getpid(), iter)
        else:
            x = list(candidate_x)
            start_time = time.time()
            update_ais_delay(x)
            offset_vm('soma', v_active)
            sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
            duration = equilibrate + 100.
            sim.tstop = duration
            t = np.arange(0., duration, dt)
            spike = False
            d_amp = 0.01
            amp = i_th['soma'] - 0.02
            while not spike:
                sim.modify_stim(0, amp=amp)
                sim.run(v_active)
                vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
                if amp == np.any(vm[:int(equilibrate / dt)] > -30.):
                    print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
                    return 1e9
                if np.any(vm[int(equilibrate / dt):int((equilibrate + 50.) / dt)] > -30.):
                    spike = True
                else:
                    amp += d_amp
                    if sim.verbose:
                        print 'increasing amp to %.3f' % amp
            i_th['soma'] = amp
            peak, threshold, ADP, AHP = get_spike_shape(vm)
            result = {}
            result['soma_peak'] = peak
            result['th_v'] = threshold
            th_x = np.where(vm[int(equilibrate / dt):] >= threshold)[0][0] + int(equilibrate / dt)
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
            if axon_peak_t > ais_peak_t + dt:
                result['ais_delay'] = 0.
            else:
                result['ais_delay'] = ais_peak_t + dt - axon_peak_t
            Err = 0.
            target = 'ais_delay'
            Err += ((target_val['na_ka'][target] - result[target]) / target_range['na_ka'][target]) ** 2.
            # attempt to find the minimal combination that produces the desired delay
            print 'Simulation took %i s' % (time.time() - start_time)
            print 'Process %i (Iter %i): [ais.sha_nas, ais.gbar_nas]: [%.3E, %.3E], ais_delay: %.3E, ' \
                  'soma_peak: %.1f, threshold: %.1f' % (os.getpid(), iter, x[0], x[1], result['ais_delay'], peak,
                                                        threshold)
            print 'Process %i: Error: %.4E' % (os.getpid(), Err)
            if Err == 0.:
                return x
            if Err < min_Err:
                min_Err = Err
                best_x = list(x)
        iter += 1
        sys.stdout.flush()
    return best_x


equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.02
amp = 0.3
th_dvdt = 10.
v_init = -67.
v_active = -61.

cell = CA1_Pyr(morph_filename=morph_filename, mech_filename=mech_filename, full_spines=spines)
if spines is False:
    cell.correct_for_spines()
cell.set_terminal_branch_na_gradient()

axon_seg_locs = [seg.x for seg in cell.axon[2].sec]

rec_locs = {'soma': 0., 'ais': 1., 'axon': axon_seg_locs[0]}
rec_nodes = {'soma': cell.tree.root, 'ais': cell.axon[1], 'axon': cell.axon[2]}

sim = QuickSim(duration, cvode=False, dt=dt, verbose=False)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

for description, node in rec_nodes.iteritems():
    sim.append_rec(cell, node, loc=rec_locs[description], description=description)

sim.append_rec(cell, cell.axon[0], loc=0.5, description='axon_hill')
sim.append_rec(cell, cell.axon[2], loc=0.5, description='axon_center')
sim.append_rec(cell, cell.axon[2], loc=1.0, description='axon_end')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ca_i_CadepK',
               description='soma intra_Ca')


i_holding = {'soma': 0.00}
i_th = {'soma': 0.05}

#the target values and acceptable ranges
target_val = {}
target_range = {}
target_val['na_ka'] = {'v_rest': v_init, 'th_v': -51., 'soma_peak': 40., 'ADP': 0., 'AHP': 3.,
                       'stability': 0., 'ais_delay': 0., 'slow_depo': 20., 'dend_amp': 0.6}
target_range['na_ka'] = {'v_rest': 0.25, 'th_v': .05, 'soma_peak': 2., 'ADP': 0.01, 'AHP': .01,
                         'stability': 1., 'ais_delay': 0.001, 'slow_depo': 0.5, 'dend_amp': 0.005}

x0 = {}
xmin = {}
xmax = {}

axon_gbar_nax = cell.axon[2].sec.gbar_nax

xlabels = {}
xlabels['ais_delay'] = ['ais.sha_nas', 'ais.gbar_nax']

# x0['ais_delay'] = [-3.6, 0.4]
# xmin['ais_delay'] = [-5., 1.1*axon_gbar_nax]
# xmax['ais_delay'] = [-1., 5.*axon_gbar_nax]

x0['ais_delay'] = [-6., 0.41]
xmin['ais_delay'] = [-6., 1.1*axon_gbar_nax]
xmax['ais_delay'] = [-1., 6.*axon_gbar_nax]

hist = optimize_history()
hist.xlabels = xlabels['ais_delay']

check_bounds = CheckBounds(xmin, xmax)

x1 = x0['ais_delay']

max_niter = 1000  # max number of iterations to run
niter_success = 400  # max number of iterations without significant progress before aborting optimization

take_step = Normalized_Step(x0['ais_delay'], xmin['ais_delay'], xmax['ais_delay'])
minimizer_kwargs = dict(method=null_minimizer)

"""
result = optimize.basinhopping(ais_delay_error, x1, niter=max_niter,
                               niter_success=niter_success, disp=True, interval=20,
                               minimizer_kwargs=minimizer_kwargs, take_step=take_step)
x1 = result.x

result = optimize.minimize(ais_delay_error, x1, method='Nelder-Mead', options={'fatol': 1e-4, 'xatol': 1e-3,
                                                                               'disp': True, 'maxiter': 400})


# best_x = optimize_ais_delay([-1., 1.1*axon_gbar_nax])
# best_x = optimize_ais_delay(x0['ais_delay'])
# best_x = hist.report_best()
best_x = result.x
# best_x = x1
update_ais_delay(best_x)
cell.export_mech_dict(cell.mech_filename)
"""
ais_delay_error(x1, plot=True)