__author__ = 'Grace Ng'
from specify_cells2 import *
from plot_results import *
import sys
import os
import random

# morph_filename = 'EB2-late-bifurcation.swc'
# morph_filename = 'DG_GC_355549.swc'
neurotree_filename = '121516_DGC_trees.pkl'
neurotree_dict = read_from_pkl(morph_dir+neurotree_filename)

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    # Start with a mechanism dictionary that has soma_ek = -77., soma_na_gbar = 0.04 and has set the other parameters accordingly
    mech_filename = '012416 GC optimizing excitability'

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

def update_ais_delay(x):
    """

    :param x: array [ais.sha_nas, ais.gbar_nax]
    """
    cell.modify_mech_param('ais', 'nax', 'sha', x[0])
    cell.modify_mech_param('ais', 'nax', 'gbar', x[1])

def na_ka_dend_error(x, plot=0):
    """

    :param x: array [soma.sh_nas, trunk.ka factor]
    :param plot: int
    :return: float
    """
    if not check_bounds.within_bounds(x, 'na_ka_dend'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9
    start_time = time.time()
    update_na_ka_dend(x)
    history.x_values.append(x)
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
        vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
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
    distal_dend_vm = np.interp(t, sim.tvec, sim.get_rec('distal_dend')['vec'])
    th_x = np.where(vm[int(equilibrate/dt):] >= threshold)[0][0] + int(equilibrate/dt)
    distal_dend_peak = np.max(distal_dend_vm[th_x:th_x+int(10./dt)])
    distal_dend_pre = np.mean(distal_dend_vm[th_x-int(0.2/dt):th_x-int(0.1/dt)])
    result['distal_dend_amp'] = (distal_dend_peak - distal_dend_pre) / (peak - threshold)
    if plot:
        sim.plot()
    Err = 0.
    for target in result:
        if target not in history.features:
            history.features[target] = []
        history.features[target].append(result[target])
        Err += ((target_val['na_ka'][target] - result[target])/target_range['na_ka'][target])**2.
    history.error_values.append(Err)
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [sh_nax/s, trunk.ka factor]: [%.2f, %.2f], threshold: %.1f, trunk_amp: %.1f' % \
          (os.getpid(), x[0], x[1], threshold, result['distal_dend_amp'])
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
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
    if not check_bounds.within_bounds(x, 'ais_delay'):
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
        vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
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
    ais_vm = np.interp(t, sim.tvec, sim.get_rec('ais')['vec'])
    ais_dvdt = np.gradient(ais_vm, [dt])
    axon_vm = np.interp(t, sim.tvec, sim.get_rec('axon')['vec'])
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
    for i, x_i in enumerate(x):
        err_factor = 1.
        Err += err_factor*(((x_i - xmin['ais_delay'][i])/min(abs(xmin['ais_delay'][i]), abs(xmax['ais_delay'][i])))**2.)
    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [ais.sha_nas, ais.gbar_nas]: [%.2f, %.2f], ais_delay: %.3E, soma_peak: %.1f, ' \
          'threshold: %.1f' % (os.getpid(), x[0], x[1], result['ais_delay'], peak, threshold)
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    if Err < min(history.error_values):
        history.error_values.append(float(Err))
        history.x_values.append(x)
        history.features['ais_delay'].append(result['ais_delay'])
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
    min_Err = 1e9
    max_iter = 20
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
            duration = equilibrate + 200.
            sim.tstop = duration
            t = np.arange(0., duration, dt)
            spike = False
            amp = 0.05
            while not spike:
                sim.modify_stim(0, amp=amp)
                sim.run(v_active)
                vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
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
            ais_vm = np.interp(t, sim.tvec, sim.get_rec('ais')['vec'])
            ais_dvdt = np.gradient(ais_vm, [dt])
            axon_vm = np.interp(t, sim.tvec, sim.get_rec('axon')['vec'])
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
            print 'Process %i (Iter %i): [ais.sha_nas, ais.gbar_nas]: [%.2f, %.2f], ais_delay: %.3E, ' \
                  'soma_peak: %.1f, threshold: %.1f' % (os.getpid(), iter, x[0], x[1], result['ais_delay'], peak,
                                                        threshold)
            print 'Process %i: Error: %.4E' % (os.getpid(), Err)
            if Err == 0.:
                return x
            if Err < min_Err:
                min_Err = Err
                best_x = list(x)
        iter += 1
    return best_x


equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.01
amp = 0.3
th_dvdt = 10.
v_init = -67.
# v_active = -61.
v_active = -67.

cell = DG_GC(neurotree_dict=neurotree_dict[0], mech_filename=mech_filename, full_spines=spines)

# get the thickest terminal branch > 300 um from the soma
candidate_branches = []
candidate_end_distances = []
for branch in (branch for branch in cell.apical if cell.is_terminal(branch)):
    if cell.get_distance_to_node(cell.tree.root, branch, 0.) >= 300.:
        candidate_branches.append(branch)
        candidate_end_distances.append(cell.get_distance_to_node(cell.tree.root, branch, 1.))
index = candidate_end_distances.index(max(candidate_end_distances))
distal_dend = candidate_branches[index]
distal_dend_loc = 1.

axon_seg_locs = [seg.x for seg in cell.axon[2].sec]

rec_locs = {'soma': 0., 'ais': 1., 'axon': axon_seg_locs[0], 'distal_dend': distal_dend_loc}
rec_nodes = {'soma': cell.tree.root, 'ais': cell.axon[1], 'axon': cell.axon[2], 'distal_dend': distal_dend}

sim = QuickSim(duration)  # , verbose=False)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

for description, node in rec_nodes.iteritems():
    sim.append_rec(cell, node, loc=rec_locs[description], description=description)

i_holding = {'soma': 0.09, 'primary': 0.09, 'secondary': 0.09}
compartment_objects = {'soma': cell.tree.root}

#the target values and acceptable ranges
target_val = {}
target_range = {}
target_val['na_ka'] = {'v_rest': v_init, 'th_v': -49., 'soma_peak': 40., 'distal_dend_amp': 0.5, 'ADP': 0., 'AHP': 4.,
                       'stability': 0., 'ais_delay': 0., 'slow_depo': 25.}
target_range['na_ka'] = {'v_rest': 0.25, 'th_v': .2, 'soma_peak': 2., 'distal_dend_amp': 0.01, 'ADP': 0.01, 'AHP': .2,
                         'stability': 1., 'ais_delay': 0.001, 'slow_depo': 1.}

x0 = {}
xmin = {}
xmax = {}

axon_gbar_nax = cell.axon[2].sec.gbar_nax

xlabels = {}
xlabels['ais_delay'] = ['ais.sha_nas', 'ais.gbar_nax']
xlabels['na_ka_dend'] = ['soma.sh_nas', 'trunk.ka factor']

x0['ais_delay'] = [-1.20, 1.57*axon_gbar_nax]  # Error: 29.16
xmin['ais_delay'] = [-5., 1.1*axon_gbar_nax]
xmax['ais_delay'] = [-1., 5.*axon_gbar_nax]

x0['na_ka_dend'] = [1.97, 1.21]  # Error: 13.57
xmin['na_ka_dend'] = [0.1, 1.1]
xmax['na_ka_dend'] = [4., 5.]

check_bounds = CheckBounds(xmin, xmax)
print 'axon_gbar_nax %.3f' % axon_gbar_nax

explore_niter = 600  # max number of iterations to run
polish_niter = 400
take_step = Normalized_Step(x0['na_ka_dend'], xmin['na_ka_dend'], xmax['na_ka_dend'])
minimizer_kwargs = dict(method=null_minimizer)

history = optimize_history()
history.xlabels = xlabels['na_ka_dend']

result = optimize.basinhopping(na_ka_dend_error, x0['na_ka_dend'], niter=explore_niter, niter_success=explore_niter,
                               disp=True, interval=20, minimizer_kwargs=minimizer_kwargs, take_step=take_step)
polished_result = optimize.minimize(na_ka_dend_error, result.x, method='Nelder-Mead', options={'ftol': 1e-5,
                                                    'disp': True, 'maxiter': polish_niter})
print polished_result
best_x = history.report_best()
update_na_ka_dend(best_x)

best_x = optimize_ais_delay([-1., 1.1*axon_gbar_nax])
update_ais_delay(best_x)
cell.export_mech_dict(cell.mech_filename)

