__author__ = 'Grace Ng'
from specify_cells2 import *
import os
import sys
from ipyparallel import interactive
import mkl

"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by a value for the amplitude of
a somatic current injection to test spike shape and stability.
"""

mkl.set_num_threads(1)

# morph_filename = 'EB2-late-bifurcation.swc'
# morph_filename = 'DG_GC_355549.swc'
neurotree_filename = '121516_DGC_trees.pkl'
neurotree_dict = read_from_pkl(morph_dir+neurotree_filename)

rec_filename = str(time.strftime('%m%d%Y', time.gmtime()))+'_'+str(time.strftime('%H%M%S', time.gmtime()))+\
               '_pid'+str(os.getpid())+'_sim_output'

x = None  # Placeholder for parameters pushed from controller.

""""
#[soma.Ca factor, soma.CadepK factor]
x['spike_adaptation'] = [0.1, 0.1]
"""

xmin = {}
xmax = {}

soma_na_gbar = 0.04
# [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor, soma.sh_nas/x, axon.gkdrbar factor,
#  dend.gkabar factor]
xmin['na_ka_stability'] = [0.01, 0.01, 1., 2., 0.1, 1., 1.]
xmax['na_ka_stability'] = [0.075, 0.05, 5., 5., 6., 2., 5.]

check_bounds = CheckBounds(xmin, xmax)

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    mech_filename = '030217 GC optimizing excitability'


@interactive
def get_spike_shape(vm):
    """

    :param t: array
    :param vm: array
    :return: tuple of float: (v_peak, th_v, ADP, AHP)
    """
    vm = vm[int((equilibrate+1.)/dt):]
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
def update_mech_dict():
    update_na_ka_stability(x)
    cell.export_mech_dict(cell.mech_filename)


@interactive
def update_na_ka_stability(x):
    """

    :param x: array [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor, soma.sh_nas/x,
                    axon.gkdrbar factor, dend.gkabar factor]
    """
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[1])
    cell.modify_mech_param('soma', 'kap', 'gkabar', x[0])
    slope = (x[6] - 1.) * x[0] / 300.
    cell.modify_mech_param('soma', 'nas', 'gbar', soma_na_gbar)
    cell.modify_mech_param('soma', 'nas', 'sh', x[4])
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
    cell.modify_mech_param('axon', 'nax', 'gbar', soma_na_gbar * x[3])
    cell.modify_mech_param('axon_hill', 'kap', 'gkabar', origin='soma')
    cell.modify_mech_param('axon_hill', 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param('axon_hill', 'nax', 'sh', x[4])
    for sec_type in ['ais', 'axon']:
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', x[1] * x[5])
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', x[0] * x[2])
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='axon_hill')


@interactive
def compute_spike_shape_features(local_x=None, plot=0):
    """
    :param local_x: array [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor, soma.sh_nas/x,
                    axon.gkdrbar factor, dend.gkabar factor]
    :return: float
    """
    if local_x is None:
        local_x = x
    if not check_bounds.within_bounds(local_x, 'na_ka_stability'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return None
    start_time = time.time()
    update_na_ka_stability(local_x)
    soma_vm = offset_vm('soma', v_active)
    result = {'v_rest': soma_vm}
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
    duration = equilibrate + 200.
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    d_amp = 0.01
    amp = i_th['soma'] - d_amp
    while not spike:
        sim.modify_stim(0, amp=amp)
        sim.run(v_active)
        vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
        if amp == 0.05 and np.any(vm[:int(equilibrate/dt)] > -30.):
            print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
            return None
        if np.any(vm[int(equilibrate/dt):int((equilibrate+50.)/dt)] > -30.):
            spike = True
        else:
            amp += d_amp
            if sim.verbose:
                print 'increasing amp to %.3f' % amp
    i_th['soma'] = amp
    peak, threshold, ADP, AHP = get_spike_shape(vm)
    dend_vm = np.interp(t, sim.tvec, sim.get_rec('dend')['vec'])
    th_x = np.where(vm[int(equilibrate / dt):] >= threshold)[0][0] + int(equilibrate / dt)
    dend_peak = np.max(dend_vm[th_x:th_x + int(10. / dt)])
    dend_pre = np.mean(dend_vm[th_x - int(0.2 / dt):th_x - int(0.1 / dt)])
    result['dend_amp'] = (dend_peak - dend_pre) / (peak - threshold)
    print 'Process %i took %.1f s to find spike rheobase at amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    if plot:
        sim.plot()
    result['v_th'] = threshold
    result['ADP'] = ADP
    result['AHP'] = AHP
    result['amp'] = amp
    return result

@interactive
def compute_spike_stability_features(amp, local_x=None, plot=0):
    sim.parameters['amp'] = amp
    if local_x is None:
        local_x = x
    if not check_bounds.within_bounds(local_x, 'na_ka_stability'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return 1e9
    start_time = time.time()
    update_na_ka_stability(local_x)
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=100.)
    duration = equilibrate + 200.
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    stability = 0.
    result = {}
    sim.modify_stim(0, amp=amp)
    sim.run(v_active)
    if plot:
        sim.plot()
    vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
    axon_vm = np.interp(t, sim.tvec, sim.get_rec('Axon Vm')['vec'])
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    v_before = np.max(vm[int((equilibrate - 50.)/dt):int((equilibrate - 1.)/dt)])
    v_after = np.max(vm[-int(50./dt):-1])
    stability += abs((v_before - v_rest) + (v_after - v_rest))
    v_min_late = np.min(vm[int((equilibrate + 80.)/dt):int((equilibrate + 99.)/dt)])
    result['stability'] = stability
    result['v_min_late'] = v_min_late
    print 'Process %i took %.1f s to test spike stability with amp: %.3f' % (os.getpid(), time.time()-start_time, amp)
    return result, axon_vm, t
"""
amp=0.6
axon_th=20.
start=250
stop=370
"""
def sim_spike_times(amp, axon_th, start, stop, local_x=None):
    """

    :param amp: float
    :param axon_th: float. This is the threshold above which the program will look for spike peaks
    :param start: float
    :param stop: float
    :param local_x: list
    :return:
    """
    if local_x is None:
        local_x = x
    result, axon_vm, t = compute_spike_stability_features(amp, local_x)
    axon_vm = axon_vm[int(start/dt):int(stop/dt)]
    t = t[int(start / dt):int(stop / dt)]
    #th_x is a list of the indices in axon_vm where the axon voltage is above the given threshold
    th_x = np.where(axon_vm > axon_th)[0]
    spike_times = []
    index = 0
    while index < len(th_x):
        peak_range = [th_x[index]]
        #look for indices of th_x that are consecutive; these indices are part of the same spike
        while (index+1) < len(th_x):
            if th_x[index+1] == (th_x[index]+1):
                peak_range.append(th_x[index+1])
                index +=1
            else:
                break
        if len(peak_range) > 1:
            v_peak = np.max(axon_vm[peak_range[0]:peak_range[-1]])
            peak_range_ind = np.where(axon_vm[peak_range[0]:peak_range[-1]] == v_peak)[0][0]
            x_peak = peak_range[peak_range_ind]
        else:
            x_peak = peak_range[0]
        spike_times.append(t[x_peak])
        index += 1
    return spike_times

def adapt_index(spike_times):
    """

    :param spike_times: list of the times at which there are spike peaks
    :return: adi is high for high spike adaptation (large differences in lengths of interspike intervals)
    """
    if len(spike_times) < 4:
        return None
    adi = 0
    count = 0
    isi = []
    for i in range(len(spike_times)-1):
        isi.append(spike_times[i+1] - spike_times[i])
    for i in range(len(isi)-1):
        adi += 1.0 * (isi[i+1] - isi[i]) / (isi[i+1] + isi[i])
        count += 1
    adi /= count
    return adi

def adjust_spike_number(target_spikes, axon_th, start, stop):
    spike_num = 0
    amp = 0.25
    while spike_num < target_spikes:
        spike_times = sim_spike_times(amp, axon_th, start, stop)
        spike_num = len(spike_times)
        amp += 0.05

def spike_adaptation_error(spike_num):
    amp = adjust_spike_number(spike_num)


@interactive
def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
        sim.export_to_file(f)


equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.01
amp = 0.3
th_dvdt = 10.
v_init = -77.
v_active = -77.

cell = DG_GC(neurotree_dict=neurotree_dict[0], mech_filename=mech_filename, full_spines=spines)
if spines is False:
    cell.correct_for_spines()

# get the thickest apical dendrite ~200 um from the soma
candidate_branches = []
candidate_diams = []
candidate_locs = []
for branch in cell.apical:
    if ((cell.get_distance_to_node(cell.tree.root, branch, 0.) >= 200.) &
            (cell.get_distance_to_node(cell.tree.root, branch, 1.) > 300.)):
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

rec_locs = {'soma': 0., 'dend': dend_loc}
rec_nodes = {'soma': cell.tree.root, 'dend': dend}

sim = QuickSim(duration, verbose=False)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

for description, node in rec_nodes.iteritems():
    sim.append_rec(cell, node, loc=rec_locs[description], description=description)

i_holding = {'soma': 0.}
i_th = {'soma': 0.05}

if type(cell.mech_dict['apical']['kap']['gkabar']) == list:
    orig_ka_dend_slope = \
        (element for element in cell.mech_dict['apical']['kap']['gkabar'] if 'slope' in element).next()['slope']
else:
    orig_ka_dend_slope = cell.mech_dict['apical']['kap']['gkabar']['slope']

orig_ka_soma_gkabar = cell.mech_dict['soma']['kap']['gkabar']['value']
orig_ka_dend_gkabar = orig_ka_soma_gkabar + orig_ka_dend_slope * 300.

x = [0.0308, 0.0220, 4.667, 4.808, 4.032, 1.297, 1.023]

sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ina_nas',
               description='Soma nas_i')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ik_kap',
               description='Soma kap_i')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ik_kdr',
               description='Soma kdr_i')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ik_km2',
               description='Soma km_i')
sim.append_rec(cell, cell.axon[2], loc=0.5, description='Axon Vm')

cell.modify_mech_param('soma', 'Ca', 'gcamult', .1)
cell.modify_mech_param('soma', 'CadepK', 'gcakmult', 0.1)
"""
for sec_type in ['axon_hill', 'ais', 'axon']:
    cell.modify_mech_param(sec_type, 'Ca', 'gcamult', origin='soma')
    cell.modify_mech_param(sec_type, 'CadepK', 'gcakmult', origin='soma')
"""

sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_i_Ca',
               description='Soma Ca_i')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ca_i_Ca',
               description='Soma Ca_intra_Ca')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_i_CadepK',
               description='Soma CadepK_i')

#Target adaptation index: 0.118989 for ~6 spikes