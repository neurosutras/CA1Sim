__author__ = 'Grace Ng'
from specify_cells2 import *
import os
import sys
from ipyparallel import interactive
import mkl

"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by a value for the amplitude of
a somatic current injection to test spike adaptation.
"""

mkl.set_num_threads(1)

# morph_filename = 'EB2-late-bifurcation.swc'
# morph_filename = 'DG_GC_355549.swc'
neurotree_filename = '121516_DGC_trees.pkl'
neurotree_dict = read_from_pkl(morph_dir+neurotree_filename)

rec_filename = str(time.strftime('%m%d%Y', time.gmtime()))+'_'+str(time.strftime('%H%M%S', time.gmtime()))+\
               '_pid'+str(os.getpid())+'_sim_output'


# placeholder for optimization parameter, must be pushed to each engine on each iteration
# x: array [soma.gCa factor, soma.gCadepK factor, soma.gKm, axon(ais, hill).gKm factor]
x = None  # Placeholder for parameters pushed from controller.

xmin = {}
xmax = {}

#[soma.gCa factor, soma.gCadepK factor, soma.gKm, axon(ais, hill).gKm factor]
xmin['spike_adaptation'] = [0., 0., 0., 0.]
xmax['spike_adaptation'] = [5., 5., 1., 10.]
#Need to update these!

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
def update_mech_dict():
    update_spike_adaptation(x)
    cell.export_mech_dict(cell.mech_filename)

@interactive
def update_spike_adaptation(x):
    """

    :param x: array [soma.gCa factor, soma.gCadepK factor, soma.gKm, axon(ais, hill).gKm factor]
    """

    cell.modify_mech_param('soma', 'Ca', 'gcamult', x[0])
    cell.modify_mech_param('soma', 'CadepK', 'gcakmult', x[1])
    cell.modify_mech_param('soma', 'km2', 'gkmbar', x[2])
    cell.modify_mech_param('ais', 'km2', 'gkmbar', x[2]*x[3])
    for sec_type in ['axon', 'axon_hill']:
        cell.reinitialize_subset_mechanisms(sec_type, 'km2')

@interactive
def compute_spike_stability_features(amp, stim_dur, plot=0):
    sim.parameters['amp'] = amp
    start_time = time.time()
    sim.modify_stim(0, node=cell.tree.root, loc=0., dur=stim_dur)
    duration = equilibrate + stim_dur + 100.
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

@interactive
def sim_spike_times(amp, stim_dur=500., local_x=None):
    """

    :param amp: float
    :param axon_th: float. This is the threshold above which the program will look for spike peaks
    :param start: float
    :param stop: float
    :return:
    """
    if local_x is None:
        local_x = x #For some reason, local_x is still none!
    if not check_bounds.within_bounds(local_x, 'spike_adaptation'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return None
    update_spike_adaptation(local_x)
    spike_stab_result, axon_vm, t = compute_spike_stability_features(amp, stim_dur)
    axon_vm = axon_vm[int(spikes_start/dt):int(spikes_stop/dt)]
    t = t[int(spikes_start / dt):int(spikes_stop / dt)]
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
    result = {}
    result['spike_times'] = spike_times
    result['amp'] = amp
    return result

@interactive
def adjust_spike_number(target_spikes, local_x=None):
    if local_x is None:
        local_x = x
    print local_x
    if not check_bounds.within_bounds(local_x, 'spike_adaptation'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return None
    spike_num = 0
    d_amp = 0.01
    amp = i_th['soma'] - d_amp
    while spike_num < target_spikes:
        result = sim_spike_times(amp, 100., local_x)
        spike_times = result['spike_times']
        spike_num = len(spike_times)
        amp += 0.01
        if target_spikes == 1 and amp > 2.:
            return None
    return result


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

#The time interval (in ms) during the recording during which the program looks for axonal spikes
spikes_start = 250.
spikes_stop = 350.
#Threshold above which an axonal spike will be counted
axon_th = 10.

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

