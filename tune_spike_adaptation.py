__author__ = 'Grace Ng'
from specify_cells2 import *
import os
import sys

"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by a value for the amplitude of
a somatic current injection to test spike shape and stability.
"""

# morph_filename = 'EB2-late-bifurcation.swc'
# morph_filename = 'DG_GC_355549.swc'
neurotree_filename = '121516_DGC_trees.pkl'
neurotree_dict = read_from_pkl(morph_dir+neurotree_filename)

rec_filename = str(time.strftime('%m%d%Y', time.gmtime()))+'_'+str(time.strftime('%H%M%S', time.gmtime()))+\
               '_pid'+str(os.getpid())+'_sim_output'

#[soma.gCa factor, soma.gCadepK factor, soma.gKm, axon(ais, hill).gKm factor]
xmin = []
xmax = []
#Update xmin and xmax!


#Check bounds to make sure soma gKm and axon gKm are within a reasonable range
#check_bounds = CheckBounds(xmin, xmax)


if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    mech_filename = '030217 GC optimizing excitability'



def update_mech_dict():
    update_spike_adaptation(x)
    cell.export_mech_dict(cell.mech_filename)


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

def sim_spike_times(amp, stim_dur, axon_th=10, start=250., stop=360.):
    """

    :param amp: float
    :param axon_th: float. This is the threshold above which the program will look for spike peaks
    :param start: float
    :param stop: float
    :return:
    """
    result, axon_vm, t = compute_spike_stability_features(amp, stim_dur)
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
    :return: adi is a large value for high spike adaptation (large differences in lengths of interspike intervals)
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

def adjust_spike_number(target_spikes, stim_dur, axon_th=10, start=250., stop=350.):
    spike_num = 0
    spikes_amp_corr = []
    d_amp = 0.01
    amp = i_th['soma'] - d_amp
    """
    while spike_num < 1:
        spike_times = sim_spike_times(amp, stim_dur, axon_th, start, stop)
        spike_num = len(spike_times)
        amp += 0.01
    spikes_amp_corr.append([amp, spike_num])
    """
    while spike_num < target_spikes:
        spike_times = sim_spike_times(amp, stim_dur, axon_th, start, stop)
        spike_num = len(spike_times)
        amp += 0.01
    """
    spikes_amp_corr.append([amp, spike_num])
    for amp_incr in [0.05, 0.1]:
        test_amp = amp + amp_incr
        test_spike_times = sim_spike_times(amp, stim_dur, axon_th, start, stop)
        spikes_amp_corr.append([test_amp, len(test_spike_times)])
    """
    return spike_times, amp

def spike_adaptation_error(spike_num, axon_th, start, stop):
    #Need to consider rheobase: the current to cross threshold for a single spike.
    spike_times, rheobase = adjust_spike_number(1, 100., axon_th, start, stop)
    spike_times, amp = adjust_spike_number(spike_num, 100., axon_th, start, stop)
    adi = adapt_index(spike_times)

    freq = [len(spike_times)/(stop-start)]
    curr_inj = [amp]
    #Calculate more frequencies using stim duration of 500 ms
    for amp_incr in [0.3, 0.6]:
        curr_spike_times, curr_amp = sim_spike_times(amp+amp_incr, 500., axon_th, start, stop+400.)
        freq.append(len(curr_spike_times)/(stop+400.-start))
        curr_inj.append(curr_amp)


    #Consider error from adaptation index and the error from the correlation of the spike frequency (spike_count/0.1s stim duration)
    # v. current injection amp, and the error from the desired slope of the spike frequency v. current injection
    #Calculate overall error: sum of all (target-test)^2/desired_variance

    # Target adaptation index: 0.118989 for ~6 spikes



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

# These values will now be saved in the mech dictionary that is updated by previous round of optimization
#[soma.gCa factor, soma.gCadepK factor, soma.gKm, axon(ais, hill).gKm factor]
#Need to adjust the gKm values!!
x = [0.1, 0.1, 0.1, 1]


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

