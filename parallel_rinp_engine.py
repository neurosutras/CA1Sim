__author__ = 'milsteina'
from specify_cells2 import *
import os
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which section to simulate (full sampling of sections).
"""

#morph_filename = 'EB2-late-bifurcation.swc'
#mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'

morph_filename = 'DG_GC_355549.swc'
mech_filename = '110316 DG_GC pas spines'
#mech_filename = '110316 DG_GC pas no_spines'

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())


def test_single_section(sec_index):
    """
    :param sec_index: int
    :return: str
    """
    start_time = time.time()
    node = nodes[sec_index]
    loc = 0. if node.type == 'soma' else 0.5
    sim.modify_rec(0, node=node, loc=loc)
    sim.modify_stim(0, node=node, loc=loc)
    sim.run(v_init)
    Rinp_params = get_Rinp(sim.tvec, sim.rec_list[0]['vec'], equilibrate + delay, equilibrate + delay + stim_dur, amp, dt=0.02)
    sim.parameters['Rinp_baseline'] = Rinp_params[0]
    sim.parameters['Rinp_steady'] = Rinp_params[1]
    sim.parameters['Rinp_peak'] = Rinp_params[2]
    sim.parameters['decay_50'] = calc_decay(0.5, sim.tvec, sim.rec_list[0]['vec'], equilibrate + delay, equilibrate + delay + stim_dur)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, sec_index)
    print 'Process:', os.getpid(), 'completed Iteration:', sec_index, 'Node:', node.name, 'in', \
        time.time() - start_time, 's'
    return rec_filename

def calc_decay(decay_degree, tvec, vec, start, stop, dt = 0.02):
    """
    Calculates the time point by which the function has decayed by a certain degree

    :param decay_degree: degree of decay (ex. 0.50, 0.70)
    :param tvec: vector of time points
    :param vec: vector of recorded voltages
    :param start: time of start of voltage trace
    :param stop: time of end of voltage trace
    :param dt:
    :return:
    """
    interp_arrs = interpolate_tvec_vec(tvec, vec, stop - start, dt)
    interp_t = interp_arrs[0]
    interp_vm = interp_arrs[1]
    left = int((start - 3.) / dt)
    right = left + int(2. / dt)
    baseline = np.mean(interp_vm[left:right])

    left = int(start / dt)
    right = int(stop / dt)
    interp_t = interp_t[left:right]
    interp_vm = np.abs(interp_vm[left:right] - baseline)
    amp = np.max(interp_vm)
    t_peak = np.where(interp_vm == amp)[0][0]
    interp_vm /= amp
    interp_t -= interp_t[0]
    decay_point = np.where(interp_vm[0:t_peak] >= decay_degree)[0][0]
    return interp_t[decay_point];


equilibrate = 250.  # time to steady-state
delay = 10.
stim_dur = 750.
duration = equilibrate + delay + stim_dur + 100.
amp = -0.15
v_init = -67.

#cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell = DG_GC(morph_filename, mech_filename, full_spines=True)
cell.zero_na()

nodes = cell.soma+cell.basal+cell.trunk+cell.apical+cell.tuft+cell.axon

sim = QuickSim(duration, verbose=False)
sim.append_rec(cell, cell.tree.root, 0.)
sim.append_stim(cell, cell.tree.root, 0., amp, equilibrate+delay, stim_dur)