__author__ = 'milsteina'
from specify_cells2 import *
import os
import sys
from ipyparallel import interactive
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which section to simulate (full sampling of sections).
"""

# morph_filename = 'DG_GC_355549.swc'
neuroH5_filename = '121516_DGC_trees.pkl'
neuroH5_dict = read_from_pkl(morph_dir+neuroH5_filename)
rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())

spines = False

if len(sys.argv) > 1:
    mech_filename = str(sys.argv[1])
else:
    mech_filename = None


@interactive
def correct_for_spines_node(node, cm_correction_fraction):
    """
    If not explicitly modeling spine compartments for excitatory synapses, this method scales cm and g_pas in this
    dendritic section proportional to the number of excitatory synapses contained in the section.
    :param node: :class:'SHocNode'
    :param cm_correction_fraction: float in [0., 1.]
    """
    SA_spine = math.pi * (1.58 * 0.077 + 0.5 * 0.5)
    if 'excitatory' in node.synapse_locs:
        these_synapse_locs = np.array(node.synapse_locs['excitatory'])
        seg_width = 1. / node.sec.nseg
        for i, segment in enumerate(node.sec):
            node.sec.push()
            SA_seg = h.area(segment.x)
            h.pop_section()
            num_spines = len(np.where((these_synapse_locs >= i * seg_width) &
                                      (these_synapse_locs < (i + 1) * seg_width))[0])
            cm_correction_factor = (SA_seg + cm_correction_fraction * num_spines * SA_spine) / SA_seg
            soma_g_pas = node.sec.cell().mech_dict['soma']['pas']['g']['value']
            gpas_correction_factor = (SA_seg * node.sec(segment.x).g_pas + num_spines * SA_spine * soma_g_pas) \
                                     / (SA_seg * node.sec(segment.x).g_pas)
            node.sec(segment.x).g_pas *= gpas_correction_factor
            node.sec(segment.x).cm *= cm_correction_factor


@interactive
def correct_for_spines_cell(cm_correction_fraction):
    """
    If not explicitly modeling spine compartments for excitatory synapses, this method scales cm and g_pas in all
    dendritic sections proportional to the number of excitatory synapses contained in each section.
    :param cell: :class:'SHocCell'
    :param cm_correction_fraction: float in [0., 1.]
    """
    cell.reinit_mechanisms(reset_cable=True)
    for loop in range(2):
        for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
            for node in cell.get_nodes_of_subtype(sec_type):
                correct_for_spines_node(node, cm_correction_fraction)
                if loop == 0:
                    node.init_nseg()
                    node.reinit_diam()
        if loop == 0:
            cell.reinit_mechanisms()


@interactive
def test_single_section(sec_index, loc=None):
    """
    :param sec_index: int
    :return: str
    """
    start_time = time.time()
    node = nodes[sec_index]
    if loc is None:
        loc = 0. if node.type == 'soma' else 0.5
    sim.modify_rec(0, node=node, loc=loc)
    sim.modify_stim(0, node=node, loc=loc)
    sim.run(v_init)
    Rinp_params = get_Rinp(sim.tvec, sim.rec_list[0]['vec'], equilibrate + delay, equilibrate + delay + stim_dur, amp,
                           dt=0.02)
    sim.parameters['Rinp_baseline'] = Rinp_params[0]
    sim.parameters['Rinp_steady'] = Rinp_params[1]
    sim.parameters['Rinp_peak'] = Rinp_params[2]
    sim.parameters['decay_90'] = calc_decay(0.9, sim.tvec, sim.rec_list[0]['vec'], equilibrate + delay,
                                            equilibrate + delay + stim_dur)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, sec_index)
    print 'Process:', os.getpid(), 'completed Iteration:', sec_index, 'Node:', node.name, 'in', \
        time.time() - start_time, 's'
    if cell.is_terminal(node):
        test_single_section(sec_index, 1.)
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
    return interp_t[decay_point]


equilibrate = 250.  # time to steady-state
delay = 10.
stim_dur = 750.
duration = equilibrate + delay + stim_dur + 100.
amp = -0.15
v_init = -67.

#cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell = DG_GC(neuroH5_dict=neuroH5_dict[0], mech_filename=mech_filename, full_spines=spines)
if not spines:
    cell.correct_for_spines()
cell.zero_na()

nodes = cell.basal+cell.trunk+cell.apical+cell.tuft

sim = QuickSim(duration, verbose=False)
sim.append_rec(cell, cell.tree.root, 0.)
sim.append_stim(cell, cell.tree.root, 0., amp, equilibrate+delay, stim_dur)