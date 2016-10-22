__author__ = 'Grace Ng'
from specify_cells2 import *
import random
import os
from ipyparallel import interactive
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by string corresponding to the
section type to test R_inp.
"""

# morph_filename = 'EB2-late-bifurcation.swc'
morph_filename = 'DG_GC_355549.swc'
mech_filename = '102016 DG_GC pas no_spines'


#x: array (soma.g_pas, trunk.g_pas slope, trunk.g_pas tau)
x = [1.52E-06, 1.63E-05, 121.9, 150.] # placeholder for optimization parameter, must be pushed to each engine on each iteration
xmin = [1.0E-7, 1.0E-07, 25.]
xmax = [2.0E-5, 2.0E-05, 400.]
i_holding = {'soma': 0., 'dend': 0., 'distal_dend': 0.}


def check_bounds(x, param_name):
    """
    For optimize_polish, based on simplex algorithm, check that the current set of parameters are within the bounds.
    :param x: array
    :param param_name: str
    :return: bool
    """
    for i in range(len(x)):
        if ((xmin[i] is not None and x[i] < xmin[i]) or
                (xmax[i] is not None and x[i] > xmax[i])):
            return False
    return True


@interactive
def offset_vm(description, vm_target=None):
    """

    :param sec_type: str
    """
    if vm_target is None:
        vm_target = v_init
    sim.modify_stim(0, amp=0.)
    node = rec_nodes[description]
    loc = rec_locs[description]
    rec_dict = sim.get_rec(description)
    sim.modify_stim(1, node=node, loc=loc, amp=0.)
    rec = rec_dict['vec']
    # stopped editing here 102016
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


@interactive
def update_pas(x):
    """

    x0 = [2.28e-05, 1.58e-06, 58.4, ]
    :param x: array (soma.g_pas, dend.g_pas slope, dend.g_pas tau, dend.g_pas xhalf)
    """
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.modify_mech_param('apical', 'pas', 'g', origin='soma', slope=x[1], tau=x[2], xhalf=x[3])
    for sec_type in ['basal', 'axon_hill', 'axon', 'ais', 'trunk', 'apical', 'tuft', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')


@interactive
def section_pas_error(section):
    """
    :param section: str
    :return: float
    """
    start_time = time.time()
    sim.tstop = duration
    amp = -0.15
    cell.zero_na()
    update_pas(x)
    offset_vm(section)
    if section == 'soma':
        location = 0.
        sec_type = cell.tree.root
        dict_index = 0
    elif section == 'trunk':
        location = 1.
        sec_type = trunk
        dict_index = 1
    elif section == 'tuft':
        location = 1.
        sec_type = tuft
        dict_index = 2
    sim.modify_stim(0, node=sec_type, loc=location, amp=amp, dur=stim_dur)
    sim.run(v_init)
    result = dict([(section, get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[dict_index]['vec']), equilibrate, duration, amp)[2])])
    print 'Process:', os.getpid(), 'finished simulation for %s, Time: %.3f s, Rinp: %.2f' % (section, time.time() - start_time, result.get(section))
    return result


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

# get the thickest apical dendrite ~200 um from the soma
candidate_branches = []
candidate_diams = []
candidate_locs = []
for branch in cell.apical:
    if ((cell.get_distance_to_node(cell.tree.root, branch, 0.) < 200.) &
            (cell.get_distance_to_node(cell.tree.root, branch, 1.) > 200.)):
        candidate_branches.append(branch)
        for seg in branch.sec:
            loc = seg.x
            if cell.get_distance_to_node(cell.tree.root, branch, loc) > 200.:
                candidate_diams.append(branch.sec(loc).diam)
                candidate_locs.append(loc)
                break
index = candidate_diams.index(max(candidate_diams))
dend = candidate_branches[index]
dend_loc = candidate_locs[index]

# get a terminal branch > 300 um from the soma
candidate_branches = []
candidate_distances = []
for branch in cell.apical:
    distance = cell.get_distance_to_node(cell.tree.root, branch, 0.)
    if (distance > 300.) & cell.is_terminal(branch):
        candidate_branches.append(branch)
        candidate_distances.append(distance)
index = candidate_distances.index(max(candidate_distances))
distal_dend = candidate_branches[index]
distal_dend_loc = 0.

rec_locs = {'soma': 0., 'dend': dend_loc, 'distal_dend': distal_dend_loc}
rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'distal_dend': distal_dend}

pas_node_keys = {'soma': 'soma', 'dend': 'dend'}

sim = QuickSim(duration)  # , verbose=False)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

for description, node in rec_nodes.iteritems():
    sim.append_rec(cell, node, loc=rec_locs[description], description=description)

"""
axon_seg_locs = [seg.x for seg in cell.axon[2].sec]
sim.append_rec(cell, cell.axon[1], loc=1., description='ais')
sim.append_rec(cell, cell.axon[2], loc=axon_seg_locs[0], description='axon')
"""