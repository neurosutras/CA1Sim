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
morph_filename = 'DGC_356550_edited.swc'
mech_filename = '102016 DG_GC pas no_spines'


#x: array (soma.g_pas, trunk.g_pas slope, trunk.g_pas tau)
x = [1.52E-06, 1.63E-05, 121.9] # placeholder for optimization parameter, must be pushed to each engine on each iteration
xmin = [1.0E-6, 1.0E-06, 25.]
xmax = [2.0E-5, 2.0E-05, 200.]
i_holding = {'soma': 0., 'trunk': 0., 'tuft': 0.}


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

@interactive
def update_pas(x):
    """

    x0 = [2.28e-05, 1.58e-06, 58.4]
    :param x: array (soma.g_pas, trunk.g_pas slope, trunk.g_pas tau)
    """
    current_tuft_e_pas = -77.
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

@interactive
def section_pas_error(section):
    """
    :param section: str
    :return: float
    """
    start_time = time.time()
    sim.tstop = duration
    amp = -0.15
    cell.zero_h()
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

    """
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
        Err += ((target_val['pas'][target] - result[target]) / target_range['pas'][target]) ** 2.
    return Err
"""


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

# get the thickest apical dendrite 200 um from the soma
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


sim = QuickSim(duration)  # , verbose=False)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)
sim.append_rec(cell, cell.tree.root, loc=0.)
sim.append_rec(cell, dend, loc=dend_loc)

"""
axon_seg_locs = [seg.x for seg in cell.axon[2].sec]
sim.append_rec(cell, cell.axon[1], loc=1., description='ais')
sim.append_rec(cell, cell.axon[2], loc=axon_seg_locs[0], description='axon')
"""