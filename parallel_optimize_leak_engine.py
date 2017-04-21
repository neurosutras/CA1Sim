__author__ = 'Grace Ng'
from specify_cells2 import *
from function_lib import *
import time
import random
import os
import sys
from ipyparallel import interactive
import pprint
import mkl

mkl.set_num_threads(1)
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by string corresponding to the
section type to test R_inp.
"""

# morph_filename = 'EB2-late-bifurcation.swc'
# morph_filename = 'DG_GC_355549.swc'
neurotree_filename = '121516_DGC_trees.pkl'
neurotree_dict = read_from_pkl(morph_dir+neurotree_filename)

rec_filename = str(time.strftime('%m%d%Y', time.gmtime()))+'_'+str(time.strftime('%H%M%S', time.gmtime()))+\
               '_pid'+str(os.getpid())+'_sim_output'

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    if spines:
        mech_filename = '120116 DG_GC pas spines'
    else:
        mech_filename = '030217 GC optimizing excitability'

i_holding = {'soma': 0., 'dend': 0., 'distal_dend': 0.}

# placeholder for optimization parameter, must be pushed to each engine on each iteration
# x: array ['soma.g_pas', 'dend.g_pas slope', 'dend.g_pas tau', 'dend.gpas max_loc']
x = []


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
    update_pas_exp(x)
    cell.export_mech_dict(cell.mech_filename)

@interactive
def update_pas_exp(x):
    """

    x0 = [2.28e-05, 1.58e-06, 58.4]
    :param x: array [soma.g_pas, dend.g_pas slope, dend.g_pas tau]
    """
    if spines is False:
        cell.reinit_mechanisms(reset_cable=True)
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.modify_mech_param('apical', 'pas', 'g', origin='soma', slope=x[1], tau=x[2])
    for sec_type in ['axon_hill', 'axon', 'ais', 'apical', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')
    if spines is False:
        cell.correct_for_spines()


def print_gpas_cm_values():
    sec_types = ['apical']
    gpas_values = {s: [] for s in sec_types}
    cm_values = {s: [] for s in sec_types}
    for i in [0, 10, 20]:
        node = cell.get_nodes_of_subtype('apical')[i]
        for i, segment in enumerate(node.sec):
            node.sec.push()
            h.pop_section()
            gpas_values['apical'].append(node.sec(segment.x).g_pas)
            cm_values['apical'].append(node.sec(segment.x).cm)
    print 'g_pas: '
    pprint.pprint(gpas_values)
    print 'cm '
    pprint.pprint(cm_values)


@interactive
def get_Rinp_for_section(section, local_x=None):
    """
    Inject a hyperpolarizing step current into the specified section, and return the steady-state input resistance.
    :param section: str
    :return: dict: {str: float}
    """
    start_time = time.time()
    sim.tstop = duration
    sim.parameters['section'] = section
    sim.parameters['target'] = 'Rinp'
    sim.parameters['optimization'] = 'pas'
    sim.parameters['description'] = sim_description
    amp = -0.05
    if local_x is None:
        update_pas_exp(x)
    else:
        update_pas_exp(local_x)
    cell.zero_na()
    offset_vm(section)
    loc = rec_locs[section]
    node = rec_nodes[section]
    rec = sim.get_rec(section)
    sim.modify_stim(0, node=node, loc=loc, amp=amp, dur=stim_dur)
    sim.run(v_init)
    Rinp = get_Rinp(np.array(sim.tvec), np.array(rec['vec']), equilibrate, duration, amp)[2]
    result = {section: Rinp}
    print 'Process:', os.getpid(), 'calculated Rinp for %s in %.1f s, Rinp: %.1f' % (section, time.time() - start_time,
                                                                                    Rinp)
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

if spines:
    sim_description = 'with_spines'
else:
    sim_description = 'no_spines'

cell = DG_GC(neurotree_dict=neurotree_dict[0], mech_filename=mech_filename, full_spines=spines)
cell.zero_na()

# get the thickest apical dendrite ~200 um from the soma
candidate_branches = []
candidate_diams = []
candidate_locs = []
for branch in cell.apical:
    if ((cell.get_distance_to_node(cell.tree.root, branch, 0.) >= 250.) &
            (cell.get_distance_to_node(cell.tree.root, branch, 1.) > 300.)):
        candidate_branches.append(branch)
        for seg in branch.sec:
            loc = seg.x
            if cell.get_distance_to_node(cell.tree.root, branch, loc) > 300.:
                candidate_diams.append(branch.sec(loc).diam)
                candidate_locs.append(loc)
                break
index = candidate_diams.index(max(candidate_diams))
dend = candidate_branches[index]
dend_loc = candidate_locs[index]

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

rec_locs = {'soma': 0., 'dend': dend_loc, 'distal_dend': distal_dend_loc}
rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'distal_dend': distal_dend}

sim = QuickSim(duration, verbose=False)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

for description, node in rec_nodes.iteritems():
    sim.append_rec(cell, node, loc=rec_locs[description], description=description)
