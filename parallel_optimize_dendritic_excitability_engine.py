__author__ = 'Grace Ng'
from specify_cells2 import *
import random
import os
import sys
from ipyparallel import interactive
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by string corresponding to the
section type to test R_inp.
"""

# morph_filename = 'EB2-late-bifurcation.swc'
morph_filename = 'DG_GC_355549.swc'
mech_filename = '102016 DG_GC pas no_spines'
rec_filename = str(time.strftime('%m%d%Y', time.gmtime()))+'_'+str(time.strftime('%H%M%S', time.gmtime()))+\
               '_pid'+str(os.getpid())+'_sim_output'

# placeholder for optimization parameter, must be pushed to each engine on each iteration
# x: array (soma.g_pas, dend.g_pas slope, dend.g_pas tau, dend.g_pas xhalf)
# x = [1.52E-06, 1.63E-05, 121.9, 150.]
x = [1.52E-06, 0.]

i_holding = {'soma': 0., 'dend': 0., 'distal_dend': 0.}


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
    return initial_v_rest


@interactive
def update_pas_linear(x):
    """
    x0 = [1.52E-06, 0.]
    :param x: array [soma.g_pas, dend.g_pas slope]
    """
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.modify_mech_param('apical', 'pas', 'g', origin='soma', slope=x[1])
    for sec_type in ['basal', 'axon_hill', 'axon', 'ais', 'trunk', 'apical', 'tuft', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')


@interactive
def update_pas_exp(x):
    """

    x0 = [2.28e-05, 1.58e-06, 58.4]
    :param x: array [soma.g_pas, dend.g_pas slope, dend.g_pas tau, dend.g_pas xmax]
    """
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    maxval = x[0] + x[1] * (np.exp(x[3]/x[2]) - 1.)
    cell.modify_mech_param('apical', 'pas', 'g', origin='soma', slope=x[1], tau=x[2], max=maxval)
    for sec_type in ['basal', 'axon_hill', 'axon', 'ais', 'trunk', 'apical', 'tuft', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')


@interactive
def update_pas_sigmoid(x):
    """

    x0 = [2.28e-05, 1.58e-06, 58.4, 200.]
    :param x: array (soma.g_pas, dend.g_pas slope, dend.g_pas tau, dend.g_pas xhalf)
    """
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.modify_mech_param('apical', 'pas', 'g', origin='soma', slope=x[1], tau=x[2], xhalf=x[3])
    for sec_type in ['basal', 'axon_hill', 'axon', 'ais', 'trunk', 'apical', 'tuft', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')


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
    cell.zero_na()
    if local_x is None:
        # update_pas_sigmoid(x)
        update_pas_exp(x)
    else:
        # update_pas_sigmoid(local_x)
        update_pas_exp(local_x)
    offset_vm(section)
    loc = rec_locs[section]
    node = rec_nodes[section]
    rec = sim.get_rec(section)
    sim.modify_stim(0, node=node, loc=loc, amp=amp, dur=stim_dur)
    sim.run(v_init)
    Rinp = get_Rinp(np.array(sim.tvec), np.array(rec['vec']), equilibrate, duration, amp)[2]
    result = {section: Rinp}
    print 'Process:', os.getpid(), 'calculated Rinp for %s in %i s, Rinp: %.1f' % (section, time.time() - start_time,
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
v_init = -67.
v_active = -61.

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False

if spines:
    sim_description = 'with_spines'
else:
    sim_description = 'no_spines'

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

# get the thickest terminal branch > 300 um from the soma
candidate_branches = []
candidate_diams = []
candidate_locs = []
for branch in (branch for branch in cell.apical if cell.is_terminal(branch)):
    if ((cell.get_distance_to_node(cell.tree.root, branch, 0.) < 300.) &
            (cell.get_distance_to_node(cell.tree.root, branch, 1.) > 300.)):
        candidate_branches.append(branch)
        for seg in branch.sec:
            loc = seg.x
            if cell.get_distance_to_node(cell.tree.root, branch, loc) > 300.:
                candidate_diams.append(branch.sec(loc).diam)
                candidate_locs.append(loc)
                break
index = candidate_diams.index(max(candidate_diams))
distal_dend = candidate_branches[index]
distal_dend_loc = candidate_locs[index]

rec_locs = {'soma': 0., 'dend': dend_loc, 'distal_dend': distal_dend_loc}
rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'distal_dend': distal_dend}

sim = QuickSim(duration, verbose=False)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

for description, node in rec_nodes.iteritems():
    sim.append_rec(cell, node, loc=rec_locs[description], description=description)


"""
axon_seg_locs = [seg.x for seg in cell.axon[2].sec]
sim.append_rec(cell, cell.axon[1], loc=1., description='ais')
sim.append_rec(cell, cell.axon[2], loc=axon_seg_locs[0], description='axon')
"""