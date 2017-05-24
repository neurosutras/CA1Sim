__author__ = 'milsteina' and 'Grace Ng'
# from specify_cells2 import *
from specify_cells3 import *
import random
import os
import sys
from ipyparallel import interactive
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which synapse to stimulate (coarse sampling of the full set of spines).
"""

neurotree_filename = '121516_DGC_trees.pkl'
neurotree_dict = read_from_pkl(morph_dir+neurotree_filename)

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    # mech_filename = '042617 GC optimizing spike stability'
    mech_filename = '051917 GC optimizing EPSP'


def offset_vm(node, loc, index):
    """


    """
    dt = 0.02
    sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    offset = True
    global i_holding
    direction = None
    while offset:
        sim.modify_stim(0, node=node, loc=loc, amp=i_holding)
        sim.run(v_init)
        rec = sim.rec_list[index]['vec']
        vm = np.interp(t, sim.tvec, rec)
        v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
        if v_rest < v_init - 1.:
            i_holding += 0.005
            if sim.verbose:
                print 'increasing i_holding to %.3f' % (i_holding)
            if direction is None:
                direction = 1
            elif direction == -1:
                break
        elif v_rest > v_init + 1.:
            i_holding -= 0.005
            if sim.verbose:
                print 'decreasing i_holding to %.3f' % (i_holding)
            if direction is None:
                direction = -1
            elif direction == 1:
                break
        else:
            offset = False
    sim.tstop = duration

@interactive
def stimulate_single_synapse(syn_index):
    """
    :param syn_index: int
    :return: str
    """
    start_time = time.time()
    syn = syn_list[syn_index]
    if spines:
        spine = syn.node
        branch = spine.parent.parent
        sim.modify_rec(1, branch)
        sim.modify_rec(2, spine)
    else:
        branch = syn.node
        sim.modify_rec(1, branch)
        sim.modify_rec(2, branch, loc=syn.loc)
    sim.parameters['input_loc'] = branch.type
    syn.source.play(spike_times)
    offset_vm(branch, syn.loc, 2)
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, syn_index)
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    print 'Process:', os.getpid(), 'completed synapse index:', syn_index, ', Branch:', \
        branch.name, ', Loc: %.2f in %i s' % (syn.loc, time.time() - start_time)
    return rec_filename


equilibrate = 250.  # time to steady-state
duration = 450.
v_init = -77.
NMDA_type = 'NMDA_KIN5'
#syn_types = ['AMPA_KIN', NMDA_type]
syn_types = ['AMPA_KIN']
local_random = random.Random()
i_holding = 0.

syn_list = []
cell = DG_GC(neurotree_dict=neurotree_dict[0], mech_filename=mech_filename, full_spines=spines)
if spines is False:
    cell.correct_for_spines()

cell.zero_na()
#cell.zero_h()

local_random.seed(0)
for branch in cell.apical:
    if spines:
        if len(branch.spines) > 1:
            if branch.sec.L <= 10.:
                node = branch.spines[local_random.sample(range(0, len(branch.spines)), 1)[0]]
                syn = Synapse(cell, node, syn_types, stochastic=0)
                syn_list.append(syn)
            else:
                num_syns = min(len(branch.spines), int(branch.sec.L // 10.))  # a random synapse every 10 um
                for i in local_random.sample(range(0, len(branch.spines)), num_syns):
                    node = branch.spines[i]
                    syn = Synapse(cell, node, syn_types, stochastic=0)
                    syn_list.append(syn)
        elif branch.spines:
            node = branch.spines[0]
            # syn = Synapse(cell, node, [syn_type], stochastic=0)
            syn = Synapse(cell, node, syn_types, stochastic=0)
            syn_list.append(syn)
    else:
        syn_loc_list = branch.synapse_locs['excitatory']
        if len(syn_loc_list) > 1:
            if branch.sec.L <= 10.:
                syn_loc = local_random.sample(syn_loc_list, 1)
                #syn = Synapse(cell, node, [syn_type], stochastic=0)
                syn = Synapse(cell, branch, syn_types, loc=syn_loc, stochastic=0)
                syn_list.append(syn)
            else:
                num_syns = min(len(syn_loc_list), int(branch.sec.L//10.))  # a random synapse every 10 um
                for syn_loc in local_random.sample(syn_loc_list, num_syns):
                    syn = Synapse(cell, branch, syn_types, loc=syn_loc, stochastic=0)
                    syn_list.append(syn)
        elif syn_loc_list:
            syn_loc = syn_loc_list[0]
            syn = Synapse(cell, branch, syn_types, loc=syn_loc, stochastic=0)
            syn_list.append(syn)
# cell.modify_mech_param('apical', 'synapse', 'gmax', value=0.002, syn_type='AMPA_KIN')
cell.init_synaptic_mechanisms()

sim = QuickSim(duration)  # , verbose=False)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration


# get the thickest apical dendrite ~200 um from the soma
candidate_branches = []
candidate_diams = []
candidate_locs = []
for branch in cell.apical:
    if ((cell.get_distance_to_node(cell.tree.root, branch, 0.) >= 200.) &
            (cell.get_distance_to_node(cell.tree.root, branch, 1.) > 300.) & (not cell.is_terminal(branch))):
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
axon_seg_locs = [seg.x for seg in cell.axon[2].sec]

rec_locs = {'soma': 0., 'dend': 0.}
rec_nodes = {'soma': cell.tree.root, 'dend': dend}

sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

for description, node in rec_nodes.iteritems():
    sim.append_rec(cell, node, loc=rec_locs[description], description=description)
sim.append_rec(cell, dend, loc=0., description='synapse')

spike_times = h.Vector([equilibrate])