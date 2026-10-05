__author__ = 'Aaron D. Milstein'
from specify_cells import *
import os
import random
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which synapse to stimulate (coarse sampling of the full set of spines).
"""

morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())


def calculate_single_attenuation_ratio(syn_index):
    """
    :param syn_index: int
    :return: str
    """
    start_time = time.time()
    sim.parameters['stim_loc'] = 'spine'
    syn = spine_syn_list[syn_index]
    loc = syn.loc
    spine = syn.node
    branch = spine.parent.parent
    sim.modify_rec(0, branch, loc=loc)
    sim.parameters['input_loc'] = branch.type
    sim.modify_rec(1, spine)
    syn.source.play(spike_times)
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, syn_index*2)
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    sim.parameters['stim_loc'] = 'branch'
    syn = branch_syn_list[syn_index]
    syn.source.play(spike_times)
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, syn_index*2+1)
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    print 'Process:', os.getpid(), 'completed Iteration:', syn_index, 'Spine:', spine.index, 'Node:', \
        branch.name, 'in', time.time() - start_time, 's'
    return rec_filename


equilibrate = 250.  # time to steady-state
duration = 350.
v_init = -67.
amp = 0.03
syn_type = 'EPSC'  # 'AMPA_KIN'

spine_syn_list = []
branch_syn_list = []
cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.zero_na()

random.seed(0)
for branch in cell.basal+cell.trunk+cell.apical+cell.tuft:
    node_list = []
    if len(branch.spines) > 1:
        if branch.sec.L <= 10.:
            node = branch.spines[random.sample(range(0, len(branch.spines)), 1)[0]]
            node_list.append(node)
        else:
            num_syns = min(len(branch.spines), int(branch.sec.L//10.))  # a random synapse every 10 um
            for i in random.sample(range(0, len(branch.spines)), num_syns):
                node = branch.spines[i]
                node_list.append(node)
    elif branch.spines:
        node = branch.spines[0]
        node_list.append(node)
    for node in node_list:
        syn = Synapse(cell, node, [syn_type], stochastic=0)
        syn.netcon(syn_type).weight[0] = amp
        spine_syn_list.append(syn)
        loc = syn.loc
        syn = Synapse(cell, branch, [syn_type], stochastic=0, loc=loc)
        syn.netcon(syn_type).weight[0] = amp
        branch_syn_list.append(syn)

sim = QuickSim(duration, verbose=False)
sim.parameters['amp'] = amp
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
sim.append_rec(cell, cell.tree.root, 0.5, description='branch')
sim.append_rec(cell, cell.tree.root, 0.5, description='spine')

spike_times = h.Vector([equilibrate])