__author__ = 'milsteina'
from specify_cells import *
import random
import os
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which synapse to stimulate (coarse sampling of the full set of spines).
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'

# linear ampar conductance gradient applied to trunk and basal, inheritance applied to apical and tuft
#mech_filename = '032315 kap_kad_ampar_scale nmda kd pas no_ih no_na.pkl'

# constant ampar conductance
#mech_filename = '032315 kap_kad_ih_scale kd pas no_na.pkl'

# exponential ampar conductance gradient applied to trunk and basal, inheritance applied to apical and tuft
mech_filename = '040815 kap_kad_ih_ampar_scale nmda kd pas no_na.pkl'

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())


def stimulate_single_synapse(syn_index):
    """
    :param syn_index: int
    :return: str
    """
    start_time = time.time()
    syn = syn_list[syn_index]
    spine = syn.node
    branch = spine.parent.parent
    sim.modify_rec(2, branch)
    sim.parameters['input_loc'] = branch.type
    sim.modify_rec(3, spine)
    syn.source.play(spike_times)
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, syn_index)
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    print 'Process:', os.getpid(), 'completed Iteration:', syn_index, 'Spine:', syn.node.index, 'Node:', \
        syn.node.parent.parent.name, 'in', time.time() - start_time, 's'
    return rec_filename


equilibrate = 150.  # time to steady-state
duration = 200.
v_init = -65.
syn_type = 'AMPA_KIN'

syn_list = []
cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
"""
for node in cell.get_nodes_of_subtype('spine_head'):
    syn = Synapse(cell, node, ['AMPA_KIN'], stochastic=0)
    syn_list.append(syn)
"""
random.seed(0)
for branch in cell.basal+cell.trunk+cell.apical+cell.tuft:
    if len(branch.spines) > 1:
        if branch.sec.L <= 10.:
            node = branch.spines[random.sample(range(0, len(branch.spines)), 1)[0]]
            syn = Synapse(cell, node, [syn_type], stochastic=0)
            syn_list.append(syn)
        else:
            num_syns = min(len(branch.spines), int(branch.sec.L//10.))  # a random synapse every 10 um
            for i in random.sample(range(0, len(branch.spines)), num_syns):
                node = branch.spines[i]
                syn = Synapse(cell, node, [syn_type], stochastic=0)
                syn_list.append(syn)
    elif branch.spines:
        node = branch.spines[0]
        syn = Synapse(cell, node, [syn_type], stochastic=0)
        syn_list.append(syn)
cell.init_synaptic_mechanisms()

sim = QuickSim(duration, verbose=False)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
sim.append_rec(cell, cell.tree.root, description='soma')
trunk = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type == 'trunk' and
                                           trunk.children[1].type == 'trunk'][0]  # trunk bifurcation
sim.append_rec(cell, trunk, description='trunk')
sim.append_rec(cell, trunk, description='branch')  # placeholders for branch and spine
sim.append_rec(cell, trunk, description='spine')

spike_times = h.Vector([equilibrate])

