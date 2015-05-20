__author__ = 'milsteina'
from specify_cells import *
import random
import os
from mpi4py import MPI
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which branch path to stimulate, and a number of synapses to stimulate. Results will be collected to produce an input-
output curve, or as a basis for comparing Expected vs. Actual depolarization.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'

# exponential ampar conductance gradient applied to trunk; inheritance applied to apical and tuft; constant basal
mech_filename = '050715 pas_exp_scale kdr ka_scale ih_sig_scale ampar_exp_scale nmda - EB2'

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())


def stimulate_synapse_group((path_index, num_syns)):
    """
    :param syn_index: int
    :return: str
    """
    start_time = time.time()
    branch = spiny_branches[path_index]
    sim.parameters['path_type'] = branch.type
    sim.parameters['path_index'] = branch.index
    sim.parameters['branch_order'] = cell.get_branch_order(branch)
    if branch.type == 'trunk':
        sim.modify_rec(2, branch)
    else:
        sim.modify_rec(2, cell.get_dendrite_origin(branch))
    sim.parameters['syn_indexes'] = []
    local_random.seed(branch.index)
    stim_syn_indexes = local_random.sample(range(0, len(branch.spines)), len(branch.spines))
    for num, i in enumerate(stim_syn_indexes[:num_syns]):
        syn = branch.spines[i].synapses[0]
        sim.parameters['syn_indexes'].append(branch.spines[i].index)
        # stimulate one spine every 0.3 ms to compare to 2-photon uncaging experiments
        spike_times = h.Vector([equilibrate+num*0.3])
        syn.source.play(spike_times)
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, int(path_index*1e6+num_syns))
    for i in stim_syn_indexes:
        syn = branch.spines[i].synapses[0]
        syn.source.play(h.Vector())     # playing an empty vector turns this synapse off for future runs while keeping
                                        # the VecStim source object in existence so it can be activated again
    print 'Process:', os.getpid(), 'completed Iteration:', path_index, 'Spines:', num_syns, 'Branch:', \
        branch.name, 'in %.3f s' % (time.time() - start_time)
    return rec_filename
    #return {'pid': os.getpid(), 'rank': MPI.COMM_WORLD.Get_rank(), 'path':path_index, 'num_syns': num_syns}


equilibrate = 250.  # time to steady-state
duration = 450.
v_init = -67.
syn_types = ['AMPA_KIN', 'NMDA_KIN']

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

# Here I create a data structure containing groups of spines for which I will construct an input-output curve. Each
# unbranched stretch of dendrite that contains >= 10 spines is included. This is considered "clustered" input, as
# opposed to "distributed" input occurring anywhere along a path from terminal branch to dendrite origin.
# This is a first pass attempt to calibrate within-branch cooperativity, and may require further tweaking.

spiny_branches = [branch for branch in cell.basal+cell.trunk+cell.apical+cell.tuft if len(branch.spines) >= 10]

trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
if trunk_bifurcation:
    trunk_branches = [branch for branch in trunk_bifurcation[0].children if branch.type == 'trunk']
    # get where the thickest trunk branch gives rise to the tuft
    trunk = max(trunk_branches, key=lambda node: node.sec(0.).diam)
    trunk = (node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type
                                                                                    for child in node.children)).next()
else:
    trunk_bifurcation = [node for node in cell.trunk if 'tuft' in (child.type for child in node.children)]
    trunk = trunk_bifurcation[0]

for branch in spiny_branches:
    for spine in branch.spines:
        syn = Synapse(cell, spine, syn_types, stochastic=0)
cell.init_synaptic_mechanisms()

sim = QuickSim(duration, verbose=False)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
sim.append_rec(cell, cell.tree.root, description='soma')
sim.append_rec(cell, trunk, 0., description='trunk')
sim.append_rec(cell, trunk, 1., description='branch')  # placeholder for dendrite branch origin

spike_times = h.Vector([equilibrate])
local_random = random.Random()