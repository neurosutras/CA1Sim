__author__ = 'milsteina'
from specify_cells import *
import random
import os
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which synapse to stimulate (all spines), for comparing Expected to Actual depolarization.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'

# exponential ampar conductance gradient applied to trunk; inheritance applied to apical and tuft; constant basal
mech_filename = '050715 pas_exp_scale kdr ka_scale ih_sig_scale ampar_exp_scale nmda - EB2'

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
    sim.parameters['spine_index'] = spine.index
    syn.source.play(spike_times)
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, syn_index)
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    print 'Process:', os.getpid(), 'completed Iteration:', syn_index, 'Spine:', syn.node.index, 'Node:', \
        syn.node.parent.parent.name, 'in %.3f s' % (time.time() - start_time)
    return rec_filename


equilibrate = 250.  # time to steady-state
duration = 450.
v_init = -67.
syn_types = ['AMPA_KIN', 'NMDA_KIN']
syn_list = []
cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

random.seed(0)
for branch in cell.basal+cell.trunk+cell.apical+cell.tuft:
    for spine in branch.spines:
        syn = Synapse(cell, spine, syn_types, stochastic=0)
        syn_list.append(syn)
cell.init_synaptic_mechanisms()

sim = QuickSim(duration, verbose=False)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
sim.append_rec(cell, cell.tree.root, description='soma')

# look for a trunk bifurcation
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

sim.append_rec(cell, trunk, 0., description='trunk')
sim.append_rec(cell, trunk, description='branch')  # placeholder for branch

spike_times = h.Vector([equilibrate])