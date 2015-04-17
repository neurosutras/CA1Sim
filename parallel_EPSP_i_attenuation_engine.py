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

mech_filename = '042015 soma_pas spines - EB2.pkl'
#mech_filename = '042015 soma_pas kdr ka_scale - adjusted - EB2.pkl'

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())


def stimulate_single_synapse(syn_index):
    """
    :param syn_index: int
    :return: str
    """
    start_time = time.time()
    branch = nodes[syn_index]
    sim.modify_rec(2, branch)
    sim.parameters['input_loc'] = branch.type
    syn = branch.synapses[0]
    syn.source.play(spike_times)
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, syn_index)
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    print 'Process:', os.getpid(), 'completed Iteration:', syn_index, 'Node:', branch.name, 'in', \
        time.time() - start_time, 's'
    return rec_filename


equilibrate = 200.  # time to steady-state
duration = 300.
v_init = -77.
syn_type = 'EPSC'

#cell = CA1_Pyr(morph_filename, mech_filename, full_spines=False)
cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

nodes = cell.soma+cell.basal+cell.trunk+cell.apical+cell.tuft

random.seed(0)
for branch in nodes:
    syn = Synapse(cell, branch, [syn_type], stochastic=0)
    syn.target(syn_type).imax = 0.03

sim = QuickSim(duration, verbose=False)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
sim.append_rec(cell, cell.tree.root, description='soma', loc=0.)
trunk = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type in ['trunk','tuft'] and
                                            trunk.children[1].type in ['trunk', 'tuft']][0]  # tuft bifurcation
sim.append_rec(cell, trunk, description='trunk', loc=0.)
sim.append_rec(cell, trunk, description='branch')  # placeholder for branch

spike_times = h.Vector([equilibrate])