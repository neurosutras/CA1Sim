__author__ = 'Aaron D. Milstein'
import time
from specify_cells import *
from plot_results import *

morph_filename = 'EB022715-stitched-proofread.swc'  # EB2: 11064 spines
#morph_filename = 'Erik_Bloss_CA1_0215_Stitched_Proofread.swc'  #EB1:
mech_filename = '030415 kap_kad_ampar_scale kd ih_scale no_na.pkl'
rec_filename = '030415 kap_kad_ih_ampar_scale kd no_na - EB2 - EPSP_attenuation'


equilibrate = 150.  # time to steady-state
duration = 200.

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
for node in cell.get_nodes_of_subtype('spine_head'):
    syn = Synapse(cell, node, ['AMPA_KIN'], stochastic=0)

sim = QuickSim(duration)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
sim.append_rec(cell, cell.tree.root, description='soma')
trunk = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type == 'trunk' and
                                           trunk.children[1].type == 'trunk'][0]  # trunk bifurcation
sim.append_rec(cell, trunk, description='trunk')
sim.append_rec(cell, trunk, description='branch')  # placeholders for branch and spine
sim.append_rec(cell, trunk, description='spine')

spike_times = h.Vector([equilibrate])

f = h5py.File(data_dir+rec_filename+'.hdf5', 'w')
simiter = 0
for node in cell.basal+cell.trunk+cell.apical+cell.tuft:
    sim.modify_rec(2, node)  # sec_type of the branch is the input_loc
    sim.parameters['input_loc'] = node.type
    for spine in node.spines:
        start_time = time.time()
        sim.modify_rec(3, spine)
        syn = spine.synapses[0]
        syn.source.play(spike_times)
        print 'Run: ', simiter, ', AMPAR EPSP'
        sim.run()
        print 'Took: ', time.time() - start_time, ' sec'
        sim.export_to_file(f, simiter)
        syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
        simiter += 1
f.close()

plot_EPSP_attenuation(rec_filename)