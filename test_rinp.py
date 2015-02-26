__author__ = 'milsteina'
import time
from specify_cells import *
from plot_results import *
"""
Iterate through every section, inject hyperpolarizing current to measure input resistance.
"""

#morph_filename = 'Erik_Bloss_CA1_0215_Stitched_Proofread.swc'  # EB1
morph_filename = 'EB022715-stitched-proofread.swc'  # EB2
mech_filename = '030415 kap_kad_scale kd ih_scale na.pkl'
rec_filename = '030415 kap_kad_ih_scale kd na full_spines - EB2 - rinp'


equilibrate = 150.  # time to steady-state
duration = 250.
stim_dur = 100.
amp = -0.1

#cell = HocCell(morph_filename, mech_filename)
cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
sim = QuickSim(duration)
sim.append_rec(cell, cell.tree.root, 0.5)
sim.append_stim(cell, cell.tree.root, 0.5, amp, equilibrate, stim_dur)

f = h5py.File(data_dir+rec_filename+'.hdf5', 'w')
simiter = 0
for node in cell.basal+cell.trunk+cell.apical+cell.tuft:
    sim.modify_rec(0, node)
    sim.modify_stim(node=node)
    start_time = time.time()
    print 'Run: ', simiter, ', Section: ', node.name
    sim.run(-65.)
    print 'Took: ', time.time() - start_time, ' sec'
    sim.export_to_file(f, simiter)
    simiter += 1
f.close()

plot_Rinp(rec_filename)