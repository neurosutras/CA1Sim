__author__ = 'milsteina'
import time
from plot_results import *
from specify_cells import *
from function_lib import *
from neuron import h
import scipy as sp
"""
Iterate through every section, inject hyperpolarizing current to measure input resistance.
"""

morph_filename = 'MC120914100xC3-scaled.swc'
mech_filename = '022315 kap_scale kd ih_scale.pkl'
rec_filename = '022315 kap_scale kd ih_scale'


equilibrate = 100.  # time to steady-state
duration = 200.
amp = -0.1

cell = HocCell(morph_filename, mech_filename)
sim = QuickSim(duration)
sim.append_rec(cell, cell.tree.root, 0.5)
sim.append_stim(cell, cell.tree.root, 0.5, amp, equilibrate, stim_dur)

f = h5py.File(data_dir+rec_filename+'.hdf5')
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