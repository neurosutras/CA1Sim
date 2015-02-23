__author__ = 'Aaron D. Milstein'
import time
from plot_results import *
from specify_cells import *
from function_lib import *
from neuron import h
"""
This simulation saves output to a file, and then calls a plot function on that file.
A non-stochastic synapse is stimulated, and different synaptic mechanisms are iterated through in order to compare
their time course and amplitude.
"""

morph_filename = 'MC120914100xC3-scaled.swc'
mech_filename = '022315 kap_scale kd ih_scale no_na.pkl'
rec_filename = '022315 Calibrate_NMDA_D'

amp, equilibrate, duration = 1., 200., 500

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=False)
spike_times = h.Vector([equilibrate])

sim = QuickSim(duration)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
branch = cell.get_node_by_distance_to_soma(200, 'apical')
trunk = cell._get_node_along_path_to_root(branch, 'trunk')
cell.insert_spine(branch, 0.5)
cell._reinit_mech(branch.spines, 1)
head = branch.spines[0]
syn = Synapse(cell, head, ['EPSC'], stochastic=0)
syn.target('EPSC').tau1 = 3
syn.target('EPSC').tau2 = 70
syn.netcon('EPSC').weight[0] = amp * 0.009
syn.source.play(spike_times)
sim.append_rec(cell, cell.tree.root, 0.5)
sim.append_rec(cell, trunk, 0.5)
sim.append_rec(cell, branch, 0.5)
sim.append_rec(cell, head, 0.5)

f = h5py.File(data_dir+rec_filename+'.hdf5', 'w')

simiter = 0
sim.parameters['description'] = 'EPSC'
sim.run()
sim.export_to_file(f, simiter)

del syn
del head.synapses[0]
syn = Synapse(cell, head, ['NMDA_D'], stochastic=0)
syn.target('NMDA_D').mg = 0.1
syn.source.play(spike_times)

sim.parameters['description'] = 'NMDA_D'
gmax_factors = {1: 2, 2: 1, 3: 0.5}
old_gmax = syn.target('NMDA_D').gmax
for simiter in range(1,4):
    syn.target('NMDA_D').gmax = old_gmax * gmax_factors[simiter]
    sim.parameters['description'] = 'NMDA_D: gmax: '+str(syn.target('NMDA_D').gmax)
    sim.run()
    sim.export_to_file(f, simiter)
f.close()
plot_superimpose_conditions(rec_filename)