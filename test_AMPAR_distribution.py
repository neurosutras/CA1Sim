__author__ = 'milsteina'
from specify_cells import *
from plot_results import *

#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '031215 kap_kad_ih_ampar_scale kd no_na.pkl'

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

for node in cell.get_nodes_of_subtype('spine_head'):
    syn = Synapse(cell, node, ['AMPA_KIN'], stochastic=0)

cell.init_synaptic_mechanisms()
plot_synaptic_param_distribution(cell, 'AMPA_KIN', 'gmax')
#plot_mech_param_distribution(cell, 'kap', 'gkabar')