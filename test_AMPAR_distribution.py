__author__ = 'milsteina'
from specify_cells import *
from plot_results import *

#morph_filename = 'Erik_Bloss_CA1_0215_Stitched_Proofread.swc'  # EB1
morph_filename = 'EB022715-stitched-proofread.swc'  # EB2
#mech_filename = '030515 kap_kad_ih_ampar_scale kd with_na.pkl'
mech_filename = '031215 kap_kad_ih_ampar_scale kd no_na.pkl'

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

for node in cell.get_nodes_of_subtype('spine_head'):
    syn = Synapse(cell, node, ['AMPA_KIN'], stochastic=0)

cell.init_synaptic_mechanisms()
plot_gmax_distribution(cell)
#plot_mech_param_distribution(cell, 'kap', 'gkabar')