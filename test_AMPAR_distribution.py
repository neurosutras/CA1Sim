__author__ = 'milsteina'
import specify_cells
from specify_cells import *

morph_filename = 'EB022715-stitched-proofread.swc'  # EB2: 11067 spines
#morph_filename = 'Erik_Bloss_CA1_0215_Stitched_Proofread.swc'  # EB1: 10777 spines
mech_filename = '030515 kap_kad_ih_ampar_scale kd no_na.pkl'

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

for spine in cell.get_nodes_of_subtype('spine_head'):
    syn = Synapse(cell, spine, ['AMPA_KIN'], stochastic=0)

cell.init_synaptic_mechanisms()
"""
cell.modify_mech_param('trunk', 'synapse', 'gmax', 0.004, origin='soma', slope=1.9e-05, syn_type='AMPA_KIN')
cell.modify_mech_param('apical', 'synapse', 'gmax', origin='trunk', syn_type='AMPA_KIN')
cell.modify_mech_param('basal', 'synapse', 'gmax', 0.004, syn_type='AMPA_KIN')
cell.modify_mech_param('tuft', 'synapse', 'gmax', origin='trunk', syn_type='AMPA_KIN')
"""
syn_list = []
distances = []
gmax_vals = []

colors = ['r', 'b', 'g', 'c']
dend_types = ['basal', 'trunk', 'apical', 'tuft']

for i, sec_type in enumerate(dend_types):
    syn_list = []
    distances = []
    gmax_vals = []
    for branch in cell.get_nodes_of_subtype(sec_type):
        for spine in branch.spines:
            syn_list.extend(spine.synapses)
    for syn in syn_list:
        distances.append(cell.get_distance_to_node(cell.tree.root, syn.node.parent.parent, syn.loc))
        gmax_vals.append(syn.target('AMPA_KIN').gmax)
    plt.scatter(distances, gmax_vals, color=colors[i], label=sec_type)
plt.legend(loc='best')
plt.show()
plt.close()

print '# spines: ', len(cell.get_nodes_of_subtype('spine_head'))