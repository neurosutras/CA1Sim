__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import random

morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'

NMDA_type = 'NMDA_KIN5'
excitatory_stochastic = 1
syn_types = ['AMPA_KIN', NMDA_type]

local_random = random.Random()

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.set_terminal_branch_na_gradient()
cell.insert_inhibitory_synapses_in_subset()

# place synapses in every spine
for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
    for node in cell.get_nodes_of_subtype(sec_type):
        for spine in node.spines:
            syn = Synapse(cell, spine, syn_types, stochastic=excitatory_stochastic)
cell.init_synaptic_mechanisms()

AMPAR_Po = 0.3252

plot_mech_param_distribution(cell, 'pas', 'g', param_label='Leak conductance gradient', svg_title='051316 - cell107')
plot_mech_param_distribution(cell, 'h', 'ghbar', param_label='Ih conductance gradient', svg_title='051316 - cell107')
plot_mech_param_distribution(cell, 'nas', 'gbar', param_label='Na conductance gradient', svg_title='051316 - cell107')
plot_sum_mech_param_distribution(cell, [('kap', 'gkabar'), ('kad', 'gkabar')],
                                 param_label='A-type K conductance gradient', svg_title='051316 - cell107')
plot_synaptic_param_distribution(cell, 'AMPA_KIN', 'gmax', scale_factor=AMPAR_Po,
                                 param_label='Synaptic AMPAR gradient', svg_title='051316 - cell107')