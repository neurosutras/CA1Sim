__author__ = 'Aaron D. Milstein'
from specify_cells import *
from plot_results import *
from scipy import integrate
import random
"""
Pick num_stim random synapses to stimulate at 5 x 100 Hz, with and without NMDARs. Synapses are either trunk+distal
apical or proximal tuft alone by default.
"""

#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '041315 kap_kad_ih_ampar_scale nmda kd pas no_na.pkl'
rec_filename = '041315 test nmda_ampa_ratio - apical - EB2'

equilibrate = 150  # time to steady-state
duration = 400
v_init = -65.
num_stim = 40
spike_times = h.Vector([equilibrate+i*10 for i in range(5)])

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

sim = QuickSim(duration)  # , verbose=False)
trunk = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type in ['trunk','tuft'] and
                                            trunk.children[1].type in ['trunk', 'tuft']][0]  # tuft bifurcation
sim.append_rec(cell, trunk, description='trunk', loc=0.)
syn_types = ['AMPA_KIN', 'NMDA_KIN']

syn_list = []
border_loc = cell.get_distance_to_node(cell.tree.root, trunk, loc=1.)

for branch in cell.trunk+cell.tuft:  # apical:  # cell.tuft:
    for spine in branch.spines:
        syn = Synapse(cell, spine, syn_types, stochastic=0)
        if (branch.type == 'trunk' and
            border_loc - cell.get_distance_to_node(cell.tree.root, branch, loc=syn.loc) <= 100.) or \
            (branch.type == 'apical' and
            border_loc - cell.get_distance_to_node(cell.tree.root, cell.get_dendrite_origin(branch), 1.) <= 100.) or \
            (branch.type == 'tuft' and
            cell.get_distance_to_node(cell.get_dendrite_origin(branch), branch, loc=syn.loc) <= 100.):
                if branch.type == 'tuft':  # remove this line when probing trunk + apical synapses
                    syn_list.append(syn)
cell.init_synaptic_mechanisms()

random.seed(0)
#stim_syn_list = [syn_list[i] for i in random.sample(range(len(syn_list)), num_stim)]
stim_syn_list = [syn_list[i] for i in random.sample(range(len(syn_list)), 5)]
for i, syn in enumerate(stim_syn_list):
    syn.source.play(h.Vector([equilibrate+i*10]))
    sim.append_rec(cell, syn.node.parent.parent, loc=syn.loc, description='branch')
    sim.append_rec(cell, syn.node, description='spine')
    sim.append_rec(cell, syn.node, param='_ref_i', object=syn.target('AMPA_KIN'), description='i_AMPA',
                   ylabel='Current (nA)')
    sim.append_rec(cell, syn.node, param='_ref_i', object=syn.target('NMDA_KIN'), description='i_NMDA',
                   ylabel='Current (nA)')

sim.parameters['description'] = 'AMPA_KIN + NMDA_KIN'
sim.run(v_init)
areas = []
t = np.array(sim.tvec)
left, right = time2index(t, equilibrate-3., equilibrate-1.)
vm = np.array(sim.rec_list[0]['vec'])
baseline = np.average(vm[left:right])
vm -= baseline
vm[np.where(vm < 0)] = 0
areas.append(integrate.trapz(vm[right:], t[right:]))
with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
    sim.export_to_file(f, 0)
"""
sim.parameters['description'] = 'AMPA_KIN alone'
for syn in syn_list:
    syn.target('NMDA_KIN').gmax = 0.
sim.run(v_init)
t = np.array(sim.tvec)
left, right = time2index(t, equilibrate-3., equilibrate-1.)
vm = np.array(sim.rec_list[0]['vec'])
baseline = np.average(vm[left:right])
vm -= baseline
vm[np.where(vm < 0)] = 0
areas.append(integrate.trapz(vm[right:], t[right:]))
print 'AMPA_KIN + NMDA_KIN:', areas[0], 'AMPA_KIN alone:', areas[1], '% Area (NMDA):', (areas[0]-areas[1])/areas[0]
with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
    sim.export_to_file(f, 1)
"""
plot_superimpose_conditions(rec_filename, legend=False)