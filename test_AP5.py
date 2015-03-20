__author__ = 'Aaron D. Milstein'
from specify_cells import *
from plot_results import *
from scipy import integrate
import random
"""
Pick num_stim random synapses along the trunk to stimulate at 5 x 100 Hz, with and without NMDARs.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '031815 kap_kad_ih_ampar_scale nmda kd no_na.pkl'  # 7.7e-3
rec_filename = 'calibrate_nmda_gmax'

equilibrate = 150  # time to steady-state
duration = 550
v_init = -65.
num_stim = 40
#spike_times = h.Vector([equilibrate + delay for delay in range(0, 50, 10)])
spike_times = h.Vector([equilibrate])

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
"""
cell.modify_mech_param('soma', 'h', 'ghbar', 0.)
cell.modify_mech_param('trunk', 'h', 'ghbar', 0.)
for nodes in [cell.basal, cell.apical, cell.tuft]:
    cell._reinit_mech(nodes)
cell.modify_mech_param('trunk', 'synapse', 'mg', 0.1, syn_type='NMDA_KIN')
cell.modify_mech_param('apical', 'synapse', 'mg', 0.1, syn_type='NMDA_KIN')

# this NMDA_KIN.gmax matches ~35% NMDA/AMPA ratio in 0.1 mM Mg and ZD from Bittner et al., 2012
# but recruits very little NMDA under normal conditions, even with large depolarizations...
"""
cell.mech_dict['trunk']['synapse']['gmax'][2]['value'] = 1.45e-4
#cell.mech_dict['trunk']['synapse']['gmax'][2]['syn_type'] = 'NMDA_KIN'
#cell.mech_dict['apical']['synapse']['gmax'][1]['syn_type'] = 'NMDA_KIN'
#cell.mech_dict['basal']['synapse']['gmax'][2]['syn_type'] = 'NMDA_KIN'
cell.mech_dict['basal']['synapse']['gmax'][2]['value'] = 1.45e-4
#cell.mech_dict['tuft']['synapse']['gmax'][1]['syn_type'] = 'NMDA_KIN'
cell.mech_dict['tuft']['synapse']['gmax'][1]['value'] = 2.85e-3

sim = QuickSim(duration)  # , verbose=False)
trunk = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type == 'trunk' and
                                               trunk.children[1].type == 'trunk'][0]  # trunk bifurcation
sim.append_rec(cell, trunk, description='trunk')
syn_types = ['AMPA_KIN', 'NMDA_KIN']
spine_list = []

for branch in cell.trunk+cell.tuft:  # cell.apical
    spine_list.extend(branch.spines)

for spine in spine_list:
    syn = Synapse(cell, spine, syn_types, stochastic=0)
cell.init_synaptic_mechanisms()

spine_list = []
for branch in cell.tuft:
    spine_list.extend(branch.spines)

random.seed(0)
stim_syn_list = [spine_list[i].synapses[0] for i in random.sample(range(len(spine_list)), num_stim)]
for syn in stim_syn_list:
    syn.source.play(spike_times)

#sim.append_rec(cell, stim_syn_list[5].node, param='_ref_Rb', object=stim_syn_list[5].target('NMDA_KIN'))

sim.parameters['description'] = 'AMPA_KIN + NMDA_KIN'
sim.run(v_init)
areas = []
t = np.array(sim.tvec)
left, right = time2index(t, equilibrate-2.0, equilibrate)
vm = np.array(sim.rec_list[0]['vec'])
baseline = np.average(vm[left:right])
vm -= baseline
depolarized = np.where(vm[right:] > 0)[0]
start = depolarized[0] + right
end = depolarized[-1] + right
areas.append(integrate.trapz(vm[start:end], t[start:end]))
with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
    sim.export_to_file(f, 0)
sim.parameters['description'] = 'AMPA_KIN alone'
for syn in stim_syn_list:
    syn.target('NMDA_KIN').gmax = 0
sim.run(v_init)
t = np.array(sim.tvec)
left, right = time2index(t, equilibrate-2.0, equilibrate)
vm = np.array(sim.rec_list[0]['vec'])
baseline = np.average(vm[left:right])
vm -= baseline
depolarized = np.where(vm[right:] > 0)[0]
start = depolarized[0] + right
end = depolarized[-1] + right
areas.append(integrate.trapz(vm[start:end], t[start:end]))
print 'AMPA_KIN + NMDA_KIN:', areas[0], 'AMPA_KIN alone:', areas[1], '% Area (NMDA):', (areas[0]-areas[1])/areas[0]
with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
    sim.export_to_file(f, 1)
plot_superimpose_conditions(rec_filename)