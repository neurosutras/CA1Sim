__author__ = 'Aaron D. Milstein'
import time
from plot_results import *
from specify_cells import *
from function_lib import *
from neuron import h
"""
Specify two cells with different ion channel mechanisms, and probe with current injections to explore features of
spikes and input resistance.
"""
#morph_filename = 'EB1-early-bifurcation.swc'  # 10777 spines
morph_filename = 'EB2-late-bifurcation.swc'  # 11067 spines
#mech_filename = 'mech_dict_021920151042.pkl'
rec_filename = '022315 compare_mech_dicts'

equilibrate = 100.  # time to steady-state
stim_dur = 100.
duration = 250.
amp = -0.1
v_rest = -90.
v_init = -75.


def update_amp(sim, amp):
    for i in range(len(sim.stim_list)):
        sim.modify_stim(i, amp=amp)
    sim.run()
    sim.plot()


def print_Rinp(sim):
    for rec in sim.rec_list:
        v_rest, peak, steady = get_Rinp(np.array(sim.tvec), np.array(rec['vec']), equilibrate, equilibrate+stim_dur,
                                        amp)
        print rec['description'], ', Rinp: peak: ', peak, ', steady-state: ', steady, ', % sag: ', 1-steady/peak


cell = HocCell(morph_filename)
cell.modify_mech_param('soma', 'cable', 'Ra', 150.)
#cell.modify_mech_param('soma', 'hh3na', 'gbar', 0.04)
cell.modify_mech_param('soma', 'kap', 'gkabar', 0.048)
cell.modify_mech_param('soma', 'kdr', 'gkdrbar', 0.04)
cell.modify_mech_param('soma', 'pas', 'g', 1.5e-5)
cell.modify_mech_param('soma', 'pas', 'e', v_init)
cell.modify_mech_param('soma', 'h', 'eh', -30.)
cell.modify_mech_param('soma', 'h', 'ghbar', 5.0e-4)
cell.modify_mech_param('soma', 'h', 'vhalfl', -82.)
cell.modify_mech_param('soma', 'ions', 'ek', -90.)
for sec_type in 'axon', 'ais', 'basal', 'trunk', 'apical', 'tuft':
    cell.modify_mech_param(sec_type, 'cable', 'Ra', origin='soma')
    cell.modify_mech_param(sec_type, 'pas', 'g', origin='soma')
    cell.modify_mech_param(sec_type, 'pas', 'e', origin='soma')
    cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param(sec_type, 'ions', 'ek', origin='soma')
    if sec_type == 'axon':
#        cell.modify_mech_param(sec_type, 'hh3nax', 'gbar', 0.2)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', 0.0096)
    elif sec_type == 'ais':
#        cell.modify_mech_param(sec_type, 'hh3nax', 'gbar', 1.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='axon')
    else:
        cell.modify_mech_param(sec_type, 'h', 'eh', origin='soma')
#        cell.modify_mech_param(sec_type, 'hh3na', 'gbar', origin='soma')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', slope=3.84e-4, max=0.24)
for sec_type in 'spine_head', 'spine_neck':
    cell.modify_mech_param(sec_type, 'cable', 'Ra', origin='soma')
    cell.modify_mech_param(sec_type, 'pas', 'g', origin='soma')
    cell.modify_mech_param(sec_type, 'pas', 'e', origin='soma')
cell.modify_mech_param('basal', 'h', 'ghbar', origin='soma')
cell.modify_mech_param('basal', 'h', 'vhalfl', origin='soma')
cell.modify_mech_param('trunk', 'h', 'ghbar', origin='soma', slope=9.0e-6)
cell.modify_mech_param('trunk', 'h', 'vhalfl', -90.)  # need to implement different values for first 100 um
for sec_type in 'apical', 'tuft':
    cell.modify_mech_param(sec_type, 'h', 'ghbar', origin='trunk')
    cell.modify_mech_param(sec_type, 'h', 'vhalfl', origin='trunk')

sim = QuickSim(duration)
sim.parameters['description'] = 'kap_scale no_na'
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
#node = cell.get_node_by_distance_to_soma(300., 'trunk')
node = cell.tree.root
sim.append_stim(cell, node, 0.5, amp, equilibrate, stim_dur)
sim.append_rec(cell, node, 0.5, description='kap')
sim.run(-65.)
f = h5py.File(data_dir+rec_filename+'.hdf5', 'w')
sim.export_to_file(f, 0)

del sim

cell2 = HocCell(morph_filename)
cell2.modify_mech_param('soma', 'cable', 'Ra', 150.)
cell2.modify_mech_param('soma', 'hh3na', 'gbar', 0.04)
cell2.modify_mech_param('soma', 'hh3k', 'gbar', 0.036)
cell2.modify_mech_param('soma', 'hh3k', 'gbar2', 0.012)
cell2.modify_mech_param('soma', 'pas', 'g', 3.0e-5)
cell2.modify_mech_param('soma', 'pas', 'e', v_rest)
cell2.modify_mech_param('soma', 'h', 'eh', -30.)
cell2.modify_mech_param('soma', 'h', 'ghbar', 5.0e-4)
cell2.modify_mech_param('soma', 'h', 'vhalfl', -82.)
cell2.modify_mech_param('soma', 'ions', 'ek', -90.)
for sec_type in 'axon', 'ais', 'basal', 'trunk', 'apical', 'tuft':
    cell2.modify_mech_param(sec_type, 'cable', 'Ra', origin='soma')
    cell2.modify_mech_param(sec_type, 'pas', 'g', origin='soma')
    cell2.modify_mech_param(sec_type, 'pas', 'e', origin='soma')
    if sec_type == 'axon':
        cell2.modify_mech_param(sec_type, 'hh3nax', 'gbar', 0.2)
        cell2.modify_mech_param(sec_type, 'hh3k', 'gbar', 0.0096)
        cell2.modify_mech_param(sec_type, 'hh3k', 'gbar2', 0.0)
    elif sec_type == 'ais':
        cell2.modify_mech_param(sec_type, 'hh3nax', 'gbar', 1.)
        cell2.modify_mech_param(sec_type, 'hh3k', 'gbar', origin='axon')
        cell2.modify_mech_param(sec_type, 'hh3k', 'gbar2', origin='axon')
    else:
        cell2.modify_mech_param(sec_type, 'hh3na', 'gbar', origin='soma')
        cell2.modify_mech_param(sec_type, 'hh3k', 'gbar', origin='soma')
        cell2.modify_mech_param(sec_type, 'hh3k', 'gbar2', origin='soma')
    cell2.modify_mech_param(sec_type, 'ions', 'ek', origin='soma')
for sec_type in 'spine_head', 'spine_neck':
    cell2.modify_mech_param(sec_type, 'cable', 'Ra', origin='soma')
    cell2.modify_mech_param(sec_type, 'pas', 'g', origin='soma')
    cell2.modify_mech_param(sec_type, 'pas', 'e', origin='soma')
cell2.modify_mech_param('basal', 'h', 'ghbar', origin='soma')
cell2.modify_mech_param('basal', 'h', 'vhalfl', origin='soma')
cell2.modify_mech_param('trunk', 'h', 'ghbar', origin='soma', slope=9.0e-6)
cell2.modify_mech_param('trunk', 'h', 'vhalfl', -90.)  # need to implement different values for first 100 um
for sec_type in 'apical', 'tuft':
    cell2.modify_mech_param(sec_type, 'h', 'ghbar', origin='trunk')
    cell2.modify_mech_param(sec_type, 'h', 'vhalfl', origin='trunk')

sim = QuickSim(duration)
sim.parameters['description'] = 'hh3k'
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
#node = cell2.get_node_by_distance_to_soma(300., 'trunk')
node2 = cell2.tree.root
sim.append_stim(cell2, node2, 0.5, amp, equilibrate, stim_dur)
sim.append_rec(cell2, node2, 0.5, description='kap')
sim.run(-65.)
sim.export_to_file(f, 1)
f.close()

plot_superimpose_conditions(rec_filename)