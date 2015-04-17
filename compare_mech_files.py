__author__ = 'milsteina'
import specify_cells
from specify_cells import *
from plot_results import *

"""
Specify two cells with different ion channel mechanisms, and probe with current injections to explore features of
spikes and input resistance.
"""

#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename1 = '030515 kap_kad_ih_ampar_scale kd no_na.pkl'
mech_filename2 = '030515 kap_kad_ih_ampar_scale kd with_na.pkl'
rec_filename = '030515 compare_spines'

#mech_dict = read_from_pkl(data_dir+mech_filename)

equilibrate = 150.  # time to steady-state
stim_dur = 100.
duration = 300.
amp = 1.5 # -0.1
v_rest = -65.  # -66.66


def update_amp(sim, amp):
    for i in range(len(sim.stim_list)):
        sim.modify_stim(i, amp=amp)
    sim.run()
    sim.plot()


def print_Rinp(sim):
    for rec in sim.rec_list:
        v_rest, peak, steady = get_Rinp(np.array(sim.tvec), np.array(rec['vec']), equilibrate, equilibrate+stim_dur,
                                        amp)
        print sim.parameters['description'], ', Rinp: peak: ', peak, ', steady-state: ', steady, ', % sag: ', \
            1-steady/peak


cell1 = CA1_Pyr(morph_filename, mech_filename1, full_spines=True)

sim = QuickSim(duration)
sim.parameters['description'] = 'full spines, no na'
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
# trunk bifurcation
node = [trunk for trunk in cell1.trunk if len(trunk.children) > 1 and trunk.children[0].type == 'trunk' and
                                                                        trunk.children[1].type == 'trunk'][0]
#node = cell1.tree.root
sim.append_stim(cell1, node, 0.5, amp, equilibrate, stim_dur)
sim.append_rec(cell1, node, 0.5, description='soma')
sim.run(v_rest)
print_Rinp(sim)
f = h5py.File(data_dir+rec_filename+'.hdf5', 'w')
sim.export_to_file(f, 0)

del sim
del cell1

cell2 = CA1_Pyr(morph_filename, mech_filename2, full_spines=True)
"""
cell2.modify_mech_param('soma', 'hh3na', 'gbar', 0.04)
cell2.modify_mech_param('axon', 'hh3nax', 'gbar', 0.2)
cell2.modify_mech_param('ais', 'hh3nax', 'gbar', 1.)
for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
    cell2.modify_mech_param(sec_type, 'hh3na', 'gbar', origin='soma')
"""

sim = QuickSim(duration)
sim.parameters['description'] = 'full spines, with na'
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
node2 = [trunk for trunk in cell2.trunk if len(trunk.children) > 1 and trunk.children[0].type == 'trunk' and
                                                                        trunk.children[1].type == 'trunk'][0]
#node2 = cell2.tree.root
sim.append_stim(cell2, node2, 0.5, amp, equilibrate, stim_dur)
sim.append_rec(cell2, node2, 0.5, description='soma')
sim.run(v_rest)
print_Rinp(sim)
#f = h5py.File(data_dir+rec_filename+'.hdf5', 'w')
sim.export_to_file(f, 1)
f.close()

plot_superimpose_conditions(rec_filename)

"""
for node in cell.trunk+cell.basal+cell.apical+cell.tuft:
    for seg in node.sec:
        distance = cell.get_distance_to_node(cell.tree.root, node) + seg.x * node.sec.L
        if hasattr(seg, 'kap'):
            print node.name, 'distance: ', distance, ', kap: ', seg.kap.gkabar
        if hasattr(seg, 'kad'):
            print node.name, 'distance: ', distance, ', kad: ', seg.kad.gkabar
"""