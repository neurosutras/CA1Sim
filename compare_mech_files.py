__author__ = 'milsteina'
import specify_cells
from specify_cells import *
from plot_results import *

"""
Specify two cells with different ion channel mechanisms, and probe with current injections to explore features of
spikes and input resistance.
"""

morph_filename = 'EB2-late-bifurcation.swc'
mech_filename1 = '20220808_default_biophysics.yaml'
mech_filename2 = '20260903_altered_soma_biophysics3.yaml'
rec_filename = '20261007 test_compare_biophysics'
rec_file_path = 'data/'+rec_filename+'.hdf5'

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
        print(sim.parameters['description'], ', Rinp: peak: ', peak, ', steady-state: ', steady, ', % sag: ', \
            1-steady/peak)


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
sim.export_to_file(rec_file_path, 0)

del sim
del cell1

cell2 = CA1_Pyr(morph_filename, mech_filename2, full_spines=True)

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

sim.export_to_file(rec_file_path, 1)

plot_superimpose_conditions(rec_filename)
