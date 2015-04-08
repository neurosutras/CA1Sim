__author__ = 'Aaron D. Milstein'
from specify_cells import *
from plot_results import *
import random

#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '040815 kap_kad_ih_ampar_scale nmda kd pas no_na.pkl'
rec_filename = '040815 kap_kad_ih_ampar_scale nmda kd pas no_na - EB1 - spine AR'


equilibrate = 200.  # time to steady-state
duration = 250.
amp = 0.06
syn_type = 'AMPA_KIN'

spine_syn_list = []
branch_syn_list = []
cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

random.seed(0)
for branch in cell.basal+cell.trunk+cell.apical+cell.tuft:
    node_list = []
    if len(branch.spines) > 1:
        if branch.sec.L <= 10.:
            node = branch.spines[random.sample(range(0, len(branch.spines)), 1)[0]]
            node_list.append(node)
        else:
            num_syns = min(len(branch.spines), int(branch.sec.L//10.))  # a random synapse every 10 um
            for i in random.sample(range(0, len(branch.spines)), num_syns):
                node = branch.spines[i]
                node_list.append(node)
    elif branch.spines:
        node = branch.spines[0]
        node_list.append(node)
    for node in node_list:
        syn = Synapse(cell, node, [syn_type], stochastic=0)
        syn.netcon(syn_type).weight[0] = amp
        spine_syn_list.append(syn)
        loc = syn.loc
        syn = Synapse(cell, branch, [syn_type], stochastic=0, loc=loc)
        syn.netcon(syn_type).weight[0] = amp
        branch_syn_list.append(syn)

sim = QuickSim(duration)
sim.parameters['amp'] = amp
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
sim.append_rec(cell, cell.tree.root, 0.5)
sim.append_rec(cell, cell.tree.root, 0.5)

spike_times = h.Vector([equilibrate])

f = h5py.File(data_dir+rec_filename+'.hdf5', 'w')
simiter = 0
for node in cell.basal+cell.trunk+cell.apical+cell.tuft:
    start_time = time.time()
    sim.parameters['stim_loc'] = 'spine'
    sim.modify_rec(0, node, description='branch')
    sim.modify_rec(1, node.spines[0], description='spine')
    syn = node.spines[0].synapses[0]
    syn.source.play(spike_times)
    print 'Run: ', simiter, ', stim spine'
    sim.run()
    print 'Took: ', time.time() - start_time, ' sec'
    sim.export_to_file(f, simiter)
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    simiter += 1
    start_time = time.time()
    sim.parameters['stim_loc'] = 'branch'
    syn = node.synapses[0]
    syn.source.play(spike_times)
    print 'Run: ', simiter, ', stim branch'
    sim.run()
    print 'Took: ', time.time() - start_time, ' sec'
    sim.export_to_file(f, simiter)
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    simiter += 1
f.close()

plot_AR(rec_filename)