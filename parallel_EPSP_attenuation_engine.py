__author__ = 'milsteina'
from specify_cells import *
import random
import os
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which synapse to stimulate (coarse sampling of the full set of spines).
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'

# exponential ampar conductance gradient applied to trunk; inheritance applied to apical and tuft; constant basal
#mech_filename = '043015 pas_exp_scale kdr ka_scale ih_sig_scale ampar_exp_scale - EB2'
#mech_filename = '071715 rebalanced na_kap_kdr_pas_h - EB2 - spines'
#mech_filename = '080615 rebalanced na_ka ampa nmda - EB2'
#mech_filename = '102915 interim dendritic excitability'
mech_filename = '103015 interim dendritic excitability ampa nmda'

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())


def zero_na():
    """

    """
    for sec_type in ['axon_hill', 'ais']:
        cell.modify_mech_param(sec_type, 'nax', 'gbar', 0.)
    cell.reinitialize_subset_mechanisms('axon', 'nax')
    cell.modify_mech_param('soma', 'nas', 'gbar', 0.)
    for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
        cell.reinitialize_subset_mechanisms(sec_type, 'nas')


def zero_h():
    """

    """
    cell.modify_mech_param('soma', 'h', 'ghbar', 0.)
    cell.mech_dict['trunk']['h']['ghbar']['value'] = 0.
    cell.mech_dict['trunk']['h']['ghbar']['slope'] = 0.
    for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
        cell.reinitialize_subset_mechanisms(sec_type, 'h')


def offset_vm(node, loc, index):
    """


    """
    dt = 0.02
    sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    offset = True
    global i_holding
    count = 0
    while offset and count < 10:
        count += 1
        sim.modify_stim(0, node=node, loc=loc, amp=i_holding)
        sim.run(v_init)
        rec = sim.rec_list[index]['vec']
        vm = np.interp(t, sim.tvec, rec)
        v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
        if v_rest < v_init - 1.:
            i_holding += 0.01
            if sim.verbose:
                print 'increasing i_holding to %.3f' % (i_holding)
        elif v_rest > v_init + 1.:
            i_holding -= 0.01
            if sim.verbose:
                print 'decreasing i_holding to %.3f' % (i_holding)
        else:
            offset = False
    sim.tstop = duration


def stimulate_single_synapse(syn_index):
    """
    :param syn_index: int
    :return: str
    """
    start_time = time.time()
    syn = syn_list[syn_index]
    spine = syn.node
    branch = spine.parent.parent
    sim.modify_rec(2, branch)
    sim.parameters['input_loc'] = branch.type
    sim.modify_rec(3, spine)
    syn.source.play(spike_times)
    offset_vm(branch, syn.loc, 2)
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, syn_index)
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    print 'Process:', os.getpid(), 'completed Iteration:', syn_index, 'Spine:', syn.node.index, 'Node:', \
        syn.node.parent.parent.name, 'in %.3f s' % (time.time() - start_time)
    return rec_filename


equilibrate = 250.  # time to steady-state
duration = 450.
v_init = -67.
#syn_type = 'AMPA_KIN'
NMDA_type = 'NMDA_KIN3'  # 'NMDA_Mg'
syn_types = ['AMPA_KIN', NMDA_type]
#syn_types = ['AMPA_KIN']
local_random = random.Random()
i_holding = 0.

syn_list = []
cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

zero_na()
#zero_h()

local_random.seed(0)
for branch in cell.basal+cell.trunk+cell.apical+cell.tuft:
    if len(branch.spines) > 1:
        if branch.sec.L <= 10.:
            node = branch.spines[local_random.sample(range(0, len(branch.spines)), 1)[0]]
            #syn = Synapse(cell, node, [syn_type], stochastic=0)
            syn = Synapse(cell, node, syn_types, stochastic=0)
            syn_list.append(syn)
        else:
            num_syns = min(len(branch.spines), int(branch.sec.L//10.))  # a random synapse every 10 um
            for i in local_random.sample(range(0, len(branch.spines)), num_syns):
                node = branch.spines[i]
                #syn = Synapse(cell, node, [syn_type], stochastic=0)
                syn = Synapse(cell, node, syn_types, stochastic=0)
                syn_list.append(syn)
    elif branch.spines:
        node = branch.spines[0]
        #syn = Synapse(cell, node, [syn_type], stochastic=0)
        syn = Synapse(cell, node, syn_types, stochastic=0)
        syn_list.append(syn)
cell.init_synaptic_mechanisms()

sim = QuickSim(duration)  # , verbose=False)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
sim.append_rec(cell, cell.tree.root, description='soma')

# look for a trunk bifurcation
trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
if trunk_bifurcation:
    trunk_branches = [branch for branch in trunk_bifurcation[0].children if branch.type == 'trunk']
    # get where the thickest trunk branch gives rise to the tuft
    trunk = max(trunk_branches, key=lambda node: node.sec(0.).diam)
    trunk = (node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type
                                                                            for child in node.children)).next()
else:
    trunk_bifurcation = [node for node in cell.trunk if 'tuft' in (child.type for child in node.children)]
    trunk = trunk_bifurcation[0]
tuft = (child for child in trunk.children if child.type == 'tuft').next()
#distal_trunk = trunk
#trunk = trunk_bifurcation[0]

sim.append_rec(cell, trunk, description='trunk')
sim.append_rec(cell, trunk, description='branch')  # placeholders for branch and spine
sim.append_rec(cell, trunk, description='spine')
sim.append_stim(cell, trunk, loc=0.5, amp=0., dur=duration, delay=0., description='branch')

spike_times = h.Vector([equilibrate])

