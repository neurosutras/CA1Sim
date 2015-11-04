__author__ = 'Aaron D. Milstein'
from specify_cells import *
import random
import os
from IPython.parallel.util import interactive
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by a list of indexes
corresponding to which synapses to stimulate.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '103115 interim dendritic excitability ampa nmda_kin3'

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())

NMDA_type = 'NMDA_KIN3'
interp_dt = 0.01

@interactive
def stim_actual(spine_indexes):
    """
    Called by controller, mapped to each engine. Sequentially activates spines specified by a list of indexes and saves
    the resulting output to a file.
    :param spine_indexes: list of int
    :return: str
    """
    sim.parameters['syn_indexes'] = []
    for i, index in enumerate(spine_indexes):
        spine = spine_list[index]
        syn = spine.synapses[0]
        spike_times = h.Vector([equilibrate + 0.3 * i])
        syn.source.play(spike_times)
        sim.parameters['syn_indexes'].append(spine.index)
    start_time = time.time()
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, len(spine_indexes))
    for index in spine_indexes:
        spine = spine_list[index]
        syn = spine.synapses[0]
        syn.source.play(h.Vector())
    print 'Process: %i stimulated %i synapses in %i s' % (os.getpid(), len(spine_indexes),
                                                                          time.time() - start_time)
    return rec_filename


@interactive
def stim_expected(spine_index):
    """
    Called by controller, mapped to each engine. Activates a single spine specified by an index and saves the
    resulting output to a file.
    :param spine_index: int
    :return: str
    """
    spine = spine_list[spine_index]
    syn = spine.synapses[0]
    spike_times = h.Vector([equilibrate])
    syn.source.play(spike_times)
    sim.parameters['spine_index'] = spine.index
    start_time = time.time()
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, spine_index)
    syn.source.play(h.Vector())
    print 'Process: %i stimulated spine: %i in %i s' % (os.getpid(), spine.index,
                                                                        time.time() - start_time)
    return rec_filename


def zero_na():
    """

    """
    for sec_type in ['axon_hill', 'ais']:
        cell.modify_mech_param(sec_type, 'nax', 'gbar', 0.)
    cell.reinitialize_subset_mechanisms('axon', 'nax')
    cell.modify_mech_param('soma', 'nas', 'gbar', 0.)
    for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
        cell.reinitialize_subset_mechanisms(sec_type, 'nas')


equilibrate = 250.  # time to steady-state
duration = 450.
v_init = -67.
syn_types = ['AMPA_KIN', NMDA_type]

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
#zero_na()

# these synapses will not be used, but must be inserted for inheritance of synaptic parameters from trunk
for branch in cell.trunk:
    for spine in branch.spines:
        syn = Synapse(cell, spine, syn_types, stochastic=0)

# choose a distal apical oblique branch that has > 25 spines within 30 um, choose spines near the middle of the branch
spine_list = []
for branch in (apical for apical in cell.apical if cell.get_distance_to_node(cell.tree.root,
            cell.get_dendrite_origin(apical)) >= 100. and
                cell.get_distance_to_node(cell.get_dendrite_origin(apical), apical, loc=1.) >= 80.):
    spine_list = [spine for spine in branch.spines if
                  30. <= cell.get_distance_to_node(cell.get_dendrite_origin(branch), spine, loc=0.) <= 60.]
    if len(spine_list) > 25:
        #print 'branch', branch.name, 'has', len(spine_list), 'spines within 30 um'
        for spine in spine_list:
            syn = Synapse(cell, spine, syn_types, stochastic=0)
        break
if not spine_list:
    raise Exception('Could not find a branch that satisfies the specified requirements.')
cell.init_synaptic_mechanisms()

local_random = random.Random()
local_random.seed(branch.index)
local_random.shuffle(spine_list)
sim = QuickSim(duration, verbose=0)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration
sim.parameters['path_type'] = branch.type
sim.parameters['path_index'] = branch.index

sim.append_rec(cell, cell.tree.root, description='soma')
sim.append_rec(cell, cell.get_dendrite_origin(branch), description='trunk', loc=1.)
loc = np.median([spine.synapses[0].loc for spine in spine_list])
sim.append_rec(cell, branch, loc=loc, description='branch')
spine = spine_list[0]
sim.append_rec(cell, spine, object=spine.synapses[0].target(NMDA_type), param='_ref_g', description='NMDA_g')