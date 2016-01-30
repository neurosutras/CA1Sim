__author__ = 'milsteina'
from specify_cells import *
import random
import os
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by a list of indexes
corresponding to which synapses to stimulate. Remember to categorize output by distance from dendrite origin to soma.
"""
morph_filename = 'EB2-late-bifurcation.swc'
# mech_filename = '052915 pas_exp_scale kdr ka_scale ih_sig_scale ampar_exp_scale nmda - EB2'
# mech_filename = '080615 rebalanced na_ka ampa nmda - EB2'
# mech_filename = '103115 interim dendritic excitability ampa nmda_kin3'
mech_filename = '012816 altered intrinsic properties - ampa nmda_kin4'
rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())

NMDA_type = 'NMDA_KIN4'
ISI = 0.3


def stim_actual_group((group_index, num_spines)):
    """
    Called by controller, mapped to each engine. Activates the specified number of spines in the specified group.
    :param group_index: int
    :param num_spines: int
    :return: str
    """
    spine_group = groups_to_stim[group_index]
    loc = np.median([spine.synapses[0].loc for spine in spine_group['spines']])
    path_type = spine_group['path_type']
    path_index = spine_group['path_index']
    sim.parameters['path_index'] = path_index
    sim.parameters['path_type'] = path_type
    sim.parameters['path_category'] = spine_group['path_category']
    spine = spine_group['spines'][0]
    sim.modify_rec(2, node=spine.parent.parent, loc=loc)
    if path_type == 'trunk':
        sim.modify_rec(3, node=spine.parent.parent)
    else:
        sim.modify_rec(3, node=cell.get_dendrite_origin(spine))
    sim.modify_rec(4, node=spine, object=spine.synapses[0].target(NMDA_type), param='_ref_g')
    sim.parameters['syn_indexes'] = []
    for i, spine in enumerate(spine_group['spines'][:num_spines]):
        sim.parameters['syn_indexes'].append(spine.index)
        syn = spine.synapses[0]
        spike_times = h.Vector([equilibrate + ISI * i])
        syn.source.play(spike_times)
    start_time = time.time()
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, int(group_index*1e6+num_spines))
    print 'Process: %i took %i s to stimulate %i synapses in path %i' % (os.getpid(), time.time() - start_time,
                                                                         num_spines, path_index)
    for spine in spine_group['spines'][:num_spines]:
        syn = spine.synapses[0]
        syn.source.play(h.Vector())
    return rec_filename


def stim_single_expected((group_index, spine_index)):
    """
    Called by controller, mapped to each engine. Activates a single spine specified by an index and saves the
    resulting output to a file.
    :param group_index: int
    :param spine_index: int
    :return: str
    """
    spine_group = groups_to_stim[group_index]
    path_type = spine_group['path_type']
    loc = np.median([spine.synapses[0].loc for spine in spine_group['spines']])
    spine = spine_group['spines'][spine_index]
    syn = spine.synapses[0]
    spike_times = h.Vector([equilibrate])
    syn.source.play(spike_times)
    sim.parameters['spine_index'] = spine.index
    sim.parameters['path_index'] = spine_group['path_index']
    sim.parameters['path_type'] = path_type
    sim.parameters['path_category'] = spine_group['path_category']
    sim.modify_rec(2, node=spine.parent.parent, loc=loc)
    if path_type == 'trunk':
        sim.modify_rec(3, node=spine.parent.parent)
    else:
        sim.modify_rec(3, node=cell.get_dendrite_origin(spine))
    sim.modify_rec(4, node=spine, object=syn.target(NMDA_type), param='_ref_g')
    start_time = time.time()
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, int(group_index*1e6+spine_index))
    syn.source.play(h.Vector())
    print 'Process: %i stimulated spine: %i in %i s' % (os.getpid(), spine.index, time.time() - start_time)
    return rec_filename


equilibrate = 250.  # time to steady-state
duration = 450.
v_init = -65.
syn_types = ['AMPA_KIN', NMDA_type]

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

cell.zero_na()

# these synapses will not be used, but must be inserted for inheritance of synaptic parameters from trunk
for branch in cell.trunk:
    for spine in branch.spines:
        syn = Synapse(cell, spine, syn_types, stochastic=0)

min_num_spines = 25
grouped_spines = {sec_type: [] for sec_type in ['basal', 'apical', 'tuft']}
for sec_type in ['basal', 'apical', 'tuft']:
    path_category = 'proximal'
    for branch in (branch for branch in cell._node_dict[sec_type] if cell.get_branch_order(branch) == 1 and
                    branch.sec.L >= 30. and len(branch.spines) >= min_num_spines):
        spine_list = [spine for spine in branch.spines if cell.get_distance_to_node(cell.get_dendrite_origin(branch),
                                                                                    spine, loc=0.) <= 30.]
        grouped_spines[sec_type].append({'path_category': path_category, 'spines': spine_list,
                                    'path_index': spine_list[0].index, 'path_type': sec_type})
    path_category = 'terminal'
    for branch in (branch for branch in cell._node_dict[sec_type] if cell.is_terminal(branch) and
                    branch.sec.L >= 30. and len(branch.spines) >= min_num_spines):
        distance = cell.get_distance_to_node(cell.get_dendrite_origin(branch), branch, loc=1.)
        spine_list = [spine for spine in branch.spines if distance -
                      cell.get_distance_to_node(cell.get_dendrite_origin(branch), spine, loc=0.) <= 30.]
        grouped_spines[sec_type].append({'path_category': path_category, 'spines': spine_list,
                                    'path_index': spine_list[0].index, 'path_type': sec_type})
    path_category = 'intermediate'
    for branch in (branch for branch in cell._node_dict[sec_type] if
                   cell.get_distance_to_node(cell.get_dendrite_origin(branch), branch, loc=1.) >= 80.):
        spine_list = [spine for spine in branch.spines if
                  30. <= cell.get_distance_to_node(cell.get_dendrite_origin(branch), spine, loc=0.) <= 60.]
        if len(spine_list) >= min_num_spines:
            grouped_spines[sec_type].append({'path_category': path_category, 'spines': spine_list,
                                    'path_index': spine_list[0].index, 'path_type': sec_type})

local_random = random.Random()
groups_to_stim = []
for sec_type in ['basal', 'apical', 'tuft']:
    groups_to_stim.extend(grouped_spines[sec_type])

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

for spine_group in groups_to_stim:
    local_random.seed(spine_group['path_index'])
    local_random.shuffle(spine_group['spines'])
    for spine in spine_group['spines']:
        syn = Synapse(cell, spine, syn_types, stochastic=0)
cell.init_synaptic_mechanisms()

sim = QuickSim(duration, verbose=0)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration

sim.append_rec(cell, cell.tree.root, description='soma', loc=0.)
sim.append_rec(cell, trunk, description='trunk', loc=0.)
sim.append_rec(cell, trunk, description='branch', loc=0.5)  # placeholder for local branch
sim.append_rec(cell, trunk, description='origin', loc=1.)  # placeholder for local dendrite origin
spine = groups_to_stim[0]['spines'][0]
sim.append_rec(cell, spine, object=spine.synapses[0].target(NMDA_type), param='_ref_g', description='NMDA_g')