__author__ = 'milsteina'
from specify_cells import *
import random
import os
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by a list of indexes
corresponding to which synapses to stimulate. Remember to categorize output by distance from dendrite origin to soma.
"""
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '052915 pas_exp_scale kdr ka_scale ih_sig_scale ampar_exp_scale nmda - EB2'
rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())

NMDA_type = 'NMDA_KIN2'
ISI = 0.3


def get_clustered_spines(cell, branch_origin, spine_list, min_num, length, direction=1):
    """
    Given a list of spines, ordered by distance from branch_origin, returns a group of at least min_num spines
    clustered in a region less than the specified length. Starts either at the beginning or end of the list of spines
    (direction = 0 or 1).
    :param cell: :class:'SHocCell'
    :param branch_origin: :class:'SHocNode'
    :param spine_list: list
    :param min_num: int
    :param length: float
    :param direction: int in [0, 1]
    :return: list
    """
    # do not assign spines to more than one spine group
    sec_type = spine_list[0].parent.parent.type
    used_spines = []
    for spine_group in grouped_spines[sec_type]:
        used_spines.extend(spine_group['spines'])
    used_spines = list(set(used_spines))
    spine_list = [spine for spine in spine_list if not spine in used_spines]
    j = len(spine_list)
    if j < min_num:
        return []
    i = 0
    if direction:
        while cell.get_distance_to_node(branch_origin, spine_list[j-1].parent, 0.) - \
                cell.get_distance_to_node(branch_origin, spine_list[i].parent, 0.) > length:
            i += 1
        if j - i > min_num:
            return spine_list[i:j]
        else:
            return []
    else:
        for i in range(0, len(spine_list)-min_num+1):
            j = len(spine_list)
            while j - i >= min_num:
                if cell.get_distance_to_node(branch_origin, spine_list[j-1].parent, 0.) - \
                        cell.get_distance_to_node(branch_origin, spine_list[i].parent, 0.) > length:
                    j -= 1
                else:
                    return spine_list[i:j]
            if cell.get_distance_to_node(branch_origin, spine_list[-1].parent, 0.) - \
                    cell.get_distance_to_node(branch_origin, spine_list[i+1].parent, 0.) < length:
                return []
        return []


def stim_actual_group((group_index, num_spines)):
    """
    Called by controller, mapped to each engine. Activates random spines of increasing number until max cooperativity
    is reached.
    :param group_index: int
    :param num_spines: int
    :return: str
    """
    spine_group = groups_to_stim[group_index]
    path_type = spine_group['path_type']
    path_index = spine_group['path_index']
    sim.parameters['path_index'] = path_index
    sim.parameters['path_type'] = path_type
    sim.parameters['path_category'] = spine_group['path_category']
    spine = spine_group['spines'][0]
    sim.modify_rec(2, node=spine.parent.parent)
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
    spine = spine_group['spines'][spine_index]
    syn = spine.synapses[0]
    spike_times = h.Vector([equilibrate])
    syn.source.play(spike_times)
    sim.parameters['spine_index'] = spine.index
    sim.parameters['path_index'] = spine_group['path_index']
    sim.parameters['path_type'] = path_type
    sim.parameters['path_category'] = spine_group['path_category']
    sim.modify_rec(2, node=spine.parent.parent)
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
v_init = -67.
syn_types = ['AMPA_KIN', NMDA_type]

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

min_num_spines = 20
max_length = 30.
grouped_spines = {sec_type: [] for sec_type in ['basal', 'trunk', 'apical', 'tuft']}
for sec_type in ['basal', 'apical', 'tuft']:
    for branch in (branch for branch in cell._node_dict[sec_type] if branch.sec.L >= max_length and
                    len(branch.spines) >= min_num_spines):
        print branch.name, ', len: ', branch.sec.L, ', spines: ', len(branch.spines)
        branch_origin = cell.get_dendrite_origin(branch)
        if cell.is_terminal(branch):
            spine_list = get_clustered_spines(cell, branch_origin, branch.spines, min_num_spines, max_length, 1)
            if spine_list:
                loc1 = cell.get_distance_to_node(branch_origin, spine_list[0].parent)
                loc2 = cell.get_distance_to_node(branch_origin, spine_list[-1].parent)
                #print 'filtered spines (terminal): ', len(spine_list), ', start: ', loc1, ', end: ', loc2
                grouped_spines[sec_type].append({'path_category': 'terminal', 'spines': spine_list,
                                    'path_index': spine_list[0].index, 'path_type': spine_list[0].parent.parent.type})
            #else:
            #    print 'Branch end did not meet criterion'
        spine_list = get_clustered_spines(cell, branch_origin, branch.spines, min_num_spines, max_length, 0)
        if spine_list:
            loc1 = cell.get_distance_to_node(branch_origin, spine_list[0].parent)
            loc2 = cell.get_distance_to_node(branch_origin, spine_list[-1].parent)
            #print 'filtered spines (proximal): ', len(spine_list), ', start: ', loc1, ', end: ', loc2
            grouped_spines[sec_type].append({'path_category': 'proximal', 'spines': spine_list,
                                    'path_index': spine_list[0].index, 'path_type': spine_list[0].parent.parent.type})
        #else:
        #    print 'Branch start did not meet criterion'
trunk_paths = []
trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
if trunk_bifurcation:
    trunk_bifurcation = trunk_bifurcation[0]
    trunk_path = cell.tree.path_between_nodes(trunk_bifurcation, cell.tree.root)
    trunk_path.remove(cell.tree.root)
    trunk_path.reverse()
    trunk_paths.append(trunk_path)
    for child in [child for child in trunk_bifurcation.children if child.type == 'trunk']:
        for node in [node for node in child.children if node.type in ['tuft', 'trunk']]:
            if node.type == 'trunk':
                child = node
                break
            else:
                break
        trunk_path = cell.tree.path_between_nodes(child, trunk_bifurcation)
        trunk_path.remove(trunk_bifurcation)
        trunk_path.reverse()
        trunk_paths.append(trunk_path)
else:
    trunk_paths.append(list(cell.trunk))
sec_type = 'trunk'
for trunk_path in trunk_paths:
    spines_in_path = []
    for trunk in trunk_path:
        spines_in_path.extend(trunk.spines)
    print trunk.name, ', spines: ', len(spines_in_path)
    branch_origin = cell.tree.root
    spine_list = get_clustered_spines(cell, branch_origin, spines_in_path, min_num_spines, max_length, 1)
    if spine_list:
        loc1 = cell.get_distance_to_node(branch_origin, spine_list[0].parent)
        loc2 = cell.get_distance_to_node(branch_origin, spine_list[-1].parent)
        #print 'filtered spines (distal): ', len(spine_list), ', start: ', loc1, ', end: ', loc2
        grouped_spines[sec_type].append({'path_category': 'distal', 'spines': spine_list,
                                    'path_index': spine_list[0].index, 'path_type': spine_list[0].parent.parent.type})
    #else:
    #    print 'Branch end did not meet criterion'
    spine_list = get_clustered_spines(cell, branch_origin, spines_in_path, min_num_spines, max_length, 0)
    if spine_list:
        loc1 = cell.get_distance_to_node(branch_origin, spine_list[0].parent)
        loc2 = cell.get_distance_to_node(branch_origin, spine_list[-1].parent)
        #print 'filtered spines (proximal): ', len(spine_list), ', start: ', loc1, ', end: ', loc2
        grouped_spines[sec_type].append({'path_category': 'proximal', 'spines': spine_list,
                                    'path_index': spine_list[0].index, 'path_type': spine_list[0].parent.parent.type})
    #else:
    #    print 'Branch start did not meet criterion'


local_random = random.Random()
groups_to_stim = []
"""
groups_to_stim.extend(grouped_spines['trunk'])
local_random.seed(0)
max_num_branches = 2
for sec_type in ['basal', 'tuft']:
    terminal_groups = [spine_group for spine_group in grouped_spines[sec_type] if spine_group['path_category'] ==
                       'terminal']
    proximal_groups = [spine_group for spine_group in grouped_spines[sec_type] if spine_group['path_category'] ==
                       'proximal']
    groups_to_stim.extend([terminal_groups[i] for i in local_random.sample(range(len(terminal_groups)),
                                                                        min(max_num_branches, len(terminal_groups)))])
    groups_to_stim.extend([proximal_groups[i] for i in local_random.sample(range(len(proximal_groups)),
                                                                        min(max_num_branches, len(proximal_groups)))])
sec_type = 'apical'
for criteria in [lambda x: x <= 100., lambda x: x >= 175.]:
    terminal_groups = [spine_group for spine_group in grouped_spines[sec_type] if spine_group['path_category'] ==
                       'terminal' and criteria(cell.get_distance_to_node(cell.tree.root,
                                                                cell.get_dendrite_origin(spine_group['spines'][0])))]
    proximal_groups = [spine_group for spine_group in grouped_spines[sec_type] if spine_group['path_category'] ==
                       'proximal' and criteria(cell.get_distance_to_node(cell.tree.root,
                                                                cell.get_dendrite_origin(spine_group['spines'][0])))]
    groups_to_stim.extend([terminal_groups[i] for i in local_random.sample(range(len(terminal_groups)),
                                                                        min(max_num_branches, len(terminal_groups)))])
    groups_to_stim.extend([proximal_groups[i] for i in local_random.sample(range(len(proximal_groups)),
                                                                        min(max_num_branches, len(proximal_groups)))])
"""
for sec_type in ['basal', 'trunk', 'apical', 'tuft']:
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

sim.append_rec(cell, cell.tree.root, description='soma', loc=0.5)
sim.append_rec(cell, trunk, description='trunk', loc=0.)
sim.append_rec(cell, trunk, description='branch', loc=0.5)  # placeholder for local branch
sim.append_rec(cell, trunk, description='origin', loc=1.)  # placeholder for local dendrite origin
spine = groups_to_stim[0]['spines'][0]
sim.append_rec(cell, spine, object=spine.synapses[0].target(NMDA_type), param='_ref_g', description='NMDA_g')