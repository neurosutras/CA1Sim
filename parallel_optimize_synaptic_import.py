__author__ = 'Grace Ng'
from specify_cells3 import *
from plot_results import *
import random

"""
Iterates through every spine and activates AMPA_KIN synapses. Allows measurement of EPSP attenuation
and kinetics.

Modify a YAML file to include parameters necessary for this script. Then, set up ipcluster and run parallel_optimize_main.py
Current YAML filepath: data/optimize_synaptic_defaults.yaml
"""

context = Context()

def config_controller():
    set_constants()

def config_engine(update_params_funcs, param_names, rec_filepath, mech_file_path, neurotree_file_path, neurotree_index,
                  spines):
    """

    :param update_params_funcs: list of function references
    :param param_names: list of str
    :param rec_filepath: str
    :param mech_file_path: str
    :param neurotree_file_path: str
    :param neurotree_index: int
    :param spines: bool
    :return:
    """
    neurotree_dict = read_from_pkl(neurotree_file_path)[neurotree_index]
    param_indexes = {param_name: i for i, param_name in enumerate(param_names)}
    context.update(locals())
    set_constants()
    setup_cell()

def set_constants():
    """

    :return:
    """
    seed_offset = 7. * 2E6
    num_random_syn = 30
    num_clustered_syn = 20
    branch_list = ['branch1, branch2'] # number of branches for which clusters of synpases will be stimulated simultaneously
                                        # to calculate compound EPSP
    syn_ind_legend = {'random': range(num_random_syn)}
    for i, branch in enumerate(branch_list):
        syn_ind_legend[branch] = range(num_random_syn + i*num_clustered_syn, num_random_syn + (i+1)*num_clustered_syn)
    unitary_sim_dur = 100. # Duration (in ms) of a simulation where only one synapse is stimulated
    num_unitary_stims = 5 # Number of unitary stimulations (each separated by unitary_sim_dur) per sim
    equilibrate = 250.  # time to steady-state
    unitary_duration = equilibrate + num_unitary_stims * unitary_sim_dur
    compound_duration = 450.
    compound_stim_offset = 0.3
    dt = 0.02
    v_init = -77. #The optimize_branch_cooperativity_nmda_kin3_engine script says -67.
    NMDA_type = 'NMDA_KIN5'
    syn_types = ['AMPA_KIN', NMDA_type]
    local_random = random.Random()
    i_holding = 0.
    context.update(locals())

def setup_cell():
    """

    """
    cell = DG_GC(neurotree_dict=context.neurotree_dict, mech_file_path=context.mech_file_path, full_spines=context.spines)
    if context.spines is False:
        cell.correct_for_spines()
    cell.zero_na()
    context.local_random.seed(int(context.neurotree_index + context.seed_offset))

    # these synapses will not be used, but must be inserted for inheritance of synaptic parameters from trunk
    if context.spines:
        for branch in cell.apical:
            for spine in branch.spines:
                syn = Synapse(cell, spine, context.syn_types, stochastic=0)
    else:
        for branch in cell.apical:
            for syn_loc in branch.synapse_locs['excitatory']:
                syn = Synapse(cell, branch, context.syn_types, loc=syn_loc, stochastic=0)
    cell.init_synaptic_mechanisms()

    final_syn_list = {}
    random_syn_list = []
    for branch in cell.apical:
        if branch.synapses > 1:
            if branch.sec.L <= 10:
                candidates = [syn for syn in branch.synapses if 'AMPA_KIN' in syn._syn]
                random_syn_list.extend(context.local_random.sample(candidates, 1))
            else:
                num_syns = min(len(branch.synapses), int(branch.sec.L // 10.))
                for syn_ind in context.local_random.sample(range(len(branch.synapses)), num_syns):
                    if 'AMPA_KIN' in syn._syn:
                        random_syn_list.append(branch.synapses[syn_ind])
        elif branch.synapses:
            random_syn_list.append(branch.synapses[0])
        if context.spines:
            if len(branch.spines) > 1:
                if branch.sec.L <= 10.:
                    node = branch.spines[context.local_random.sample(range(0, len(branch.spines)), 1)[0]]
                    random_syn_list.extend([syn for syn in node.synapses if 'AMPA_KIN' in syn._syn])
                else:
                    num_syns = min(len(branch.spines), int(branch.sec.L // 10.))  # a random synapse every 10 um
                    for i in context.local_random.sample(range(0, len(branch.spines)), num_syns):
                        node = branch.spines[i]
                        random_syn_list.extend([syn for syn in node.synapses if 'AMPA_KIN' in syn._syn])
            elif branch.spines:
                node = branch.spines[0]
                random_syn_list.extend([syn for syn in node.synapses if 'AMPA_KIN' in syn._syn])
    final_syn_list['random'] = context.local_random.sample(random_syn_list, context.num_random_syn)

    # choose a distal apical oblique branch that has > 30 synapses within 30 um, choose synapses near the middle of the branch
    branch_list = [apical for apical in cell.apical if 150. < cell.get_distance_to_node(cell.tree.root, apical) < 250.
                   and apical.sec.L > 80.]
    branch = context.local_random.sample(branch_list, 1)[0]
    tested_branches = 0
    success_branches = 0
    while success_branches < 2:
        branch_syn_list = []
        if branch.synapses:
            branch_syn_list.extend([syn for syn in branch.synapses if ('AMPA_KIN' in syn._syn and
                                                                       30. <= syn.loc*branch.sec.L <= 60.)])
        if context.spines:
            for spine in branch.spines:
                branch_syn_list.extend([syn for syn in spine.synapses if ('AMPA_KIN' in syn._syn and
                                            30. <= (cell.get_distance_to_node(branch.parent, spine, loc=0.)) <= 60.)])
        if len(branch_syn_list) > context.num_clustered_syn:
            success_branches += 1
            label = 'branch' + str(success_branches)
            final_syn_list[label] = context.local_random.sample(branch_syn_list, context.num_clustered_syn)
        tested_branches += 1
        if tested_branches == len(branch_list) and success_branches < 2:
            raise Exception('Could not find a branch that satisfies the specified requirements.')
        branch = context.local_random.sample(branch_list, 1)[0]

    context.cell = cell
    context.syn_list = final_syn_list

    # get the thickest apical dendrite ~200 um from the soma
    candidate_branches = []
    candidate_diams = []
    candidate_locs = []
    for branch in cell.apical:
        if ((cell.get_distance_to_node(cell.tree.root, branch, 0.) >= 200.) &
                (cell.get_distance_to_node(cell.tree.root, branch, 1.) > 300.) & (not cell.is_terminal(branch))):
            candidate_branches.append(branch)
            for seg in branch.sec:
                loc = seg.x
                if cell.get_distance_to_node(cell.tree.root, branch, loc) > 250.:
                    candidate_diams.append(branch.sec(loc).diam)
                    candidate_locs.append(loc)
                    break
    index = candidate_diams.index(max(candidate_diams))
    dend = candidate_branches[index]

    rec_locs = {'soma': 0., 'thickest_dend': 0., 'local_branch': 0.}
    context.rec_locs = rec_locs
    rec_nodes = {'soma': cell.tree.root, 'thickest_dend': dend, 'local_branch': 0.}
    context.rec_nodes = rec_nodes

    equilibrate = context.equilibrate
    unitary_duration = context.unitary_duration
    dt = context.dt

    sim = QuickSim(unitary_duration, cvode=False, dt=dt, verbose=False)
    sim.parameters['equilibrate'] = equilibrate
    sim.parameters['duration'] = unitary_duration
    sim.parameters['spines'] = context.spines
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=unitary_duration) #Should the duration be this long for stim in offset_vm?
    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    context.sim = sim

def update_mech_dict(x, update_function, mech_file_path):
    update_function(x)
    context.cell.export_mech_dict(mech_file_path)

def get_unitary_EPSP_features(indiv, c, client_range, export=False):
    """
    Distribute simulations across available engines for testing EPSP attenuation.
    :param indiv: dict {'pop_id': pop_id, 'x': x arr, 'features': features dict}
    :param c: Client object
    :param client_range: list of ints
    :param export: False (for exporting voltage traces)
    :return: dict
    """
    dv = c[client_range]
    x = indiv['x']
    syn_groups = {}
    test_syns = []
    AP5_conditions = []
    group_types = []
    for type in context.syn_ind_legend:
        num_groups = len(context.syn_ind_legend[type]) // context.num_unitary_stims
        syn_groups[type] = [context.syn_ind_legend[type][i*context.num_unitary_stims:(i+1)*context.num_unitary_stims]
                             for i in range(num_groups)]
        extra_group = [context.syn_ind_legend[type][i] for i in range(num_groups * context.num_unitary_stims,
                                                                                len(context.syn_ind_legend[type]))]
        if extra_group:
            syn_groups[type].append(extra_group)
        test_syns.extend(syn_groups[type] * 2)
        AP5_conditions.extend(['AP5'] * len(syn_groups[type]) + ['con'] * len(syn_groups[type]))
        group_types.extend([type] * 2 * len(syn_groups[type]))
    result = dv.map_async(compute_EPSP_amp_features, [x] * len(test_syns), test_syns, AP5_conditions, group_types,
                          [export] * len(test_syns))
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result,
            'filter_features': filter_unitary_EPSP_features}

def filter_unitary_EPSP_features(get_result, old_features):
    """

    :param get_result: list of dict (each dict has the results from a particular simulation)
    :param old_features: dict (empty at this point in the optimization)
    :return: dict
    """
    random_AP5_results = {}
    random_con_results = {}
    clustered_results = {'branch1': {'AP5': {}, 'con': {}}, 'branch2': {'AP5': {}, 'con': {}}}
    for this_result in get_result:
        for syn_id, syn_result in this_result.iteritems():
            if syn_result['group_type'] == 'random':
                if syn_result['AP5_condition']:
                    random_AP5_results[syn_id] = syn_result['EPSP_amp']
                else:
                    random_con_results[syn_id] = syn_result['EPSP_amp']
            else:
                group_type = syn_result['group_type']
                AP5_condition = syn_result['AP5_condition']
                clustered_results[group_type][AP5_condition][syn_id] = syn_result
    avg_EPSP_AP5 = np.average(random_AP5_results.values())
    avg_EPSP_con = np.average(random_con_results.values())
    NMDA_contributions = []
    for syn_id in random_AP5_results.iterkeys():
        NMDA_contributions.append((random_con_results[syn_id] - random_AP5_results[syn_id])/random_con_results[syn_id])
    avg_NMDA_contr = np.average(NMDA_contributions)
    new_features = {'soma_EPSP_AP5': avg_EPSP_AP5, 'soma_EPSP_control': avg_EPSP_con,
                    'NMDA_contribution': avg_NMDA_contr, 'branch1_unitary': clustered_results['branch1'],
                    'branch2_unitary': clustered_results['branch2']}
    return new_features

def get_compound_EPSP_features(indiv, c, client_range, export=False):
    """
    Distribute simulations across available engines for testing spike stability.
    :param indiv: dict {'pop_id': pop_id, 'x': x arr, 'features': features dict}
    :param c: Client object
    :param client_range: list of ints
    :param export: False (for exporting voltage traces)
    :return: dict
    """
    dv = c[client_range]
    x = indiv['x']
    syn_groups = {}
    test_syns = []
    AP5_conditions = []
    group_types = []
    for type in context.branch_list:
        syn_groups[type] = [context.syn_ind_legend[type][0:i+1] for i in range(len(context.syn_ind_legend[type]))]
        test_syns.extend(syn_groups[type] * 2)
        AP5_conditions.extend(['AP5'] * len(syn_groups[type]) + ['con'] * len(syn_groups[type]))
        group_types.extend([type] * 2 * len(syn_groups[type]))
    result = dv.map_async(compute_branch_cooperativity_features, [x] * len(test_syns), test_syns, AP5_conditions, group_types,
                          [export] * len(test_syns))
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result,
            'filter_features': filter_unitary_EPSP_features}

def filter_compound_EPSP_features(get_result, old_features):
    """

    :param get_result: list of dict (each dict has the results from a particular simulation)
    :param old_features: dict
    :return: dict
    """
    expected_compound_EPSP = {}
    actual_compound_EPSP = {}
    sim_id_list = {}
    for group_type in context.branch_list:
        expected_compound_EPSP[group_type] = compute_expected_compound_EPSP(old_features[group_type + '_unitary'])
        actual_compound_EPSP[group_type] = {'AP5': [], 'con': []}
        sim_id_list[group_type] = {'AP5': [], 'con': []}
    for this_result in get_result:
        for num_sim, sim_id in this_result.iterkeys():
            sim_result = this_result[sim_id]
            group_type = sim_result['group_type']
            AP5_condition = sim_result['AP5_condition']
            sim_id_list[group_type][AP5_condition].append(sim_id)
            actual_compound_EPSP[group_type][AP5_condition].append(np.max(sim_result['rec']['soma']))
                # This forms an array of the compound EPSPs, where each element is the EPSP from a certain number of synapses
    for group_type in sim_id_list:
        for AP5_condition in sim_id_list[group_type]:
            actual_compound_EPSP[group_type][AP5_condition].sort(key=dict(zip(actual_compound_EPSP[group_type][AP5_condition],
                                                                              sim_id_list[group_type][AP5_condition])).get)
    peak_supralinearities = {'AP5': [], 'con': []}
    min_supralinearities = {'AP5': [], 'con': []}
    for group_type in context.branch_list:
        for AP5_condition in ['AP5', 'con']:
            actual = np.array(actual_compound_EPSP[group_type][AP5_condition])
            expected = np.array(expected_compound_EPSP[group_type][AP5_condition])
            supralinearity = (actual - expected) / expected * 100.
            peak_supralinearity = np.max(supralinearity)
            if peak_supralinearity < 0.:  # there is no gradient if integration is always sublinear
                peak_supralinearity = np.min(supralinearity)  # exaggerate error for sublinear integration
                min_supralinearity = np.min(supralinearity)
            peak_supralinearities[AP5_condition].append(peak_supralinearity)
            min_supralinearities[AP5_condition].append(min_supralinearity)
    new_features = {'peak_supralinearity_AP5': np.average(peak_supralinearities['AP5']),
                    'peak_supralinearity_con': np.average(peak_supralinearities['con']),
                    'min_supralinearity_AP5': np.average(min_supralinearities['AP5']),
                    'min_supralinearity_con': np.average(min_supralinearities['con'])}
    return new_features

def get_objectives(features, objective_names, target_val, target_range):
    """

    :param features: dict
    :return: dict
    """
    objectives = {}
    for objective_name in objective_names:
        objectives[objective_name] = ((target_val[objective_name] - features[objective_name]) /
                                                  target_range[objective_name]) ** 2.
    return features, objectives

def compute_EPSP_amp_features(x, test_syns, AP5_condition, group_type, export=False, plot=False):
    """

    :param x: arr
    :param test_syns: list of int
    :param AP5_condition: str
    :param group_type: str
    :param export: bool
    :param plot: bool
    :return: dict
    """
    start_time = time.time()
    context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    for update_func in context.update_params_funcs:
        context.cell = update_func(context.cell, x, context.param_indexes)
    # sim.cvode_state = True
    soma_vm = offset_vm('soma', context.v_init)
    synapses = [context.syn_list[syn_index] for syn_index in test_syns]
    stim_times = [h.Vector([context.equilibrate + i*context.unitary_sim_dur]) for i in range(context.num_unitary_stims)]
    for i, syn in enumerate(synapses):
        syn.source.play(stim_times[i])
        if i == 0:
            if context.spines:
                spine = syn.node
                branch = spine.parent.parent
                context.sim.parameters['input_loc'] = branch.type
            else:
                branch = syn.node
                context.sim.parameters['input_loc'] = branch.type
        if AP5_condition == 'AP5':
            syn.target(context.NMDA_type).gmax = 0.
    description = 'Unitary Stim %s: Synapses %s' %(AP5_condition, str(test_syns[0]) + '-' + str(test_syns[-1]))
    context.sim.parameters['description'] = description
    for param_name in context.param_names:
        context.sim.parameters[param_name] = x[context.param_indexes[param_name]]
    context.sim.run()
    t = np.array(context.sim.tvec[:])
    duration = context.unitary_duration
    dt = context.dt
    equilibrate = context.equilibrate
    interp_t = np.arange(0, duration, dt)
    left, right = time2index(interp_t, equilibrate-3.0, equilibrate-1.0)
    result = {}
    for rec_ind, rec in enumerate(context.sim.rec_list):
        vm = np.array(rec['vec'][:])
        interp_vm = np.interp(interp_t, t, vm)
        for i, syn in enumerate(synapses):
            start, end = time2index(interp_t, equilibrate + i*context.unitary_sim_dur, equilibrate + (i+1)*context.unitary_sim_dur)
            baseline = np.average(interp_vm[left:right])
            corrected_vm = interp_vm - baseline
            peak = np.max(corrected_vm[start:end])
            peak_index = np.where(corrected_vm == peak)[0][0]
            zero_index = np.where(corrected_vm[peak_index:] <= 0.)[0]
            if np.any(zero_index):
                corrected_vm[peak_index+zero_index[0]:] = 0.
            distance = context.cell.get_distance_to_node(context.cell.tree.root, syn.node, syn.loc)
            syn_id = test_syns[i]
            if syn_id not in result.keys():
                result[syn_id] = {'AP5_condition': AP5_condition, 'group_type': group_type, 'distance': distance,
                                  'rec': {}}
            if rec['description'] is 'soma':
                result[syn_id]['EPSP_amp'] = peak
            result[syn_id]['rec'][rec['description']] = corrected_vm[start:end]
            if rec_ind == len(context.sim.rec_list) - 1:
                corrected_t = interp_t[start:end]
                result[syn_id]['tvec'] = corrected_t - corrected_t[0] #sets each tvec to start at 0
                syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while
                                            # keeping the VecStim source object in existence so it can be activated again
                print 'Process:', os.getpid(), 'optimized synapse:', test_syns[i], 'on Node:', syn.node.name, \
                        'distance: %.2f in %.3f s' % (distance, time.time() - start_time)
    if plot:
        context.sim.plot()
    if export:
        export_sim_results()
    return result

def compute_branch_cooperativity_features(x, test_syns, AP5_condition, group_type, export=False, plot=False):
    """
    This simulation tests a value for gmax_NMDA_KIN. Results are later compared to branch cooperativity data from
    Harnett et al., 2012 (performed in TTX). One apical oblique ~150 um from the soma is chosen, up to 50 spines within a
    30 um long stretch of branch are stimulated. Actual and Expected EPSP amplitudes are compared, with a target peak
    nonlinearity of 44 +/- 6%.

    :param x: arr
    :param test_syns: list of int
    :param AP5_condition: str
    :param group_type: str
    :param export: bool
    :param plot: bool
    :return:
    """
    start_time = time.time()
    context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    for update_func in context.update_params_funcs:
        context.cell = update_func(context.cell, x, context.param_indexes)
    # sim.cvode_state = True
    synapses = [context.syn_list[syn_index] for syn_index in test_syns]
    stim_times = [h.Vector([context.equilibrate + i * context.compound_stim_offset]) for i in range(len(test_syns))]
    for i, syn in enumerate(synapses):
        syn.source.play(stim_times[i])
        if i == 0:
            if context.spines:
                spine = syn.node
                branch = spine.parent.parent
                context.sim.modify_rec(2, branch, 0.)
                context.sim.parameters['input_loc'] = branch.type
            else:
                branch = syn.node
                context.sim.modify_rec(2, branch, 0.)
                context.sim.parameters['input_loc'] = branch.type
        if AP5_condition == 'AP5':
            syn.target(context.NMDA_type).gmax = 0.
    description = 'Compound Stim %s: Synapses %s' %(AP5_condition, str(test_syns[0]) + '-' + str(test_syns[-1]))
    context.sim.parameters['description'] = description
    context.sim.run(context.v_init)

    t = np.array(context.sim.tvec[:])
    duration = context.compound_duration
    dt = context.dt
    equilibrate = context.equilibrate
    interp_t = np.arange(0, duration, dt)
    left, right = time2index(interp_t, equilibrate - 3.0, equilibrate - 1.0)
    start, end = time2index(interp_t, equilibrate, duration)
    sim_id = test_syns[-1]
    corrected_t = interp_t[start:end]
    result = {sim_id: {'AP5_condition': AP5_condition, 'group_type': group_type, 'tvec': corrected_t - corrected_t[0],
                       'rec': {}}}
    for rec_ind, rec in enumerate(context.sim.rec_list):
        vm = np.array(rec['vec'][:])
        interp_vm = np.interp(interp_t, t, vm)
        baseline = np.average(interp_vm[left:right])
        corrected_vm = interp_vm - baseline
        peak = np.max(corrected_vm[start:end])
        peak_index = np.where(corrected_vm == peak)[0][0]
        zero_index = np.where(corrected_vm[peak_index:] <= 0.)[0]
        if np.any(zero_index):
            corrected_vm[peak_index + zero_index[0]:] = 0.
        result[sim_id]['rec'][rec['description']] = corrected_vm[start:end]
        for syn in synapses:
            syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while
                                        # keeping the VecStim source object in existence so it can be activated again
    print 'Process: %i stimulated %i synapses in %i s with gmax: %.3E' % (os.getpid(), len(test_syns), time.time() - start_time,
                                                                          x[context.param_indexes['NMDA.gmax']])
    if plot:
        context.sim.plot()
    if export:
        export_sim_results()
    return result

def compute_expected_compound_EPSP(unitary_branch_results):
    """

    :param unitary_branch_features: dict {'AP5': {syn_id: {'rec': {loc: trace}}
    :return: dict {'AP5': {sim_id: {loc: summed_traces, ...}, sim_id: ...}, 'con': {sim_id: {loc: summed_traces, ..}}}
    """
    expected_compound_traces = {}
    expected_compound_EPSP = {}
    for AP5_condition in ['AP5', 'con']:
        expected_compound_traces[AP5_condition] = {}
        expected_compound_EPSP[AP5_condition] = []
        summed_traces = {}
        total_syns = len(unitary_branch_results[AP5_condition])
        for num_syn, syn_id in enumerate(unitary_branch_results[AP5_condition].keys()):
            syn_result = unitary_branch_results[AP5_condition][syn_id]
            for loc in unitary_branch_results[syn_id]['rec'].iterkeys():
                if loc not in summed_traces:
                    offset_inds = math.ceil(context.compound_stim_offset / context.unitary_sim_dur *
                                            len(syn_result['rec'][loc])) # Number of indices in an interpolated time vector
                                                                        # corresponding to the offset time between stims
                    unitary_sim_inds = len(syn_result['rec'][loc]) # Number of indices in an interpolated time vector
                                                                    # corresponding to length of a unitary simulation
                    summed_traces[loc] = np.zeros(unitary_sim_inds + offset_inds * (total_syns - 1))
                start_ind = num_syn * offset_inds
                summed_traces[loc][start_ind:start_ind+unitary_sim_inds] += syn_result['rec'][loc]
            expected_compound_traces[AP5_condition][syn_id] = summed_traces # This stores the summed traces of all the synapses
                                                                            # up to and including the synapse with this syn_id
            expected_compound_EPSP[AP5_condition].append(np.max(summed_traces['soma']))
    return expected_compound_EPSP

def offset_vm(description, vm_target=None):
    """

    :param description: str
    :param vm_target: float
    """
    if vm_target is None:
        vm_target = context.v_init
    context.sim.modify_stim(0, amp=0.)
    node = context.rec_nodes[description]
    loc = context.rec_locs[description]
    rec_dict = context.sim.get_rec(description)
    context.sim.modify_stim(0, node=node, loc=loc, amp=0.)
    rec = rec_dict['vec']
    offset = True

    equilibrate = context.equilibrate
    dt = context.dt
    duration = context.unitary_duration

    context.sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    context.sim.modify_stim(0, amp=context.i_holding[description])
    context.sim.run(vm_target)
    vm = np.interp(t, context.sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < vm_target - 0.5:
        context.i_holding[description] += 0.01
        while offset:
            if context.sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (context.i_holding[description], description)
            context.sim.modify_stim(0, amp=context.i_holding[description])
            context.sim.run(vm_target)
            vm = np.interp(t, context.sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest < vm_target - 0.5:
                context.i_holding[description] += 0.01
            else:
                offset = False
    elif v_rest > vm_target + 0.5:
        context.i_holding[description] -= 0.01
        while offset:
            if context.sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (context.i_holding[description], description)
            context.sim.modify_stim(0, amp=context.i_holding[description])
            context.sim.run(vm_target)
            vm = np.interp(t, context.sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                context.i_holding[description] -= 0.01
            else:
                offset = False
    context.sim.tstop = duration
    return v_rest

def update_AMPA_NMDA(cell, x, param_indexes):
    """

    :param x: array. Default: value=1.711e-03, slope=6.652e-05, tau=7.908e+01
    :return:
    """
    if context.spines is False:
        cell.correct_for_spines()
    cell.modify_mech_param('apical', 'synapse', 'gmax', value=x[param_indexes['AMPA.g0']], origin='soma', syn_type='AMPA_KIN',
                                   slope=x[param_indexes['AMPA.slope']], tau=x[param_indexes['AMPA.tau']])
    cell.modify_mech_param('apical', 'synapse', 'gmax', value=x[param_indexes['AMPA.g0']], syn_type='AMPA_KIN',
                           custom={'method': 'custom_inheritance_by_nonterm_branchord', 'branch_cutoff': 4}, replace=False)
    cell.modify_mech_param('apical', 'synapse', 'gmax', value=x[param_indexes['NMDA.gmax']], syn_type=context.NMDA_type)
    cell.modify_mech_param('apical', 'synapse', 'gamma', value=x[param_indexes['NMDA.gamma']], syn_type=context.NMDA_type)
    cell.modify_mech_param('apical', 'synapse', 'Kd', value=x[param_indexes['NMDA.Kd']], syn_type=context.NMDA_type)
    cell.modify_mech_param('apical', 'synapse', 'kin_scale', value=x[param_indexes['NMDA.Kd']], syn_type=context.NMDA_type)
    cell.init_synaptic_mechanisms()
    return cell

def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(context.rec_filepath, 'a') as f:
        context.sim.export_to_file(f)