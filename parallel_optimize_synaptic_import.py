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

def setup_module_from_file(param_file_path='data/optimize_synaptic_defaults.yaml', rec_file_path=None,
                           export_file_path=None, verbose=False):
    """

    :param param_file_path: str (path to a yaml file)
    :return:
    """
    params_dict = read_from_yaml(param_file_path)
    param_gen = params_dict['param_gen']
    param_names = params_dict['param_names']
    default_params = params_dict['default_params']
    x0 = params_dict['x0']
    bounds = [params_dict['bounds'][key] for key in param_names]
    feature_names = params_dict['feature_names']
    objective_names = params_dict['objective_names']
    target_val = params_dict['target_val']
    target_range = params_dict['target_range']
    optimization_title = params_dict['optimization_title']
    kwargs = params_dict['kwargs']  # Extra arguments
    kwargs['verbose'] = verbose
    update_params = params_dict['update_params']
    update_params_funcs = []
    for update_params_func_name in update_params:
        func = globals().get(update_params_func_name)
        if not callable (func):
            raise Exception('Multi-Objective Optimization: update_params: %s is not a callable function.'
                            % (update_params_func_name))
    if rec_file_path is None:
        rec_file_path = 'data/sim_output' + datetime.datetime.today().strftime('%m%d%Y%H%M') + \
                   '_pid' + str(os.getpid()) + '.hdf5'
    if export_file_path is None:
        export_file_path = 'data/%s_%s_%s_optimization_exported_traces.hdf5' % \
                           (datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)
    context.update(locals())
    config_engine(update_params_funcs, param_names, default_params, rec_file_path, export_file_path, **kwargs)

def config_controller(export_file_path):
    context.update(locals())
    set_constants()

def config_engine(update_params_funcs, param_names, default_params, rec_file_path, export_file_path, mech_file_path,
                  neurotree_file_path, neurotree_index, spines, verbose=False):
    """

    :param update_params_funcs: list of function references
    :param param_names: list of str
    :param default_params: dict
    :param rec_file_path: str
    :param export_file_path: str
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
    setup_cell(verbose)

def set_constants():
    """

    :return:
    """
    seed_offset = 7. * 2E6
    num_random_syn = 3 #30
    num_clustered_syn = 2 #20
    branch_list = ['branch1', 'branch2'] # number of branches for which clusters of synpases will be stimulated simultaneously
                                        # to calculate compound EPSP
    group_type_list = ['random'] + branch_list
    syn_ind_legend = {'random': range(num_random_syn)}
    for i, branch in enumerate(branch_list):
        syn_ind_legend[branch] = range(num_random_syn + i*num_clustered_syn, num_random_syn + (i+1)*num_clustered_syn)
    AP5_cond_list = ['AP5', 'con']
    unitary_isi = 125. # Inter-stimulus interval for unitary sims
    num_unitary_stims = 3 #5 # Number of unitary stims (each separated by isi) in a simulation
    equilibrate = 250.  # time to steady-state
    unitary_duration = equilibrate + num_unitary_stims * unitary_isi
    compound_duration = 450.
    compound_isi = 0.3
    stim_dur = 150.
    th_dvdt = 10.
    dt = 0.02
    v_init = -77. #The optimize_branch_cooperativity_nmda_kin3_engine script says -67.
    v_active = -77.
    NMDA_type = 'NMDA_KIN5'
    syn_types = ['AMPA_KIN', NMDA_type]
    local_random = random.Random()
    i_holding = {'soma': 0.}
    i_th = {'soma': 0.1}
    context.update(locals())

def setup_cell(verbose=False):
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

    final_syn_dict = {}

    # Choose random synapses, which will be used to calculate the average unitary EPSP
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
    final_syn_dict['random'] = context.local_random.sample(random_syn_list, context.num_random_syn)

    # Choose clustered synapses from two distal apical oblique branches to find compound EPSP. Each branch must have
    # > 30 synapses within 30 um, choose synapses near the middle of the branch.
    branch_list = [apical for apical in cell.apical if 150. < cell.get_distance_to_node(cell.tree.root, apical) < 250.
                   and apical.sec.L > 80.]
    branch = context.local_random.sample(branch_list, 1)[0]
    tested_branches = 0
    success_branches = 0
    while success_branches < len(context.branch_list):
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
            label = context.branch_list[success_branches-1]
            final_syn_dict[label] = context.local_random.sample(branch_syn_list, context.num_clustered_syn)
        tested_branches += 1
        if tested_branches == len(branch_list) and success_branches < len(context.branch_list):
            raise Exception('Could not find a branch that satisfies the specified requirements.')
        branch = context.local_random.sample(branch_list, 1)[0]

    context.cell = cell
    syn_list = []
    for group_type in context.group_type_list:
        syn_list.extend(final_syn_dict[group_type])
    context.syn_list = syn_list

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

    rec_locs = {'soma': 0., 'dend': 0., 'local_branch': 0.}
    context.rec_locs = rec_locs
    rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'local_branch': dend}
    context.rec_nodes = rec_nodes

    equilibrate = context.equilibrate
    unitary_duration = context.unitary_duration
    dt = context.dt
    stim_dur = context.stim_dur

    sim = QuickSim(unitary_duration, cvode=False, dt=dt, verbose=verbose)
    sim.parameters['equilibrate'] = equilibrate
    sim.parameters['duration'] = unitary_duration
    sim.parameters['spines'] = context.spines
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=context.compound_duration)
    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    context.sim = sim

    context.spike_output_vec = h.Vector()
    cell.spike_detector.record(context.spike_output_vec)

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
        for AP5_condition in context.AP5_cond_list:
            test_syns.extend(syn_groups[type])
            AP5_conditions.extend([AP5_condition] * len(syn_groups[type]))
            group_types.extend([type] * len(syn_groups[type]))
    result = dv.map_async(compute_EPSP_amp_features, [x] * len(test_syns), test_syns, AP5_conditions, group_types,
                          [export] * len(test_syns))
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result,
            'filter_features': filter_unitary_EPSP_features}

def filter_unitary_EPSP_features(get_result, old_features, export):
    """

    :param get_result: list of dict (each dict has the results from a particular simulation)
    :param old_features: dict (empty at this point in the optimization)
    :param export: bool
    :return: dict
    """
    random_AP5_results = {}
    random_con_results = {}
    clustered_results = {}
    for branch in context.branch_list:
        clustered_results[branch] = {}
        for AP5_cond in context.AP5_cond_list:
            clustered_results[branch][AP5_cond] = {}
    for this_result in get_result:
        for syn_id, syn_result in this_result.iteritems():
            if syn_result['group_type'] == 'random':
                if syn_result['AP5_condition'] == 'AP5':
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
                    'NMDA_contribution': avg_NMDA_contr}
    for branch in context.branch_list:
        new_features[branch+'_unitary'] = clustered_results[branch]
    if export:
        new_export_file_path = context.export_file_path.replace('.hdf5', '_processed.hdf5')
        with h5py.File(new_export_file_path, 'a') as f:
            f.create_group('actual_compound_traces')
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
        for AP5_condition in context.AP5_cond_list:
            test_syns.extend(syn_groups[type])
            AP5_conditions.extend([AP5_condition] * len(syn_groups[type]))
            group_types.extend([type] * len(syn_groups[type]))
    result = dv.map_async(compute_branch_cooperativity_features, [x] * len(test_syns), test_syns, AP5_conditions, group_types,
                          [export] * len(test_syns))
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result,
            'filter_features': filter_compound_EPSP_features}

def filter_compound_EPSP_features(get_result, old_features, export):
    """

    :param get_result: list of dict (each dict has the results from a particular simulation)
    :param old_features: dict
    :param export: bool
    :return: dict
    """
    expected_compound_traces = {}
    expected_compound_EPSP = {}
    actual_compound_traces = {}
    actual_compound_EPSP = {}
    sim_id_list = {}
    for group_type in context.branch_list:
        expected_compound_traces[group_type], expected_compound_EPSP[group_type] = \
            get_expected_compound_features(old_features[group_type + '_unitary'])
        actual_compound_traces[group_type] = {}
        actual_compound_EPSP[group_type] = {}
        sim_id_list[group_type] = {}
        for AP5_cond in context.AP5_cond_list:
            actual_compound_traces[group_type][AP5_cond] = {}
            actual_compound_EPSP[group_type][AP5_cond] = []
            sim_id_list[group_type][AP5_cond] = []
    for this_result in get_result:
        for num_sim, sim_id in enumerate(this_result.keys()):
            sim_result = this_result[sim_id]
            group_type = sim_result['group_type']
            AP5_condition = sim_result['AP5_condition']
            sim_id_list[group_type][AP5_condition].append(sim_id)
            actual_compound_traces[group_type][AP5_condition][sim_id] = {'rec': {}}
            for loc in sim_result['rec']:
                actual_compound_traces[group_type][AP5_condition][sim_id]['rec'][loc] = sim_result['rec'][loc]
            actual_compound_traces[group_type][AP5_condition][sim_id]['tvec'] = sim_result['tvec']
            actual_compound_EPSP[group_type][AP5_condition].append(np.max(sim_result['rec']['soma']))
                # This forms an array of the compound EPSPs, where each element is the EPSP from a certain number of synapses
    for group_type in sim_id_list:
        for AP5_condition in sim_id_list[group_type]:
            actual_compound_EPSP[group_type][AP5_condition].sort(key=dict(zip(actual_compound_EPSP[group_type][AP5_condition],
                                                                              sim_id_list[group_type][AP5_condition])).get)
    peak_supralinearities = {}
    min_supralinearities = {}
    for AP5_cond in context.AP5_cond_list:
        peak_supralinearities[AP5_cond] = []
        min_supralinearities[AP5_cond] = []
    for group_type in context.branch_list:
        for AP5_condition in context.AP5_cond_list:
            actual = np.array(actual_compound_EPSP[group_type][AP5_condition])
            expected = np.array(expected_compound_EPSP[group_type][AP5_condition])
            supralinearity = (actual - expected) / expected * 100.
            peak_supralinearity = np.max(supralinearity)
            if peak_supralinearity < 0.:  # there is no gradient if integration is always sublinear
                peak_supralinearity = np.min(supralinearity)  # exaggerate error for sublinear integration
                min_supralinearity = np.min(supralinearity)
            else:
                peak_index = np.where(supralinearity == peak_supralinearity)[0][0]
                if peak_index == 0:
                    min_supralinearity = supralinearity[0]
                else:
                    min_supralinearity = np.min(supralinearity[:peak_index])
            peak_supralinearities[AP5_condition].append(peak_supralinearity)
            min_supralinearities[AP5_condition].append(min_supralinearity)
    new_features = {'peak_supralinearity_AP5': np.average(peak_supralinearities['AP5']),
                    'peak_supralinearity_con': np.average(peak_supralinearities['con']),
                    'min_supralinearity_AP5': np.average(min_supralinearities['AP5']),
                    'min_supralinearity_con': np.average(min_supralinearities['con'])}
    if export:
        new_export_file_path = context.export_file_path.replace('.hdf5', '_processed.hdf5')
        with h5py.File(new_export_file_path, 'a') as f:
            f.create_group('actual_compound_traces')
            f.create_group('expected_compound_traces')
            for i, trace_dict in enumerate([actual_compound_traces, expected_compound_traces]):
                if i == 0:
                    trace_group = f['actual_compound_traces']
                else:
                    trace_group = f['expected_compound_traces']
                for group_type in trace_dict:
                    group = trace_group.create_group(group_type)
                    for AP5_condition in trace_dict[group_type]:
                        AP5_group = group.create_group(AP5_condition)
                        for sim_id in trace_dict[group_type][AP5_condition]:
                            sim_group = AP5_group.create_group(str(sim_id))
                            sim_group.create_group('rec')
                            sim_group.create_dataset('tvec', compression='gzip', compression_opts=9,
                                                               data=trace_dict[group_type][AP5_condition][sim_id]['tvec'])
                            for loc in trace_dict[group_type][AP5_condition][sim_id]['rec']:
                                sim_group['rec'].create_dataset(loc, compression='gzip', compression_opts=9,
                                                            data=trace_dict[group_type][AP5_condition][sim_id]['rec'][loc])
    return new_features

def get_stability_features(indiv, c, client_range, export=False):
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
    result = dv.map_async(compute_stability_features, [x], [export])
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result}

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
    test_syns = np.array(test_syns)
    start_time = time.time()
    context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    for update_func in context.update_params_funcs:
        context.cell = update_func(context.cell, x, context.param_indexes, context.default_params)
    # sim.cvode_state = True
    soma_vm = offset_vm('soma', context.v_init)
    synapses = [context.syn_list[syn_index] for syn_index in test_syns]
    for i, syn in enumerate(synapses):
        syn.source.play(h.Vector([context.equilibrate + i*context.unitary_isi]))
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
    description = 'Unitary Stim %s: Synapses %s' % (AP5_condition, str(test_syns[0]) + '-' + str(test_syns[-1]))
    context.sim.parameters['description'] = description
    for param_name in context.param_names:
        context.sim.parameters[param_name] = x[context.param_indexes[param_name]]
    duration = context.unitary_duration
    context.sim.tstop = duration
    context.sim.run()
    t = np.array(context.sim.tvec)
    dt = context.dt
    equilibrate = context.equilibrate
    interp_t = np.arange(0, duration, dt)
    left, right = time2index(interp_t, equilibrate-3.0, equilibrate-1.0)
    #baseline should change for each interval
    result = {}
    for rec_ind, rec in enumerate(context.sim.rec_list):
        vm = np.array(rec['vec'])
        interp_vm = np.interp(interp_t, t, vm)
        for i, syn in enumerate(synapses):
            start, end = time2index(interp_t, equilibrate + i*context.unitary_isi, equilibrate + (i+1)*context.unitary_isi)
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
                print 'Process:', os.getpid(), 'optimized %s synapse' %(group_type), test_syns[i], 'on Node', \
                    syn.node.name, 'with distance: %.2f' % (distance)
    print 'Process:', os.getpid(), 'completed in %.3f s' % (time.time() - start_time)
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
    test_syns = np.array(test_syns)
    context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    for update_func in context.update_params_funcs:
        context.cell = update_func(context.cell, x, context.param_indexes, context.default_params)
    # sim.cvode_state = True
    synapses = [context.syn_list[syn_index] for syn_index in test_syns]
    for i, syn in enumerate(synapses):
        syn.source.play(h.Vector([context.equilibrate + i * context.compound_isi]))
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
    duration = context.compound_duration
    context.sim.tstop = duration
    context.sim.run(context.v_init)

    t = np.array(context.sim.tvec)
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
        vm = np.array(rec['vec'])
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

def compute_stability_features(x, export=False, plot=False):
    """
    :param local_x: array [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor, soma.sh_nas/x,
                    axon.gkbar factor, dend.gkabar factor]
    :param plot: bool
    :return: float
    """
    start_time = time.time()
    context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    for update_func in context.update_params_funcs:
        context.cell = update_func(context.cell, x, context.param_indexes, context.default_params)
    # sim.cvode_state = True

    v_active = context.v_active
    equilibrate = context.equilibrate
    dt = context.dt
    i_th = context.i_th

    soma_vm = offset_vm('soma', v_active)
    stim_dur = 150.
    context.sim.modify_stim(0, node=context.cell.tree.root, loc=0., dur=stim_dur)
    duration = equilibrate + stim_dur
    context.sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    d_amp = 0.01
    amp = max(0., i_th['soma'] - 0.02)
    while not spike:
        context.sim.modify_stim(0, amp=amp)
        context.sim.run(v_active)
        vm = np.interp(t, context.sim.tvec, context.sim.get_rec('soma')['vec'])
        if np.any(vm[:int(equilibrate/dt)] > -30.):
            print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
            return None
        if np.any(vm[int(equilibrate/dt):int((equilibrate+50.)/dt)] > -30.):
            spike = True
        elif amp >= 0.4:
            print 'Process %i: Aborting - rheobase outside target range' % (os.getpid())
            return None
        else:
            amp += d_amp
            if context.sim.verbose:
                print 'increasing amp to %.3f' % amp
    context.sim.parameters['amp'] = amp
    context.sim.parameters['description'] = 'spike shape'
    i_th['soma'] = amp
    spike_times = context.cell.spike_detector.get_recordvec().to_python()
    peak, threshold, ADP, AHP = get_spike_shape(vm, spike_times)
    dend_vm = np.interp(t, context.sim.tvec, context.sim.get_rec('dend')['vec'])
    th_x = np.where(vm[int(equilibrate / dt):] >= threshold)[0][0] + int(equilibrate / dt)
    if len(spike_times) > 1:
        end = min(th_x + int(10. / dt), int((spike_times[1] - 5.)/dt))
    else:
        end = th_x + int(10. / dt)
    dend_peak = np.max(dend_vm[th_x:end])
    dend_pre = np.mean(dend_vm[int((equilibrate - 3.) / dt):int((equilibrate - 1.) / dt)])
    result = {'dend_amp': (dend_peak - dend_pre) / (peak - soma_vm)}
    print 'Process %i took %.1f s to find spike rheobase at amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    if plot:
        context.sim.plot()
    if export:
        export_sim_results()
    return result

def get_expected_compound_features(unitary_branch_results):
    """

    :param unitary_branch_features: dict {'AP5': {syn_id: {'rec': {loc: trace}}
    :return: dict {'AP5': {sim_id: {loc: summed_traces, ...}, sim_id: ...}, 'con': {sim_id: {loc: summed_traces, ..}}}
    """
    expected_compound_traces = {}
    expected_compound_EPSP = {}
    baseline_len = 0 # int(10. / context.dt)
    unitary_len = int(context.unitary_isi / context.dt)
    compound_isi_len = int(context.compound_isi / context.dt)
    actual_trace_len = int((context.compound_duration - context.equilibrate) / context.dt) + baseline_len
    for AP5_condition in context.AP5_cond_list:
        expected_compound_traces[AP5_condition] = {}
        expected_compound_EPSP[AP5_condition] = []
        summed_traces = {}
        total_syns = len(unitary_branch_results[AP5_condition])
        start_ind = baseline_len
        for num_syn, syn_id in enumerate(unitary_branch_results[AP5_condition]):
            syn_result = unitary_branch_results[AP5_condition][syn_id]
            for loc in unitary_branch_results[AP5_condition][syn_id]['rec'].iterkeys():
                if loc not in summed_traces:
                    summed_traces[loc] = np.zeros(actual_trace_len)

                summed_traces[loc][start_ind:start_ind+unitary_len] += syn_result['rec'][loc][baseline_len:]
            expected_compound_traces[AP5_condition][syn_id] = {}
            expected_compound_traces[AP5_condition][syn_id]['rec'] = copy.deepcopy(summed_traces) # This stores the summed traces of
                                                                                    # all the synapses up to and including
                                                                                    # the synapse with this syn_id
            expected_compound_EPSP[AP5_condition].append(np.max(summed_traces['soma']))
            start_ind += compound_isi_len
    return expected_compound_traces, expected_compound_EPSP

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
    context.sim.modify_stim(1, node=node, loc=loc, amp=0.)
    rec = rec_dict['vec']
    offset = True

    equilibrate = context.equilibrate
    dt = context.dt
    duration = context.compound_duration

    context.sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    context.sim.modify_stim(1, amp=context.i_holding[description])
    context.sim.run(vm_target)
    vm = np.interp(t, context.sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < vm_target - 0.5:
        context.i_holding[description] += 0.01
        while offset:
            if context.sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (context.i_holding[description], description)
            context.sim.modify_stim(1, amp=context.i_holding[description])
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
            context.sim.modify_stim(1, amp=context.i_holding[description])
            context.sim.run(vm_target)
            vm = np.interp(t, context.sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                context.i_holding[description] -= 0.01
            else:
                offset = False
    context.sim.tstop = duration
    return v_rest

def get_spike_shape(vm, spike_times):
    """

    :param vm: array
    :return: tuple of float: (v_peak, th_v, ADP, AHP)
    """
    equilibrate = context.equilibrate
    dt = context.dt
    th_dvdt = context.th_dvdt

    start = int((equilibrate+1.)/dt)
    vm = vm[start:]
    dvdt = np.gradient(vm, dt)
    th_x = np.where(dvdt > th_dvdt)[0]
    if th_x.any():
        th_x = th_x[0] - int(1.6/dt)
    else:
        th_x = np.where(vm > -30.)[0][0] - int(2./dt)
    th_v = vm[th_x]
    v_before = np.mean(vm[th_x-int(0.1/dt):th_x])
    v_peak = np.max(vm[th_x:th_x+int(5./dt)])
    x_peak = np.where(vm[th_x:th_x+int(5./dt)] == v_peak)[0][0]
    if len(spike_times) > 1:
        end = max(th_x+x_peak, int((spike_times[1] - 5.) / dt) - start)
    else:
        end = len(vm)
    v_AHP = np.min(vm[th_x+x_peak:end])
    x_AHP = np.where(vm[th_x+x_peak:end] == v_AHP)[0][0]
    AHP = v_before - v_AHP
    # if spike waveform includes an ADP before an AHP, return the value of the ADP in order to increase error function
    rising_x = np.where(dvdt[th_x+x_peak:th_x+x_peak+x_AHP] > 0.)[0]
    if rising_x.any():
        v_ADP = np.max(vm[th_x+x_peak+rising_x[0]:th_x+x_peak+x_AHP])
        ADP = v_ADP - th_v
    else:
        ADP = 0.
    return v_peak, th_v, ADP, AHP


def update_AMPA_NMDA(cell, x, param_indexes, default_params):
    """

    :param x: array. Default: value=1.711e-03, slope=6.652e-05, tau=7.908e+01
    :return:
    """
    if context.spines is False:
        cell.correct_for_spines()
    cell.modify_mech_param('apical', 'synapse', 'gmax', value=find_param_value('AMPA.g0', x, param_indexes, default_params),
                           origin='soma', syn_type='AMPA_KIN', slope=find_param_value('AMPA.slope', x, param_indexes,
                                                                                      default_params),
                           tau=find_param_value('AMPA.tau', x, param_indexes, default_params))
    cell.modify_mech_param('apical', 'synapse', 'gmax', value=find_param_value('AMPA.g0', x, param_indexes, default_params),
                           syn_type='AMPA_KIN', custom={'method': 'custom_inheritance_by_nonterm_branchord',
                                                        'branch_cutoff': 4}, replace=False)
    cell.modify_mech_param('apical', 'synapse', 'gmax', value=find_param_value('NMDA.gmax', x, param_indexes, default_params),
                           syn_type=context.NMDA_type)
    cell.modify_mech_param('apical', 'synapse', 'gamma', value=find_param_value('NMDA.gamma', x, param_indexes, default_params),
                           syn_type=context.NMDA_type)
    cell.modify_mech_param('apical', 'synapse', 'Kd', value=find_param_value('NMDA.Kd', x, param_indexes, default_params),
                           syn_type=context.NMDA_type)
    cell.modify_mech_param('apical', 'synapse', 'kin_scale', value=find_param_value('NMDA.kin_scale', x, param_indexes,
                                                                                    default_params), syn_type=context.NMDA_type)
    cell.init_synaptic_mechanisms()
    for sec_type in ['apical']:
        cell.modify_mech_param(sec_type, 'nas', 'gbar', find_param_value('dend.gbar_nas', x, param_indexes, default_params))
        cell.modify_mech_param(sec_type, 'nas', 'gbar', origin='parent',
                               slope=find_param_value('dend.gbar_nas slope', x, param_indexes, default_params),
                               min=find_param_value('dend.gbar_nas min', x, param_indexes, default_params),
                               custom={'method': 'custom_gradient_by_branch_ord', 'branch_order':
                                   find_param_value('dend.gbar_nas bo', x, param_indexes, default_params)}, replace=False)
        cell.modify_mech_param(sec_type, 'nas', 'gbar', origin='parent',
                               slope=find_param_value('dend.gbar_nas slope', x, param_indexes, default_params),
                               min=find_param_value('dend.gbar_nas min', x, param_indexes, default_params),
                               custom={'method': 'custom_gradient_by_terminal'}, replace=False)
    return cell

def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(context.rec_file_path, 'a') as f:
        context.sim.export_to_file(f)

def plot_exported_synaptic_features(processed_export_file_path):
    """

    :param processed_export_file_path: str
    :return:
    """
    with h5py.File(processed_export_file_path, 'r') as f:
        trace_type_list = [trace_type for trace_type in f]
        trace_group_list = [trace_group for trace_group in f.itervalues()]
        num_colors = context.num_clustered_syn * len(trace_type_list) * len(context.branch_list) * len(context.AP5_cond_list)
        colors = list(cm.rainbow(np.linspace(0, 1, num_colors)))
        plot_num = 0
        for branch in context.branch_list:
            for AP5_cond in context.AP5_cond_list:
                fig, axes = plt.subplots(1, len(trace_group_list))
                for g, trace_group in enumerate(trace_group_list):
                    for sim_id in trace_group[branch][AP5_cond].keys():
                        axes[g].scatter(trace_group[branch][AP5_cond][sim_id]['tvec'],
                                           trace_group[branch][AP5_cond][sim_id]['rec']['soma'],
                                           label='Syn %d' %int(sim_id))
                        plot_num += 1
                    axes[g].set_xlabel('Time (ms)')
                    axes[g].set_ylabel('Vm (mv)')
                    axes[g].set_title('%s %s %s' %(branch, AP5_cond, trace_type_list[g].split('_')[0]))
                    axes[g].legend(loc='best', frameon=False, framealpha=0.5)
                clean_axes(axes)
                fig.tight_layout()
        plt.show()
        plt.close()