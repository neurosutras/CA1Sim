__author__ = 'Grace Ng'
from specify_cells3 import *
from plot_results import *
import random
from numpy.polynomial import Polynomial
from moopgen import *

"""
Iterates through every spine and activates AMPA_KIN synapses. Allows measurement of EPSP attenuation
and kinetics.

Modify a YAML file to include parameters necessary for this script. Then, set up ipcluster and run parallel_optimize_main.py
Current YAML filepath: data/optimize_synaptic_defaults.yaml
"""

context = Context()


def setup_module_from_file(param_file_path='data/optimize_synaptic_defaults.yaml', output_dir='data',
                           rec_file_path=None, export_file_path=None, verbose=True, disp=True):
    """
    :param param_file_path: str (.yaml file path)
    :param output_dir: str (dir path)
    :param rec_file_path: str (.hdf5 file path)
    :param export_file_path: str (.hdf5 file path)
    :param verbose: bool
    :param disp: bool
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
                            % update_params_func_name)
        update_params_funcs.append(func)
    if rec_file_path is None:
        rec_file_path = output_dir + '/sim_output' + datetime.datetime.today().strftime('%m%d%Y%H%M') + \
                   '_pid' + str(os.getpid()) + '.hdf5'
    if export_file_path is None:
        export_file_path = output_dir + '/%s_%s_%s_optimization_exported_traces.hdf5' % \
                           (datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)
    x0_array = param_dict_to_array(x0, param_names)
    context.update(locals())
    context.update(kwargs)
    config_engine(update_params_funcs, param_names, default_params, rec_file_path, export_file_path, output_dir, disp,
                  **kwargs)


def config_controller(export_file_path, **kwargs):
    """

    :param export_file_path: str
    """
    context.update(locals())
    set_constants()


def config_engine(update_params_funcs, param_names, default_params, rec_file_path, export_file_path, output_dur, disp,
                  mech_file_path, neuroH5_file_path, neuroH5_index, spines, **kwargs):
    """
    :param update_params_funcs: list of function references
    :param param_names: list of str
    :param default_params: dict
    :param rec_file_path: str
    :param export_file_path: str
    :param output_dur: str (dir path)
    :param disp: bool
    :param mech_file_path: str
    :param neuroH5_file_path: str
    :param neuroH5_index: int
    :param spines: bool
    """
    neuroH5_dict = read_from_pkl(neuroH5_file_path)[neuroH5_index]
    param_indexes = {param_name: i for i, param_name in enumerate(param_names)}
    context.update(locals())
    context.update(kwargs)
    set_constants()
    setup_cell(**kwargs)


def set_constants():
    """

    :return:
    """
    seed_offset = 7. * 2E6
    num_random_syn = 5  # 30
    num_clustered_syn = 10  # 20
    branch_list = ['branch1', 'branch2'] # number of branches for which clusters of synpases will be stimulated simultaneously
                                        # to calculate compound EPSP
    group_type_list = ['random'] + branch_list
    syn_ind_legend = {'random': range(num_random_syn)}
    for i, branch in enumerate(branch_list):
        syn_ind_legend[branch] = range(num_random_syn + i*num_clustered_syn, num_random_syn + (i+1)*num_clustered_syn)
    AP5_cond_list = ['AP5', 'con']
    unitary_isi = 125. # Inter-stimulus interval for unitary sims
    num_unitary_stims = 5  # Number of unitary stims (each separated by isi) in a simulation
    equilibrate = 250.  # time to steady-state
    unitary_duration = equilibrate + num_unitary_stims * unitary_isi
    compound_duration = 450.
    compound_isi = 1.1
    trace_baseline = 10.
    stim_dur = 150.
    duration = max(equilibrate + stim_dur, unitary_duration, compound_duration)
    th_dvdt = 10.
    dt = 0.02
    v_init = -77.
    v_active = -77.
    NMDA_type = 'NMDA_KIN5'
    syn_types = ['AMPA_KIN', NMDA_type]
    local_random = random.Random()
    i_holding = {'soma': 0.}
    i_th = {'soma': 0.1}
    context.update(locals())


def setup_cell(verbose=False, cvode=False, daspk=False, **kwargs):
    """

    :param verbose: bool
    :param cvode: bool
    :param daspk: bool
    """
    cell = DG_GC(neuroH5_dict=context.neuroH5_dict, mech_file_path=context.mech_file_path,
                 full_spines=context.spines)
    context.cell = cell

    context.local_random.seed(int(context.neuroH5_index + context.seed_offset))

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
    candidate_branch_list = [apical for apical in cell.apical if
                             50. < cell.get_distance_to_node(cell.tree.root, apical) < 150. and
                             apical.sec.L > 80.]

    context.local_random.shuffle(candidate_branch_list)
    success_branches = 0
    for branch in candidate_branch_list:
        branch_syn_list = []
        if branch.synapses:
            branch_syn_list.extend([syn for syn in branch.synapses if ('AMPA_KIN' in syn._syn and
                                                                       30. <= syn.loc*branch.sec.L <= 60.)])
        if context.spines:
            for spine in branch.spines:
                branch_syn_list.extend([syn for syn in spine.synapses if ('AMPA_KIN' in syn._syn and
                                            30. <= (cell.get_distance_to_node(branch.parent, spine, loc=0.)) <= 60.)])
        if len(branch_syn_list) > context.num_clustered_syn:
            label = context.branch_list[success_branches-1]
            final_syn_dict[label] = context.local_random.sample(branch_syn_list, context.num_clustered_syn)
            success_branches += 1
        if success_branches > 1:
            break
    if success_branches < 2:
        raise Exception('Could not find enough branches that satisfy the specified requirements.')

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
    dend_loc = candidate_locs[index]
    rec_locs = {'soma': 0., 'dend': dend_loc, 'local_branch': 0.5}
    context.rec_locs = rec_locs
    rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'local_branch': dend}
    context.rec_nodes = rec_nodes

    equilibrate = context.equilibrate
    unitary_duration = context.unitary_duration
    dt = context.dt
    stim_dur = context.stim_dur

    sim = QuickSim(unitary_duration, cvode=True, daspk=daspk, dt=dt, verbose=verbose)
    sim.parameters['equilibrate'] = equilibrate
    sim.parameters['spines'] = context.spines
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=context.duration)
    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    context.sim = sim

    context.spike_output_vec = h.Vector()
    cell.spike_detector.record(context.spike_output_vec)


def update_mech_dict(x, mech_file_path):
    """

    :param x: array
    :param mech_file_path: str
    """
    for update_func in context.update_params_funcs:
        update_func(x, context)
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


def filter_unitary_EPSP_features(get_result, old_features, export=False):
    """

    :param get_result: list of dict (each dict has the results from a particular simulation)
    :param old_features: dict (empty at this point in the optimization)
    :param export: bool
    :return: dict
    """
    random_results = {}
    random_traces = {}
    for AP5_cond in context.AP5_cond_list:
        random_results[AP5_cond] = {}
        random_traces[AP5_cond] = []
    clustered_results = {}
    for branch in context.branch_list:
        clustered_results[branch] = {}
        for AP5_cond in context.AP5_cond_list:
            clustered_results[branch][AP5_cond] = {}
    for this_result in get_result:
        for syn_id, syn_result in this_result.iteritems():
            AP5_cond = syn_result['AP5_condition']
            group_type = syn_result['group_type']
            if group_type == 'random':
                random_results[AP5_cond][syn_id] = syn_result['EPSP_amp']
                random_traces[AP5_cond].append(syn_result['rec']['soma'])
            else:
                clustered_results[group_type][AP5_cond][syn_id] = syn_result
    avg_EPSP_AP5 = np.mean(random_results['AP5'].values())
    avg_EPSP_con = np.mean(random_results['con'].values())
    NMDA_contributions = []
    for syn_id in random_results['AP5'].iterkeys():
        NMDA_contributions.append((random_results['con'][syn_id] - random_results['AP5'][syn_id]) /
                                  random_results['con'][syn_id])
    avg_NMDA_contr = np.mean(NMDA_contributions)
    new_features = {'soma_EPSP_AP5': avg_EPSP_AP5, 'soma_EPSP_control': avg_EPSP_con,
                    'NMDA_contribution': avg_NMDA_contr}
    for branch in context.branch_list:
        new_features[branch+'_unitary'] = clustered_results[branch]
    if export:
        new_export_file_path = context.export_file_path.replace('.hdf5', '_processed.hdf5')
        with h5py.File(new_export_file_path, 'a') as f:
            average_traces = f.create_group('average_unitary_traces')
            for AP5_cond in context.AP5_cond_list:
                average_traces.create_group(AP5_cond, compression='gzip', compression_opts=9,
                                            data=np.mean(random_traces[AP5_cond], 0))
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


def filter_compound_EPSP_features(get_result, old_features, export=False):
    """

    :param get_result: list of dict (each dict has the results from a particular simulation)
    :param old_features: dict
    :param export: bool
    :return: dict
    """
    expected_compound_traces = {}
    expected_compound_EPSP_amp = {}
    actual_compound_traces = {}
    actual_compound_EPSP_amp = {}
    syn_group_id_list = {}
    for group_type in context.branch_list:
        expected_compound_traces[group_type], expected_compound_EPSP_amp[group_type] = \
            get_expected_compound_features(old_features[group_type + '_unitary'])
        actual_compound_traces[group_type] = {}
        actual_compound_EPSP_amp[group_type] = {}
        syn_group_id_list[group_type] = {}
        for AP5_cond in context.AP5_cond_list:
            actual_compound_traces[group_type][AP5_cond] = {}
            actual_compound_EPSP_amp[group_type][AP5_cond] = []
            syn_group_id_list[group_type][AP5_cond] = []
    for this_result in get_result:
        for num_sim, syn_group_id in enumerate(this_result):
            sim_result = this_result[syn_group_id]
            group_type = sim_result['group_type']
            AP5_condition = sim_result['AP5_condition']
            syn_group_id_list[group_type][AP5_condition].append(syn_group_id)
            actual_compound_traces[group_type][AP5_condition][syn_group_id] = {'rec': {}}
            for rec_loc in sim_result['rec']:
                actual_compound_traces[group_type][AP5_condition][syn_group_id]['rec'][rec_loc] = \
                    sim_result['rec'][rec_loc]
            actual_compound_traces[group_type][AP5_condition][syn_group_id]['tvec'] = sim_result['tvec']
            actual_compound_EPSP_amp[group_type][AP5_condition].append(np.max(sim_result['rec']['soma']))
                # This forms an array of the compound EPSPs, where each element is the EPSP from a certain number of synapses
    for group_type in syn_group_id_list:
        for AP5_condition in syn_group_id_list[group_type]:
            actual_compound_EPSP_amp[group_type][AP5_condition].sort(
                key=dict(zip(actual_compound_EPSP_amp[group_type][AP5_condition],
                             syn_group_id_list[group_type][AP5_condition])).get)
    integration_gain = {}
    initial_gain = {}
    for AP5_cond in context.AP5_cond_list:
        integration_gain[AP5_cond] = []
        initial_gain[AP5_cond] = []
    for group_type in context.branch_list:
        for AP5_condition in context.AP5_cond_list:
            actual = np.array(actual_compound_EPSP_amp[group_type][AP5_condition])
            expected = np.array(expected_compound_EPSP_amp[group_type][AP5_condition])
            ratio = np.divide(actual, expected)
            this_initial_gain = np.mean(ratio[:2])  # The ratio should be close to 1 for the first few synapses.
            fit = Polynomial.fit(expected, actual, 1, domain=(-1, 1))
            this_integration_gain = fit.coef[1]
            integration_gain[AP5_condition].append(this_integration_gain)
            initial_gain[AP5_condition].append(this_initial_gain)
    new_features = {'integration_gain_AP5': np.mean(integration_gain['AP5']),
                    'integration_gain_con': np.mean(integration_gain['con']),
                    'initial_gain_AP5': np.mean(initial_gain['AP5']),
                    'initial_gain_con': np.mean(initial_gain['con'])}
    if export:
        new_export_file_path = context.export_file_path.replace('.hdf5', '_processed.hdf5')
        with h5py.File(new_export_file_path, 'a') as f:
            for trace_dict, trace_group_label in zip([actual_compound_traces, expected_compound_traces],
                                                     ['actual_compound_traces', 'expected_compound_traces']):
                trace_group = f.create_group(trace_group_label)
                for group_type in trace_dict:
                    group = trace_group.create_group(group_type)
                    for AP5_condition in trace_dict[group_type]:
                        AP5_group = group.create_group(AP5_condition)
                        for syn_group_id in trace_dict[group_type][AP5_condition]:
                            syn_group = AP5_group.create_group(str(syn_group_id))
                            syn_group.create_group('rec')
                            syn_group.create_dataset('tvec', compression='gzip', compression_opts=9,
                                                     data=trace_dict[group_type][AP5_condition][syn_group_id]['tvec'])
                            for rec_loc in trace_dict[group_type][AP5_condition][syn_group_id]['rec']:
                                syn_group['rec'].create_dataset(rec_loc, compression='gzip', compression_opts=9, data=
                                    trace_dict[group_type][AP5_condition][syn_group_id]['rec'][rec_loc])
            for EPSP_amp_dict, group_label in zip([actual_compound_EPSP_amp, expected_compound_EPSP_amp],
                                               ['actual_compound_EPSP_amp', 'expected_compound_EPSP_amp']):
                EPSP_amp_group = f.create_group(group_label)
                for group_type in trace_dict:
                    group = EPSP_amp_group.create_group(group_type)
                    for AP5_condition in trace_dict[group_type]:
                        group.create_dataset(AP5_condition, compression='gzip', compression_opts=9,
                                                 data=EPSP_amp_dict[group_type][AP5_condition])
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
    :param objective_names: list of str
    :param target_val: dict of float
    :param target_range: dict of float
    :return: tuple of dict
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

    context.cell.reinit_mechanisms(from_file=True)
    if not context.spines:
        context.cell.correct_g_pas_for_spines()
    for update_func in context.update_params_funcs:
        update_func(x, context)

    test_syns = np.array(test_syns)
    soma_vm = offset_vm('soma', context.v_init)
    synapses = [context.syn_list[syn_index] for syn_index in test_syns]
    for i, syn in enumerate(synapses):
        syn.source.play(h.Vector([context.equilibrate + i * context.unitary_isi]))
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
    context.sim.run(context.v_init)
    t = np.array(context.sim.tvec)
    dt = context.dt
    equilibrate = context.equilibrate
    interp_t = np.arange(0, duration, dt)
    trace_baseline = context.trace_baseline
    corrected_t = interp_t[int((equilibrate - trace_baseline)/dt):int((equilibrate + context.unitary_isi) / dt)] - \
                  equilibrate - trace_baseline
    result = {}
    for i, syn in enumerate(synapses):
        start = int((equilibrate + i * context.unitary_isi) / dt)
        end = start + int(context.unitary_isi / dt)
        trace_start = int(start - trace_baseline / dt)
        baseline_start, baseline_end = int(start - 3. / dt), int(start - 1. / dt)
        distance = context.cell.get_distance_to_node(context.cell.tree.root, syn.node, syn.loc)
        syn_id = test_syns[i]
        if syn_id not in result:
            result[syn_id] = {'AP5_condition': AP5_condition, 'group_type': group_type, 'distance': distance,
                              'rec': {}}
        for rec in context.sim.rec_list:
            vm = np.array(rec['vec'])
            interp_vm = np.interp(interp_t, t, vm)
            baseline = np.mean(interp_vm[baseline_start:baseline_end])
            corrected_vm = interp_vm - baseline
            peak = np.max(corrected_vm[start:end])
            peak_index = np.where(corrected_vm == peak)[0][0]
            zero_index = np.where(corrected_vm[peak_index:] <= 0.)[0]
            if np.any(zero_index):
                corrected_vm[peak_index+zero_index[0]:] = 0.
            if rec['description'] == 'soma':
                result[syn_id]['EPSP_amp'] = peak
            result[syn_id]['rec'][rec['description']] = corrected_vm[trace_start:end]
            if 'tvec' not in result[syn_id]:
                result[syn_id]['tvec'] = corrected_t  # sets each tvec to start at -trace_baseline
        syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while
                                    # keeping the VecStim source object in existence so it can be activated again
    print 'Process:', os.getpid(), 'Synapse group: %s stimulated synapses: %s' % (group_type, test_syns), \
        'on Node', syn.node.name, 'with distance: %.2f in %.3f s' % (distance, time.time() - start_time)
    if plot:
        context.sim.plot()
    if export:
        export_sim_results()
    return result


def compute_branch_cooperativity_features(x, test_syns, AP5_condition, group_type, export=False, plot=False):
    """
    This simulation tests a value for gmax_NMDA_KIN. Results are later compared to branch cooperativity data from
    Krueppel et al., 2011. One apical oblique ~150 um from the soma is chosen, up to 20 spines within a 30 um long
    stretch of branch are stimulated. Actual and Expected EPSP amplitudes are compared.

    :param x: arr
    :param test_syns: list of int
    :param AP5_condition: str
    :param group_type: str
    :param export: bool
    :param plot: bool
    :return: dict
    """
    start_time = time.time()

    context.cell.reinit_mechanisms(from_file=True)
    if not context.spines:
        context.cell.correct_g_pas_for_spines()
    for update_func in context.update_params_funcs:
        update_func(x, context)

    test_syns = np.array(test_syns)
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
    trace_baseline = context.trace_baseline
    corrected_t = interp_t[int((equilibrate - trace_baseline) / dt):] - trace_baseline
    baseline_start, baseline_end = int((equilibrate - 3.0) / dt), int((equilibrate - 1.0) / dt)
    trace_start = int((equilibrate - trace_baseline) / dt)
    trace_end = int(duration / dt)
    syn_group_id = test_syns[-1]
    result = {syn_group_id: {'AP5_condition': AP5_condition, 'group_type': group_type, 'tvec': corrected_t, 'rec': {}}}
    for rec in context.sim.rec_list:
        vm = np.array(rec['vec'])
        interp_vm = np.interp(interp_t, t, vm)
        baseline = np.mean(interp_vm[baseline_start:baseline_end])
        corrected_vm = interp_vm - baseline
        peak = np.max(corrected_vm[trace_start:trace_end])
        peak_index = np.where(corrected_vm == peak)[0][0]
        zero_index = np.where(corrected_vm[peak_index:] <= 0.)[0]
        if np.any(zero_index):
            corrected_vm[peak_index + zero_index[0]:] = 0.
        result[syn_group_id]['rec'][rec['description']] = corrected_vm[trace_start:trace_end]
    for syn in synapses:
        syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while
                                    # keeping the VecStim source object in existence so it can be activated again
    print 'Process: %i, Synapse group: %s, Condition: %s, stimulated %i synapses in %i s with gmax: %.3E' % \
          (os.getpid(), group_type, AP5_condition, len(test_syns), time.time() - start_time,
           x[context.param_indexes['NMDA.gmax']])
    if plot:
        context.sim.plot()
    if export:
        export_sim_results()
    return result


def compute_stability_features(x, export=False, plot=False):
    """
    :param x: array
    :param export: bool
    :param plot: bool
    :return: dict
    """
    start_time = time.time()

    context.cell.reinit_mechanisms(from_file=True)
    if not context.spines:
        context.cell.correct_g_pas_for_spines()
    for update_func in context.update_params_funcs:
        update_func(x, context)

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
    result = {'dend_AP_amp': (dend_peak - dend_pre) / (peak - soma_vm)}
    print 'Process %i took %.1f s to find spike rheobase at amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    if plot:
        context.sim.plot()
    if export:
        export_sim_results()
    return result


def get_expected_compound_features(unitary_branch_results):
    """

    :param unitary_branch_features: dict {'AP5': {syn_id: {'rec': {loc: trace}}
    :return: dict {'AP5': {syn_group_id: {loc: summed_traces, ...}, syn_group_id: ...},
                   'con': {syn_group_id: {loc: summed_traces, ..}}}
    """
    expected_compound_traces = {}
    expected_compound_EPSP_amp = {}
    syn_group_id_list = {}
    baseline_len = int(context.trace_baseline / context.dt)
    unitary_len = int(context.unitary_isi / context.dt)
    compound_isi_len = int(context.compound_isi / context.dt)
    actual_trace_len = int((context.compound_duration - context.equilibrate) / context.dt) + baseline_len
    for AP5_condition in context.AP5_cond_list:
        expected_compound_traces[AP5_condition] = {}
        expected_compound_EPSP_amp[AP5_condition] = []
        syn_group_id_list[AP5_condition] = []
        summed_traces = {}
        start_ind = baseline_len
        for num_syn, syn_id in enumerate(unitary_branch_results[AP5_condition]):
            syn_result = unitary_branch_results[AP5_condition][syn_id]
            for rec_loc in unitary_branch_results[AP5_condition][syn_id]['rec']:
                if rec_loc not in summed_traces:
                    summed_traces[rec_loc] = np.zeros(actual_trace_len)
                summed_traces[rec_loc][start_ind:start_ind+unitary_len] += syn_result['rec'][rec_loc][baseline_len:]
            expected_compound_traces[AP5_condition][syn_id] = {}
            expected_compound_traces[AP5_condition][syn_id]['rec'] = copy.deepcopy(summed_traces) # This stores the summed traces of
                                                                                    # all the synapses up to and including
                                                                                    # the synapse with this syn_id
            expected_compound_EPSP_amp[AP5_condition].append(np.max(summed_traces['soma']))
            syn_group_id_list[AP5_condition].append(syn_id)
            start_ind += compound_isi_len
        expected_compound_EPSP_amp[AP5_condition].sort(key=dict(zip(expected_compound_EPSP_amp[AP5_condition],
                                                                syn_group_id_list[AP5_condition])).get)
    return expected_compound_traces, expected_compound_EPSP_amp


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
    duration = context.duration

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
    :param spike_times: array
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
        end = max(th_x + x_peak + int(2./dt), int((spike_times[1] - 4.) / dt) - start)
    else:
        end = len(vm)
    v_AHP = np.min(vm[th_x+x_peak:end])
    x_AHP = np.where(vm[th_x+x_peak:end] == v_AHP)[0][0]
    AHP = v_before - v_AHP
    # if spike waveform includes an ADP before an AHP, return the value of the ADP in order to increase error function
    ADP = 0.
    rising_x = np.where(dvdt[th_x+x_peak+1:th_x+x_peak+x_AHP-1] > 0.)[0]
    if rising_x.any():
        v_ADP = np.max(vm[th_x+x_peak+1+rising_x[0]:th_x+x_peak+x_AHP])
        pre_ADP = np.mean(vm[th_x+x_peak+1+rising_x[0] - int(0.1/dt):th_x+x_peak+1+rising_x[0]])
        ADP += v_ADP - pre_ADP
    falling_x = np.where(dvdt[th_x + x_peak + x_AHP + 1:end] < 0.)[0]
    if falling_x.any():
        v_ADP = np.max(vm[th_x + x_peak + x_AHP + 1: th_x + x_peak + x_AHP + 1 + falling_x[0]])
        ADP += v_ADP - v_AHP
    return v_peak, th_v, ADP, AHP


def update_AMPA_NMDA(x, local_context=None):
    """

    :param x: array. Default: value=1.711e-03, slope=6.652e-05, tau=7.908e+01
    :param local_context: :class:'Context'
    """
    if local_context is None:
        local_context = context
    cell = local_context.cell
    param_indexes = local_context.param_indexes
    cell.modify_mech_param('apical', 'synapse', 'gmax',
                           value=x[param_indexes['AMPA.g0']], origin='soma', syn_type='AMPA_KIN',
                           slope=x[param_indexes['AMPA.slope']], tau=x[param_indexes['AMPA.tau']])
    cell.modify_mech_param('apical', 'synapse', 'gmax', origin='parent', syn_type='AMPA_KIN',
                           custom={'method': 'custom_inherit_by_branch_order', 'branch_order':
                               x[param_indexes['AMPA.gmax bo']]}, replace=False)
    cell.modify_mech_param('apical', 'synapse', 'gmax', value=x[param_indexes['NMDA.gmax']],
                           syn_type=local_context.NMDA_type)
    cell.modify_mech_param('apical', 'synapse', 'gamma', value=x[param_indexes['NMDA.gamma']],
                           syn_type=local_context.NMDA_type)
    cell.modify_mech_param('apical', 'synapse', 'Kd', value=x[param_indexes['NMDA.Kd']],
                           syn_type=local_context.NMDA_type)
    cell.modify_mech_param('apical', 'synapse', 'kin_scale', value=x[param_indexes['NMDA.kin_scale']],
                           syn_type=local_context.NMDA_type)
    # cell.init_synaptic_mechanisms()
    for sec_type in ['apical']:
        cell.modify_mech_param(sec_type, 'nas', 'gbar', x[param_indexes['dend.gbar_nas']])
        cell.modify_mech_param(sec_type, 'nas', 'gbar', origin='parent', slope=x[param_indexes['dend.gbar_nas slope']],
                               min=x[param_indexes['dend.gbar_nas min']],
                               custom={'method': 'custom_gradient_by_branch_ord', 'branch_order':
                                   x[param_indexes['dend.gbar_nas bo']]}, replace=False)
        cell.modify_mech_param(sec_type, 'nas', 'gbar', origin='parent',
                               slope=x[param_indexes['dend.gbar_nas slope']], min=x[param_indexes['dend.gbar_nas min']],
                               custom={'method': 'custom_gradient_by_terminal'}, replace=False)


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
        trace_type_list = ['actual_compound_traces', 'expected_compound_traces']
        trace_group_list = [f[trace_type] for trace_type in trace_type_list]
        num_colors = context.num_clustered_syn * len(trace_type_list) * len(context.branch_list) * len(context.AP5_cond_list)
        colors = list(cm.rainbow(np.linspace(0, 1, num_colors)))
        for branch in context.branch_list:
            for AP5_cond in context.AP5_cond_list:
                fig, axes = plt.subplots(1, len(trace_group_list))
                for g, trace_group in enumerate(trace_group_list):
                    for syn_group_id in trace_group[branch][AP5_cond].keys():
                        axes[g].scatter(trace_group[branch][AP5_cond][syn_group_id]['tvec'],
                                           trace_group[branch][AP5_cond][syn_group_id]['rec']['soma'],
                                           label='Syn %d' %int(syn_group_id))
                    axes[g].set_xlabel('Time (ms)')
                    axes[g].set_ylabel('Vm (mv)')
                    axes[g].set_title('%s %s %s' %(branch, AP5_cond, trace_type_list[g].split('_')[0]))
                    axes[g].legend(loc='best', frameon=False, framealpha=0.5)
                clean_axes(axes)
                fig.tight_layout()
        for branch in context.branch_list:
            for AP5_cond in context.AP5_cond_list:
                fig, axes = plt.subplots(1)
                axes.scatter(f['expected_compound_EPSP_amp'][branch][AP5_cond],
                             f['actual_compound_EPSP_amp'][branch][AP5_cond])
                axes.set_xlabel('Expected Compound EPSP Amp')
                axes.set_ylabel('Actual Compound EPSP Amp')
                axes.set_title('%s %s EPSPs' %(branch, AP5_cond))
                clean_axes(axes)
                fig.tight_layout()
        plt.show()
        plt.close()
