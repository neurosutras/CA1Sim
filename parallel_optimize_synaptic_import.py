__author__ = 'Grace Ng'
from specify_cells3 import *
from plot_results import *

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
    seed_offset = 1
    equilibrate = 250.  # time to steady-state
    duration = 450.
    dt = 0.02
    v_init = -77.
    NMDA_type = 'NMDA_KIN5'
    syn_types = ['AMPA_KIN']
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
    context.local_random.seed(context.seed)
    final_syn_list = {}
    random_syn_list = []
    for branch in cell.apical:
        if context.spines:
            if len(branch.spines) > 1:
                if branch.sec.L <= 10.:
                    node = branch.spines[context.local_random.sample(range(0, len(branch.spines)), 1)[0]]
                    syn = Synapse(cell, node, [context.syn_type], stochastic=0)
                    random_syn_list.append(syn)
                else:
                    num_syns = min(len(branch.spines), int(branch.sec.L // 10.))  # a random synapse every 10 um
                    for i in context.local_random.sample(range(0, len(branch.spines)), num_syns):
                        node = branch.spines[i]
                        syn = Synapse(cell, node, context.syn_types, stochastic=0)
                        random_syn_list.append(syn)
            elif branch.spines:
                node = branch.spines[0]
                syn = Synapse(cell, node, context.syn_types, stochastic=0)
                random_syn_list.append(syn)
        else:
            syn_loc_list = branch.synapse_locs['excitatory']
            if len(syn_loc_list) > 1:
                if branch.sec.L <= 10.:
                    syn_loc = context.local_random.sample(syn_loc_list, 1)
                    syn = Synapse(cell, branch, context.syn_types, loc=syn_loc, stochastic=0)
                    random_syn_list.append(syn)
                else:
                    num_syns = min(len(syn_loc_list), int(branch.sec.L // 10.))  # a random synapse every 10 um
                    for syn_loc in context.local_random.sample(syn_loc_list, num_syns):
                        syn = Synapse(cell, branch, context.syn_types, loc=syn_loc, stochastic=0)
                        random_syn_list.append(syn)
            elif syn_loc_list:
                syn_loc = syn_loc_list[0]
                syn = Synapse(cell, branch, context.syn_types, loc=syn_loc, stochastic=0)
                random_syn_list.append(syn)
    random_syn_list = context.local_random.sample()
    cell.init_synaptic_mechanisms()
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

    rec_locs = {'soma': 0., 'dend': 0.}
    context.rec_locs = rec_locs
    rec_nodes = {'soma': cell.tree.root, 'dend': dend}
    context.rec_nodes = rec_nodes

    equilibrate = context.equilibrate
    duration = context.duration
    dt = context.dt

    sim = QuickSim(duration, cvode=False, dt=dt, verbose=False)
    sim.parameters['equilibrate'] = equilibrate
    sim.parameters['duration'] = duration
    sim.parameters['spines'] = context.spines
    sim.append_stim(cell, cell.tree.root, loc  =0., amp=0., delay=0., dur=duration)
    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    sim.append_rec(cell, dend, loc=0., description='synapse')
    context.sim = sim
    context.spike_times = h.Vector([equilibrate])

def update_mech_dict(x, update_function, mech_file_path):
    update_function(x)
    context.cell.export_mech_dict(mech_file_path)

def get_EPSP_amp_features(indiv, c, client_range, export=False):
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
    # need to somehow get length of synapse list
    setup_cell()
    num_synapses = len(context.syn_list)
    result = dv.map_async(compute_EPSP_amp_features, [x] * num_synapses, range(num_synapses), [export] * num_synapses)
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result}

def get_branch_cooperativity_features(indiv, c, client_range, export=False):
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
    # need to somehow get length of synapse list
    setup_cell()
    num_synapses = len(context.syn_list)
    result = dv.map_async(compute_branch_cooperativity_features, [x] * num_synapses, range(num_synapses), [export] * num_synapses)
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result}

def get_objectives(features, objective_names, target_val, target_range):
    """

    :param features: dict
    :return: dict
    """
    if features is None: #No rheobase value found
        objectives = None
    else:
        objectives = {}
        for objective in objective_names:
            objectives[objective] = 0.
        rheobase = features['rheobase']
        for target in ['v_rest', 'v_th', 'ADP', 'AHP', 'spont_firing', 'rebound_firing', 'vm_stability', 'ais_delay',
                       'slow_depo', 'dend_amp', 'soma_peak', 'th_count', 'na_gradient']:
            # don't penalize AHP or slow_depo less than target
            if not ((target == 'AHP' and features[target] < target_val[target]) or
                        (target == 'slow_depo' and features[target] < target_val[target])):
                objectives[target] = ((target_val[target] - features[target]) / target_range[target]) ** 2.
        for i, this_adi in enumerate(features['adi']):
            if this_adi is not None and features['exp_adi'] is not None:
                objectives['adi'] += ((this_adi - features['exp_adi'][i]) / (0.01 * features['exp_adi'][i])) ** 2.
        features.pop('exp_adi')
        all_adi = []
        for adi in features['adi']:
            if adi is not None:
                all_adi.append(adi)
        features['adi'] = np.mean(all_adi)
        num_increments = context.num_increments
        i_inj_increment = context.i_inj_increment
        target_f_I = [context.experimental_f_I_slope * np.log((rheobase + i_inj_increment * (i + 1)) / rheobase)
                      for i in range(num_increments)]
        f_I_residuals = [(features['f_I'][i] - target_f_I[i]) for i in range(num_increments)]
        features['f_I_residuals'] = np.mean(f_I_residuals)
        for i in range(num_increments):
            objectives['f_I_slope'] += (f_I_residuals[i] / (0.01 * target_f_I[i])) ** 2.
        I_inj = [np.log((rheobase + i_inj_increment * (i + 1)) / rheobase) for i in range(num_increments)]
        slope, intercept, r_value, p_value, std_err = stats.linregress(I_inj, features['f_I'])
        features['f_I_slope'] = slope
        features.pop('f_I')
    return features, objectives

def compute_EPSP_amp_features(x, syn_index, syn_type, export=False, plot=False):
    """
    :param local_x: array [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor, soma.sh_nas/x,
                    axon.gkbar factor, dend.gkabar factor]
    :param plot: bool
    :return: float
    """
    start_time = time.time()
    if context.prev_job_type == 'EPSP_amp':
        context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    else:
        init_synaptic_engine()
        setup_cell()
    update_na_ka_stability(x)
    syn = context.syn_list[syn_index]
    syn.source.play(context.spike_times)
    for i in range(len(x)):
        setattr(syn.target(syn_type), context.param_names[i], x[i]) #Should this go in update?
    # sim.cvode_state = True
    if context.spines:
        spine = syn.node
        branch = spine.parent.parent
        context.sim.modify_rec(1, branch)
        context.sim.modify_rec(2, spine)
    else:
        branch = syn.node
        context.sim.modify_rec(1, branch)
        context.sim.modify_rec(2, branch, loc=syn.loc)
    context.sim.parameters['input_loc'] = branch.type
    context.sim.parameters['syn_type'] = syn_type
    for index, param_name in enumerate(context.param_names):
        context.sim.parameters[param_name] = x[index]
    context.sim.run()
    t = np.array(context.sim.tvec)
    vm = np.array(context.sim.rec_list[0]['vec'])

    duration = context.duration
    dt = context.dt
    equilibrate = context.equilibrate

    interp_t = np.arange(0, duration, dt)
    interp_vm = np.interp(interp_t, t, vm)
    left, right = time2index(interp_t, equilibrate-3.0, equilibrate-1.0)
    baseline = np.average(interp_vm[left:right])
    start, end = time2index(interp_t, equilibrate, duration)
    amp = np.max(interp_vm[start:end]) - baseline
    result = {'EPSP_amp': amp}
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    distance = context.cell.get_distance_to_node(context.cell.tree.root, syn.node, syn.loc)
    print 'Process:', os.getpid(), 'optimized synapse:', syn_index, 'on Node:', syn.node.name, 'distance: %.2f in ' \
                                                '%.3f s, x: %.2E, after %i iterations with Err: %.2E' % \
                                                (distance, time.time() - start_time, result.x[0], result.nfev, result.fun)
    if plot:
        context.sim.plot()
    if export:
        export_sim_results()
    context.prev_job_type = 'EPSP_amp'
    return result

def compute_branch_cooperativity_features(x, syn_index, export=False, plot=False):
    """
    :param local_x: array [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor, soma.sh_nas/x,
                    axon.gkbar factor, dend.gkabar factor]
    :param plot: bool
    :return: float
    """
    start_time = time.time()
    if context.prev_job_type == 'EPSP_attn':
        context.cell.reinit_mechanisms(reset_cable=True, from_file=True)
    else:
        init_synaptic_engine()
        setup_cell()
    update_na_ka_stability(x)
    # sim.cvode_state = True


    if plot:
        context.sim.plot()
    if export:
        export_sim_results()
    context.prev_job_type = 'EPSP_attn'
    return result



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

def update_na_ka_stability(x):
    """

    :param x: array ['soma.gbar_nas', 'dend.gbar_nas', 'axon.gbar_nax', 'ais.gbar_nax', 'soma.gkabar', 'dend.gkabar',
                       'soma.gkdrbar', 'axon.gkbar', 'soma.sh_nas/x', 'ais.sha_nas', 'soma.gCa factor',
                       'soma.gCadepK factor', 'soma.gkmbar', 'ais.gkmbar']
    """
    param_indexes = context.param_indexes
    context.cell.modify_mech_param('soma', 'nas', 'gbar', x[param_indexes['soma.gbar_nas']])
    context.cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[param_indexes['soma.gkdrbar']])
    context.cell.modify_mech_param('soma', 'kap', 'gkabar', x[param_indexes['soma.gkabar']])
    slope = (x[param_indexes['dend.gkabar']] - x[param_indexes['soma.gkabar']]) / 300.
    context.cell.modify_mech_param('soma', 'nas', 'sh', x[param_indexes['soma.sh_nas/x']])
    for sec_type in ['apical']:
        context.cell.reinitialize_subset_mechanisms(sec_type, 'nas')
        context.cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        context.cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        context.cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        context.cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=x[param_indexes['soma.gkabar']]+slope*75., replace=False)
        context.cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300.,
                               value=x[param_indexes['soma.gkabar']]+slope*300., replace=False)
        context.cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
        context.cell.modify_mech_param(sec_type, 'nas', 'sha', 5.)
        context.cell.modify_mech_param(sec_type, 'nas', 'gbar', x[param_indexes['dend.gbar_nas']])
    context.cell.set_terminal_branch_na_gradient()
    context.cell.reinitialize_subset_mechanisms('axon_hill', 'kap')
    context.cell.reinitialize_subset_mechanisms('axon_hill', 'kdr')
    context.cell.modify_mech_param('ais', 'kdr', 'gkdrbar', origin='soma')
    context.cell.modify_mech_param('ais', 'kap', 'gkabar', x[param_indexes['axon.gkbar']])
    context.cell.modify_mech_param('axon', 'kdr', 'gkdrbar', origin='ais')
    context.cell.modify_mech_param('axon', 'kap', 'gkabar', origin='ais')
    context.cell.modify_mech_param('axon_hill', 'nax', 'sh', x[param_indexes['soma.sh_nas/x']])
    context.cell.modify_mech_param('axon_hill', 'nax', 'gbar', context.soma_na_gbar)
    context.cell.modify_mech_param('axon', 'nax', 'gbar', x[param_indexes['axon.gbar_nax']])
    for sec_type in ['ais', 'axon']:
        context.cell.modify_mech_param(sec_type, 'nax', 'sh', origin='axon_hill')
    context.cell.modify_mech_param('soma', 'Ca', 'gcamult', x[param_indexes['soma.gCa factor']])
    context.cell.modify_mech_param('soma', 'CadepK', 'gcakmult', x[param_indexes['soma.gCadepK factor']])
    context.cell.modify_mech_param('soma', 'km3', 'gkmbar', x[param_indexes['soma.gkmbar']])
    context.cell.modify_mech_param('ais', 'km3', 'gkmbar', x[param_indexes['ais.gkmbar']])
    context.cell.modify_mech_param('axon_hill', 'km3', 'gkmbar', origin='soma')
    context.cell.modify_mech_param('axon', 'km3', 'gkmbar', origin='ais')
    context.cell.modify_mech_param('ais', 'nax', 'sha', x[param_indexes['ais.sha_nas']])
    context.cell.modify_mech_param('ais', 'nax', 'gbar', x[param_indexes['ais.gbar_nax']])
    if context.spines is False:
        context.cell.correct_for_spines()
    context.cell.set_terminal_branch_na_gradient()

def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(context.rec_filepath, 'a') as f:
        context.sim.export_to_file(f)