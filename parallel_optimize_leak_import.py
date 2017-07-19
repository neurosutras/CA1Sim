__author__ = 'Grace Ng'
from specify_cells3 import *
from plot_results import *

"""
Aims for spike initiation at initial segment by increasing nax density and decreasing activation V1/2 relative to soma,
axon_hill, and axon compartments. Extend linear kap gradient into basals and obliques, aim for 60% spike attenuation
at bifurcation of trunk and tuft.

Optimizes gbar_nax/nas/sh/sha, gkabar_kap/d, gkdrbar for target na spike threshold, AHP amp, and vm stability

Import this script into parallel_optimize_main. Then, set up ipcluster and run parallel_optimize_main.py
"""
# param_file_path = 'data/optimize_leak_defaults.yaml'

equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.02
# dt = 0.002
amp = 0.3
th_dvdt = 10.
v_init = -77.
v_active = -77.
i_holding = {'soma': 0., 'dend': 0., 'distal_dend': 0.}
soma_ek = -77.

@interactive
def init_engine(engine_param_names, engine_mech_file_path, engine_neurotree_dict, engine_spines, ind):
    """

    :param engine_param_names: list
    :param engine_mech_file_path: str (path)
    :param engine_neurotree_dict: dict
    :param engine_spines: bool
    :return:
    """
    global module
    module = feat_module_refs[ind]
    global prev_job_type
    if prev_job_type != 'leak':
        global param_names
        param_names = engine_param_names
        global mech_file_path
        mech_file_path = engine_mech_file_path
        global neurotree_dict
        neurotree_dict = engine_neurotree_dict
        global spines
        spines = engine_spines
        global param_indexes
        param_indexes = {param_name: i for i, param_name in enumerate(param_names)}

@interactive
def setup_cell():
    """

    """
    global cell
    cell = DG_GC(neurotree_dict=neurotree_dict, mech_file_path=mech_file_path, full_spines=spines)

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

    # get the most distal terminal branch > 300 um from the soma
    candidate_branches = []
    candidate_end_distances = []
    for branch in (branch for branch in cell.apical if cell.is_terminal(branch)):
        if cell.get_distance_to_node(cell.tree.root, branch, 0.) >= 300.:
            candidate_branches.append(branch)
            candidate_end_distances.append(cell.get_distance_to_node(cell.tree.root, branch, 1.))
    index = candidate_end_distances.index(max(candidate_end_distances))
    distal_dend = candidate_branches[index]
    distal_dend_loc = 1.

    global rec_locs
    rec_locs = {'soma': 0., 'dend': dend_loc, 'distal_dend': distal_dend_loc}
    global rec_nodes
    rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'distal_dend': distal_dend}

    equilibrate = module.equilibrate
    stim_dur = module.stim_dur
    duration = module.duration
    dt = module.dt

    global sim
    sim = QuickSim(duration, cvode=False, dt=dt, verbose=False)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)
    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    sim.parameters['spines'] = spines

    global spike_output_vec
    spike_output_vec = h.Vector()
    cell.spike_detector.record(spike_output_vec)

@interactive
def update_mech_dict(x, update_function, mech_file_path):
    update_function(x)
    cell.export_mech_dict(mech_file_path)

@interactive
def get_Rinp_features(indiv, c, client_range, param_names, mech_file_path, neurotree_dict, spines, ind,
                           feat_module_ref, export=False):
    """
    Distribute simulations across available engines for testing spike stability.
    :param indiv: dict {'pop_id': pop_id, 'x': x arr, 'features': features dict}
    :param client_range: list of ints
    :param export: False (for exporting voltage traces)
    :return: dict
    """
    dv = c[client_range]
    dv.map_sync(feat_module_ref.init_engine, [param_names] * len(client_range), [mech_file_path] * len(client_range),
                [neurotree_dict] * len(client_range), [spines] * len(client_range),
                [ind] * len(client_range))
    sec_list = ['soma', 'dend', 'distal_dend']
    x = indiv['x']
    result = dv.map_async(feat_module_ref.get_Rinp_for_section, sec_list, [x] * len(sec_list), [export] * len(sec_list))
    return {'pop_id': indiv['pop_id'], 'client_range': client_range, 'async_result': result}

@interactive
def get_objectives(module, features, objective_names, target_val, target_range):
    """

    :param features: dict
    :return: dict
    """
    objectives = {}
    for feature_name in ['soma R_inp', 'dend R_inp']:
        objective_name = feature_name
        objectives[objective_name] = ((target_val[objective_name] - features[feature_name]) /
                                                  target_range[objective_name]) ** 2.
    this_feature = features['distal_dend R_inp'] - features['dend R_inp']
    objective_name = 'distal_dend R_inp'
    if this_feature < 0.:
        objectives[objective_name] = (this_feature / target_range['dend R_inp']) ** 2.
    else:
        objectives[objective_name] = 0.
    return features, objectives

@interactive
def get_Rinp_for_section(section, x, export=False):
    """
    Inject a hyperpolarizing step current into the specified section, and return the steady-state input resistance.
    :param section: str
    :param x: array
    :param export: bool
    :return: dict: {str: float}
    """
    start_time = time.time()
    global prev_job_type
    if prev_job_type == 'leak':
        cell.reinit_mechanisms(reset_cable=True, from_file=True)
    else:
        module.setup_cell()
    module.update_pas_exp(x)
    cell.zero_na()

    duration = module.duration
    stim_dur = module.stim_dur
    equilibrate = module.equilibrate
    v_init = module.v_init

    sim.tstop = duration
    sim.parameters['section'] = section
    sim.parameters['target'] = 'Rinp'
    sim.parameters['optimization'] = 'pas'
    amp = -0.05
    module.offset_vm(section)
    loc = rec_locs[section]
    node = rec_nodes[section]
    rec = sim.get_rec(section)
    sim.modify_stim(0, node=node, loc=loc, amp=amp, dur=stim_dur)
    sim.run(v_init)
    Rinp = get_Rinp(np.array(sim.tvec), np.array(rec['vec']), equilibrate, duration, amp)[2]
    result = {}
    result[section+' R_inp'] = Rinp
    print 'Process:', os.getpid(), 'calculated Rinp for %s in %.1f s, Rinp: %.1f' % (section, time.time() - start_time,
                                                                                    Rinp)
    if export:
        module.export_sim_results()
    prev_job_type = 'leak'
    return result

@interactive
def offset_vm(description, vm_target=None):
    """

    :param description: str
    :param vm_target: float
    """
    if vm_target is None:
        vm_target = module.v_init
    sim.modify_stim(0, amp=0.)
    node = rec_nodes[description]
    loc = rec_locs[description]
    rec_dict = sim.get_rec(description)
    sim.modify_stim(1, node=node, loc=loc, amp=0.)
    rec = rec_dict['vec']
    offset = True

    equilibrate = module.equilibrate
    dt = module.dt
    i_holding = module.i_holding
    duration = module.duration

    sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    sim.modify_stim(1, amp=i_holding[description])
    sim.run(vm_target)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < vm_target - 0.5:
        i_holding[description] += 0.01
        while offset:
            if sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(1, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest < vm_target - 0.5:
                i_holding[description] += 0.01
            else:
                offset = False
    elif v_rest > vm_target + 0.5:
        i_holding[description] -= 0.01
        while offset:
            if sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(1, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                i_holding[description] -= 0.01
            else:
                offset = False
    sim.tstop = duration
    return v_rest

@interactive
def update_pas_exp(x):
    """

    x0 = ['soma.g_pas': 2.28e-05, 'dend.g_pas slope': 1.58e-06, 'dend.g_pas tau': 58.4]
    :param x: array [soma.g_pas, dend.g_pas slope, dend.g_pas tau]
    """
    if spines is False:
        cell.reinit_mechanisms(reset_cable=True)
    cell.modify_mech_param('soma', 'pas', 'g', x[param_indexes['soma.g_pas']])
    cell.modify_mech_param('apical', 'pas', 'g', origin='soma', slope=x[param_indexes['dend.g_pas slope']],
                           tau=x[param_indexes['dend.g_pas tau']])
    for sec_type in ['axon_hill', 'axon', 'ais', 'apical', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')
    if spines is False:
        cell.correct_for_spines()

@interactive
def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f)