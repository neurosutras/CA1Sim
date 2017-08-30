__author__ = 'Aaron D. Milstein' and 'Grace Ng'
from specify_cells2 import *
import random
import os
import sys
from ipyparallel import interactive
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which synapse to optimize (coarse sampling of the full set of spines).
"""
neurotree_filename = '121516_DGC_trees.pkl'
neurotree_dict = read_from_pkl(morph_dir+neurotree_filename)

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    mech_filename = '042617 GC optimizing spike stability'


@interactive
def offset_vm(description, vm_target=None):
    """

    :param description: str
    :param vm_target: float
    """
    if vm_target is None:
        vm_target = v_init
    node = rec_nodes[description]
    loc = rec_locs[description]
    rec_dict = sim.get_rec(description)
    sim.modify_stim(0, node=node, loc=loc)
    rec = rec_dict['vec']
    offset = True
    sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    sim.modify_stim(0, amp=i_holding[description])
    sim.run(vm_target)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < vm_target - 0.5:
        i_holding[description] += 0.01
        while offset:
            if sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(0, amp=i_holding[description])
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
            sim.modify_stim(0, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                i_holding[description] -= 0.01
            else:
                offset = False
    sim.tstop = duration
    return v_rest


def epsp_amp_error(x, syn, syn_index):
    """
    Function called by optimize.minimize. Sets specified synaptic point_process parameters, runs a simulation
    stimulating one synapse, and calculates error based on distance from target amplitude of resulting somatic
    EPSP.
    :param x: list of parameters
    :param syn: :class:'Synapse'
    :return: float: error
    """
    for i in range(len(x)):
        setattr(syn.target(syn_type), param_names[i], x[i])
    start_time = time.time()
    if spines:
        spine = syn.node
        branch = spine.parent.parent
        sim.modify_rec(1, branch)
        sim.modify_rec(2, spine)
    else:
        branch = syn.node
        sim.modify_rec(1, branch)
        sim.modify_rec(2, branch, loc=syn.loc)
    sim.parameters['input_loc'] = branch.type
    sim.parameters['syn_type'] = syn_type
    for index, param_name in enumerate(param_names):
        sim.parameters[param_name] = x[index]
    sim.run()
    t = np.array(sim.tvec)
    vm = np.array(sim.rec_list[0]['vec'])
    interp_t = np.arange(0, duration, dt)
    interp_vm = np.interp(interp_t, t, vm)
    left, right = time2index(interp_t, equilibrate-3.0, equilibrate-1.0)
    baseline = np.average(interp_vm[left:right])
    start, end = time2index(interp_t, equilibrate, duration)
    amp = np.max(interp_vm[start:end]) - baseline
    result = {'EPSP_amp': amp}
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print 'Process:', os.getpid(), 'Synapse:', syn_index, 'Node:', syn.node.name, 'Time: %.3f s, x: ' \
                                                            '%.2E, Amp: %.3f, Error: %.2E' % (time.time() - start_time,
                                                            x[0], amp, Err)
    return Err

@interactive
def optimize_single_synapse(syn_index):
    """
    Called by controller, mapped to each engine. Runs optimization procedure for a single synapse, returns the optimized
    parameters, distance of the spine from the soma, and the sec_type of the associated dendritic branch.
    :param syn_index: str
    :return: dict
    """
    start_time = time.time()
    soma_vm = offset_vm('soma', v_init)
    syn = syn_list[syn_index]
    syn.source.play(spike_times)
    #result = optimize.minimize(epsp_amp_error, x0, method='L-BFGS-B', args=(syn,), options={'ftol': 1e-3},
    #                           bounds=xbounds)
    # options={'maxfun': 25}
    result = optimize.minimize(epsp_amp_error, x0, method='Nelder-Mead', args=(syn, syn_index), options={'xatol': 1e-7,
                                                                                    'fatol': 1e-3, 'maxiter': 20})
    epsp_amp_error(result.x, syn, syn_index)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, syn_index)
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    distance = cell.get_distance_to_node(cell.tree.root, syn.node, syn.loc)
    print 'Process:', os.getpid(), 'optimized synapse:', syn_index, 'on Node:', syn.node.name, 'distance: %.2f in ' \
                                                '%.3f s, x: %.2E, after %i iterations with Err: %.2E' % \
                                                (distance, time.time() - start_time, result.x[0], result.nfev, result.fun)
    return rec_filename
    #param_vals = [p for p in result.x]
    #return {'rec_filename': rec_filename, 'distance': distance, 'result': param_vals, 'sec_type': syn.node.type}


equilibrate = 250.  # time to steady-state
duration = 300.
v_init = -77.
dt = 0.02
syn_type = 'AMPA_KIN'
param_names = ['gmax']
param_ylabels = ['Peak Conductance (uS)']
local_random = random.Random()

syn_list = []
cell = DG_GC(neurotree_dict=neurotree_dict[0], mech_filename=mech_filename, full_spines=spines)
if not spines:
    cell.correct_for_spines()

cell.zero_na()
local_random.seed(0)

for branch in cell.apical:
    if spines:
        if len(branch.spines) > 1:
            if branch.sec.L <= 10.:
                node = branch.spines[local_random.sample(range(0, len(branch.spines)), 1)[0]]
                syn = Synapse(cell, node, [syn_type], stochastic=0)
                syn_list.append(syn)
            else:
                num_syns = min(len(branch.spines), int(branch.sec.L // 10.))  # a random synapse every 10 um
                for i in local_random.sample(range(0, len(branch.spines)), num_syns):
                    node = branch.spines[i]
                    syn = Synapse(cell, node, [syn_type], stochastic=0)
                    syn_list.append(syn)
        elif branch.spines:
            node = branch.spines[0]
            syn = Synapse(cell, node, [syn_type], stochastic=0)
            syn_list.append(syn)
    else:
        syn_loc_list = branch.synapse_locs['excitatory']
        if len(syn_loc_list) > 1:
            if branch.sec.L <= 10.:
                syn_loc = local_random.sample(syn_loc_list, 1)
                syn = Synapse(cell, branch, [syn_type], loc=syn_loc, stochastic=0)
                syn_list.append(syn)
            else:
                num_syns = min(len(syn_loc_list), int(branch.sec.L//10.))  # a random synapse every 10 um
                for syn_loc in local_random.sample(syn_loc_list, num_syns):
                    syn = Synapse(cell, branch, [syn_type], loc=syn_loc, stochastic=0)
                    syn_list.append(syn)
        elif syn_loc_list:
            syn_loc = syn_loc_list[0]
            syn = Synapse(cell, branch, [syn_type], loc=syn_loc, stochastic=0)
            syn_list.append(syn)

cell.init_synaptic_mechanisms()
sim = QuickSim(duration, verbose=0)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['duration'] = duration

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
axon_seg_locs = [seg.x for seg in cell.axon[2].sec]

rec_locs = {'soma': 0., 'dend': 0.}
rec_nodes = {'soma': cell.tree.root, 'dend': dend}

sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

for description, node in rec_nodes.iteritems():
    sim.append_rec(cell, node, loc=rec_locs[description], description=description)
sim.append_rec(cell, dend, loc=0., description='synapse')

spike_times = h.Vector([equilibrate])

i_holding = {'soma': 0.}


#the target values and acceptable ranges
target_val = {'EPSP_amp': 0.6}
target_range = {'EPSP_amp': 0.01}

#the initial guess
# x = [gmax]
x0 = [0.002]
