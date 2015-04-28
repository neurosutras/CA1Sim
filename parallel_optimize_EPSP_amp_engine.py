__author__ = 'Aaron D. Milstein'
from specify_cells import *
import random
import os
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which synapse to optimize (coarse sampling of the full set of spines).
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'

mech_filename = '042315 pas_ih_exp_scale kdr ka_scale - EB2.pkl'


def epsp_amp_error(x, syn):
    """
    Function called by optimize.minimize. Sets specified synaptic point_process parameters, runs a simulation
    stimulating one spine synapse, and calculates error based on distance from target amplitude of resulting somatic
    EPSP.
    :param x: list of parameters
    :param syn: :class:'Synapse'
    :return: float: error
    """
    for i in range(len(x)):
        setattr(syn.target(syn_type), param_names[i], x[i])
    start_time = time.time()
    sim.run()
    t = np.array(sim.tvec)
    vm = np.array(sim.rec_list[0]['vec'])
    interp_t = np.arange(0, duration, 0.001)
    interp_vm = np.interp(interp_t, t, vm)
    left, right = time2index(interp_t, equilibrate-3.0, equilibrate-1.0)
    baseline = np.average(interp_vm[left:right])
    start, end = time2index(interp_t, equilibrate, duration)
    amp = np.max(interp_vm[start:end]) - baseline
    result = {'EPSP_amp': amp}
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print 'Process:', os.getpid(), 'Spine:', syn.node.index, 'Node:', syn.node.parent.parent.name, 'Time:', \
        time.time() - start_time, 's, x:', x, 'Amp:', amp, 'Error:', Err
    return Err


def optimize_single_synapse(syn_index):
    """
    Called by controller, mapped to each engine. Runs optimization procedure for a single spine, returns the optimized
    parameters, distance of the spine from the soma, and the sec_type of the associated dendritic branch.
    :param syn_index: str
    :return: dict
    """
    start_time = time.time()
    syn = syn_list[syn_index]
    syn.source.play(spike_times)
    #result = optimize.minimize(epsp_amp_error, x0, method='L-BFGS-B', args=(syn,), options={'ftol': 1e-3},
    #                           bounds=xbounds)
    # options={'maxfun': 25}
    result = optimize.minimize(epsp_amp_error, x0, method='Nelder-Mead', args=(syn,), options={'ftol': 1e-3})
    syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
    print 'Process:', os.getpid(), 'optimized Spine:', syn.node.index, 'on Node:', syn.node.parent.parent.name, 'in', \
        time.time() - start_time, 's. x:', result.x, 'after', result.nfev, 'iterations'
    distance = cell.get_distance_to_node(cell.tree.root, syn.node.parent.parent, syn.loc)
    param_vals = [p for p in result.x]
    return {'distance': distance, 'result': param_vals, 'sec_type': syn.node.parent.parent.type}


equilibrate = 200.  # time to steady-state
duration = 250.
v_init = -65.
syn_type = 'AMPA_KIN'
param_names = ['gmax']
param_ylabels = ['Peak Conductance (uS)']

syn_list = []
cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
random.seed(0)
for branch in cell.basal+cell.trunk+cell.apical+cell.tuft:
    if len(branch.spines) > 1:
        if branch.sec.L <= 10.:
            node = branch.spines[random.sample(range(0, len(branch.spines)), 1)[0]]
            syn = Synapse(cell, node, [syn_type], stochastic=0)
            syn_list.append(syn)
        else:
            num_syns = min(len(branch.spines), int(branch.sec.L//10.))  # a random synapse every 10 um
            for i in random.sample(range(0, len(branch.spines)), num_syns):
                node = branch.spines[i]
                syn = Synapse(cell, node, [syn_type], stochastic=0)
                syn_list.append(syn)
    elif branch.spines:
        node = branch.spines[0]
        syn = Synapse(cell, node, [syn_type], stochastic=0)
        syn_list.append(syn)
cell.init_synaptic_mechanisms()
sim = QuickSim(duration, verbose=0)
sim.append_rec(cell, cell.tree.root, 0.5, description='soma')
spike_times = h.Vector([equilibrate])

#the target values and acceptable ranges
target_val = {'EPSP_amp': 0.2}
target_range = {'EPSP_amp': 0.02}

#the initial guess
# x = [gmax]
x0 = [0.0005]

# the bounds
xmin = [0.000001]
xmax = [1.]

# rewrite the bounds in the way required by L-BFGS-B
xbounds = [(low, high) for low, high in zip(xmin, xmax)]