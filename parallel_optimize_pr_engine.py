__author__ = 'Aaron D. Milstein'
from specify_cells import *
import random
import os
from IPython.parallel.util import interactive
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which synapse to optimize (coarse sampling of the full set of spines).
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
#mech_filename = '052915 pas_exp_scale kdr ka_scale ih_sig_scale ampar_exp_scale nmda - EB2'
mech_filename = '072815 optimized basal ka_scale dend_sh_ar_nas - ampa_scale - EB2'

param_names = ['n', 'P0', 'f', 'tau_F', 'd1', 'tau_D1']
x = []  # placeholder for optimization parameters, must be pushed to each engine at each basinhopping step
num_stims = 3
interp_dt = 0.05

@interactive
def sim_stim_train(ISI):
    """
    Called by controller, mapped to each engine. Simultaneously activates n synapses in the syn_list at the specified
    ISI and returns the mean peak amplitude. If ISI == 300, this method parses the recorded vector into three units and
    returns the mean unit size across runs.
    :param ISI: int
    :return: np.array
    """
    #x = x['x']
    #x[0] = int(1000 * x[0])
    param_dict = {}
    for i in range(len(x)):
        param_dict[param_names[i]] = x[i]
    param_dict['n'] = int(param_dict['n'] * 1000)
    #return param_dict
    local_random.seed(0)
    duration = equilibrate + ISI * (num_stims - 1) + 101
    interp_t = np.arange(0, duration, interp_dt)
    sim.tstop = duration
    spike_times = h.Vector([equilibrate + ISI * i for i in range(num_stims)])
    stim_syn_list = [syn_list[i] for i in local_random.sample(range(len(syn_list)), param_dict['n'])]
    for syn in stim_syn_list:
        syn.source.play(spike_times)
        for i, param_name in enumerate(param_names):
            if not param_name == 'n':
                setattr(syn.target('Pr'), param_name, x[i])
    start_time = time.time()
    sim.run(v_init)
    for syn in stim_syn_list:
        syn.source.play(h.Vector())
    t = np.array(sim.tvec)
    left, right = time2index(t, equilibrate-2.0, equilibrate)
    vm = np.array(sim.rec_list[0]['vec'])
    baseline = np.average(vm[left:right])
    vm -= baseline
    rec = np.interp(interp_t, t, vm)
    print 'Process:', os.getpid(), 'ISI:', ISI, 'synapses:', param_dict['n'], 'took', time.time() - start_time, 's'
    left, right = time2index(interp_t, equilibrate-2.0, duration)
    return {ISI: rec[left:right]}


def restore_random_sequence_locations():
    """
    Restores the random object for each synapse to the initial position of its sequence. Should be done in between
    basinhopping steps.
    """
    for i, syn in enumerate(syn_list):
        syn.randObj.seq(rand_seq_locs[i])


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
duration = 550.
v_init = -67.
syn_types = ['AMPA_KIN', 'NMDA_KIN2']

syn_list = []
rand_seq_locs = []
cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
zero_na()

local_random = random.Random()
for branch in cell.trunk+cell.apical:
    for node in branch.spines:
        syn = Synapse(cell, node, syn_types, stochastic=1)
        syn_list.append(syn)
        rand_seq_locs.append(syn.randObj.seq())
cell.init_synaptic_mechanisms()
sim = QuickSim(duration, verbose=0)

# look for a trunk bifurcation
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
sim.append_rec(cell, trunk, description='trunk', loc=0.)