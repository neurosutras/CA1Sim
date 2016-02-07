__author__ = 'Aaron D. Milstein'
from specify_cells import *
import random
import os
from IPython.parallel.util import interactive
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which synapse to optimize (coarse sampling of the full set of spines).
"""
# morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
# mech_filename = '052915 pas_exp_scale kdr ka_scale ih_sig_scale ampar_exp_scale nmda - EB2'
# mech_filename = '080615 rebalanced na_ka ampa nmda - EB2'
# mech_filename = '103115 interim dendritic excitability ampa nmda_kin3'
# mech_filename = '012816 altered intrinsic properties - ampa nmda_kin4'
mech_filename = '020516 altered km2 rinp - ampa nmda_kin5'
# rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())

param_names = ['f', 'tau_F', 'd1', 'tau_D1']
x = []  # placeholder for optimization parameters, must be pushed to each engine at each iteration
num_stims = 3
interp_dt = 0.02
P0 = 0.2  # set as static to reduce difficulty of optimization


@interactive
def sim_stim_train_dynamics(ISI):
    """
    Called by controller, mapped to each engine. Simultaneously activates n synapses in the syn_list at the specified
    ISI with the dynamics stored in x and returns the resulting EPSP waveform.
    :param ISI: int
    :return: np.array
    """
    N = len(syn_list.syns_to_stim)
    this_Pr = Pr(P0, x[0], x[1], x[2], x[3])
    stim_time_array = [[] for i in range(len(syn_list.syns_to_stim))]
    if ISI == 10:
        duration = equilibrate + ISI * (num_stims - 1) + 110. + 101.
        stim_times = [equilibrate + ISI * i for i in range(num_stims + 2)]
        stim_times.append(stim_times[-1] + 110.)
    else:
        duration = equilibrate + ISI * (num_stims - 1) + 101.
        stim_times = [equilibrate + ISI * i for i in range(num_stims)]
    # P_list = []
    for stim_time in stim_times:
        P = this_Pr.stim(stim_time)
        # P_list.append(P)
        for j in local_random.sample(range(N), int(P * N)):
            stim_time_array[j].append(stim_time)
    stim_time_vector_array = [h.Vector(stim_time_list) for stim_time_list in stim_time_array]
    interp_t = np.arange(0., duration, interp_dt)
    sim.tstop = duration
    for i, syn in enumerate(syn_list.syns_to_stim):
        syn.source.play(stim_time_vector_array[i])
    start_time = time.time()
    sim.run(v_init)
    del stim_time_vector_array
    for syn in syn_list.syns_to_stim:
        syn.source.play(h.Vector())
    t = np.array(sim.tvec)
    left, right = time2index(t, equilibrate-2.0, equilibrate)
    vm = np.array(sim.rec_list[0]['vec'])
    baseline = np.average(vm[left:right])
    vm -= baseline
    rec = np.interp(interp_t, t, vm)
    print 'Process:', os.getpid(), 'ISI:', ISI, 'synapses:', N, 'took', time.time() - start_time, 's'
    left, right = time2index(interp_t, equilibrate-2.0, duration)
    return {ISI: rec[left:right]}


@interactive
def sim_stim_train_basal(x):
    """
    Called by controller, mapped to each engine. Simultaneously activates N synapses in the syn_list at 300 ms ISI and
    returns the resulting EPSP waveform.
    :param x: array of float
    :return: np.array
    """
    N = int(x[0] * 10000.)
    ISI = 300.
    duration = equilibrate + ISI * (num_stims - 1) + 101.
    syn_list.choose_syns_to_stim(N)
    N = len(syn_list.syns_to_stim)
    stim_time_array = [[] for i in range(N)]
    for stim_time in [equilibrate + ISI * i for i in range(num_stims)]:
        for j in local_random.sample(range(N), int(P0 * N)):
            stim_time_array[j].append(stim_time)
    stim_time_vector_array = [h.Vector(stim_time_list) for stim_time_list in stim_time_array]
    interp_t = np.arange(0., duration, interp_dt)
    sim.tstop = duration
    for i, syn in enumerate(syn_list.syns_to_stim):
        syn.source.play(stim_time_vector_array[i])
    start_time = time.time()
    sim.run(v_init)
    del stim_time_vector_array
    for syn in syn_list.syns_to_stim:
        syn.source.play(h.Vector())
    t = np.array(sim.tvec)
    left, right = time2index(t, equilibrate-2.0, equilibrate)
    vm = np.array(sim.rec_list[0]['vec'])
    baseline = np.average(vm[left:right])
    vm -= baseline
    rec = np.interp(interp_t, t, vm)
    print 'Process:', os.getpid(), 'ISI:', ISI, 'synapses:', N, 'took', time.time() - start_time, 's'
    left, right = time2index(interp_t, equilibrate-2.0, duration)
    return rec[left:right]


class Pr(object):
    """
    This object contains internal variables to track the evolution in time of parameters governing synaptic release
    probability, used during optimization, and then exported to pr.mod for use during patterned input simulations.
    """
    def __init__(self, P0, f, tau_F, d, tau_D):
        self.P0 = P0
        self.f = f
        self.tau_F = tau_F
        self.d = d
        self.tau_D = tau_D
        self.P = P0
        self.tlast = None
        self.F = 1.
        self.D = 1.

    def stim(self, stim_time):
        """
        Evolve the dynamics until the current stim_time, report the current P, and update the internal parameters.
        :param stim_time: float
        :return: float
        """
        if self.tlast is not None:
            self.F = 1. + (self.F - 1.) * np.exp(-(stim_time - self.tlast) / self.tau_F)
            self.D = 1. - (1. - self.D) * np.exp(-(stim_time - self.tlast) / self.tau_D)
            self.P = min(1., self.P0 * self.F * self.D)
        self.tlast = stim_time
        self.F += self.f
        self.D *= self.d
        return self.P


class SynList(object):
    """
    This object contains internal variables that will allow a subset of synapses to be chosen according to various
    criterion, and can easily be modified by calling internal methods within functions during optimization.
    """
    def __init__(self, cell, start=100., end=300., sec_type_list=['trunk', 'apical']):
        """

        :param cell: :class: 'SHocCell'
        :param start: float
        :param end: float
        :param sec_type_list: list of str
        """
        self.cell = cell
        self.internal_random = random.Random()
        self.internal_random.seed(0)
        self.syns_to_stim = []
        self.sec_type_list = sec_type_list
        self.all_syns_in_range = {sec_type: [] for sec_type in self.sec_type_list}
        for sec_type in self.sec_type_list:
            for node in self.cell.get_nodes_of_subtype(sec_type):
                if (start is None) or self.within_distance_range(node, start, end):
                    for spine in node.spines:
                        for syn in spine.synapses:
                            self.all_syns_in_range[sec_type].append(syn)
        self.num_syns_in_range = {sec_type: len(self.all_syns_in_range[sec_type]) for sec_type in self.sec_type_list}
        self.fraction_syns_in_range = {sec_type: float(self.num_syns_in_range[sec_type]) /
                                            float(np.sum(self.num_syns_in_range.values())) for sec_type in
                                            self.sec_type_list}

    def choose_syns_to_stim(self, N):
        """
        Based on the proportion of synapses in the dendritic sections pre-initialized in the SynList object, choose N
        random synapses.
        :param N: int
        """
        self.internal_random.seed(0)
        self.syns_to_stim = []
        for sec_type in self.sec_type_list:
            for syn in self.internal_random.sample(self.all_syns_in_range[sec_type],
                                                   int(N * self.fraction_syns_in_range[sec_type])):
                self.syns_to_stim.append(syn)

    def within_distance_range(self, node, start, end):
        """
        Check if dendrite origin of specified node is within the specified distance from the soma.
        :param node: :class:'SHocNode'
        :param start: float
        :param end: float
        :return: boolean
        """
        sec_type = node.type
        if sec_type == 'trunk':
            distance = self.cell.get_distance_to_node(self.cell.tree.root, node, loc=0.)
        else:
            distance = self.cell.get_distance_to_node(self.cell.tree.root, self.cell.get_dendrite_origin(node), loc=1.)
        if end is None:
            end = np.Infinity
        return start <= distance <= end


equilibrate = 250.  # time to steady-state
duration = 550.
v_init = -67.
syn_types = ['AMPA_KIN', 'NMDA_KIN5']

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.zero_na()

local_random = random.Random()
local_random.seed(os.getpid())
for branch in cell.trunk+cell.apical:
    for node in branch.spines:
        syn = Synapse(cell, node, syn_types, stochastic=0)  # stochasticity will be implemented in python rather than
                                                            # in neuron (hoc) during optimization
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
syn_list = SynList(cell)