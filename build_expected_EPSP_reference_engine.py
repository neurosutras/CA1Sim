__author__ = 'milsteina'
from specify_cells import *
import random
import os
import time
"""

"""
morph_filename = 'EB2-late-bifurcation.swc'
# mech_filename = '103115 interim dendritic excitability ampa nmda_kin3'
# mech_filename = '112915_less_excitable'
# mech_filename = '012316 alternate km kinetics'
# mech_filename = '012816 altered intrinsic properties - ampa nmda_kin4'
mech_filename = '020516 altered km2 rinp - ampa nmda_kin5'

rec_filename = 'expected_ref'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())


class EngineContainer(object):
    """
    This object contains internal variables that will allow the Controller to initialize each engine with a cell seed,
    and the desired number of excitatory and inhibitory synapses, by calling internal methods.
    """
    def __init__(self, cell):
        """

        :param cell: 'SHocCell'
        """
        self.cell = cell
        self.local_random = random.Random()
        self.all_exc_syns = {sec_type: [] for sec_type in ['basal', 'trunk', 'apical', 'tuft']}
        self.all_inh_syns = {sec_type: [] for sec_type in ['soma', 'basal', 'trunk', 'apical', 'tuft']}

        # place synapses in every spine
        for sec_type in self.all_exc_syns:
            for node in self.cell.get_nodes_of_subtype(sec_type):
                for spine in node.spines:
                    syn = Synapse(self.cell, spine, syn_types, stochastic=0)
                    self.all_exc_syns[sec_type].append(syn)
        self.cell.init_synaptic_mechanisms()

        # collate inhibitory synapses
        for sec_type in self.all_inh_syns:
            for node in self.cell.get_nodes_of_subtype(sec_type):
                for syn in node.synapses:
                    if 'GABA_A_KIN' in syn._syn:
                        self.all_inh_syns[sec_type].append(syn)
        self.stim_exc_syn_list = []

    def distribute_synapses(self, seed, num_exc_syns, num_inh_syns):
        """
        This code has been wrapped into an internal method in order to allow the Controller to choose the seed and
        number of synapses after initializing each engine.
        :param seed: int
        :param num_exc_syns: int
        :param num_inh_syns: int
        :return: boolean
        """
        self.local_random.seed(seed)
        # get the fraction of total spines contained in each sec_type
        total_exc_syns = {sec_type: len(self.all_exc_syns[sec_type]) for sec_type in ['basal', 'trunk', 'apical',
                                                                                      'tuft']}
        fraction_exc_syns = {sec_type: float(total_exc_syns[sec_type]) / float(np.sum(total_exc_syns.values())) for
                             sec_type in ['basal', 'trunk', 'apical', 'tuft']}
        stim_exc_syns = {'CA3': [], 'ECIII': []}
        stim_inh_syns = {'perisomatic': [], 'apical dendritic': [], 'distal apical dendritic': [],
                              'tuft feedforward': [], 'tuft feedback': []}
        peak_locs = {'CA3': [], 'ECIII': []}
        for sec_type in self.all_exc_syns:
            for i in self.local_random.sample(range(len(self.all_exc_syns[sec_type])),
                                              int(num_exc_syns*fraction_exc_syns[sec_type])):
                syn = self.all_exc_syns[sec_type][i]
                if sec_type == 'tuft':
                    stim_exc_syns['ECIII'].append(syn)
                else:
                    stim_exc_syns['CA3'].append(syn)

        # get the fraction of inhibitory synapses contained in each sec_type
        total_inh_syns = {sec_type: len(self.all_inh_syns[sec_type]) for sec_type in ['soma', 'basal', 'trunk',
                                                                                      'apical', 'tuft']}
        fraction_inh_syns = {sec_type: float(total_inh_syns[sec_type]) / float(np.sum(total_inh_syns.values())) for
                             sec_type in ['soma', 'basal', 'trunk', 'apical', 'tuft']}
        num_inh_syns = min(num_inh_syns, int(np.sum(total_inh_syns.values())))

        for sec_type in self.all_inh_syns:
            for i in self.local_random.sample(range(len(self.all_inh_syns[sec_type])),
                                              int(num_inh_syns*fraction_inh_syns[sec_type])):
                syn = self.all_inh_syns[sec_type][i]
                if syn.node.type == 'tuft':
                    if self.cell.is_terminal(syn.node):
                        # GABAergic synapses on terminal tuft branches are about 25% feedforward
                        group = self.local_random.choice(['tuft feedforward', 'tuft feedback', 'tuft feedback',
                                                          'tuft feedback'])
                    else:
                        # GABAergic synapses on intermediate tuft branches are about 50% feedforward
                        group = self.local_random.choice(['tuft feedforward', 'tuft feedback'])
                elif syn.node.type == 'trunk':
                    distance = self.cell.get_distance_to_node(self.cell.tree.root, syn.node, syn.loc)
                    if distance <= 50.:
                        group = 'perisomatic'
                    elif distance <= 150.:
                        group = 'apical dendritic'
                    else:
                        group = self.local_random.choice(['apical dendritic', 'distal apical dendritic',
                                                          'distal apical dendritic'])
                elif syn.node.type == 'basal':
                    distance = self.cell.get_distance_to_node(self.cell.tree.root, syn.node, syn.loc)
                    group = 'perisomatic' if distance <= 50. and not self.cell.is_terminal(syn.node) else \
                        'apical dendritic'
                elif syn.node.type == 'soma':
                    group = 'perisomatic'
                elif syn.node.type == 'apical':
                    distance = self.cell.get_distance_to_node(self.cell.tree.root,
                                                              self.cell.get_dendrite_origin(syn.node), loc=1.)
                    if distance <= 150.:
                        group = 'apical dendritic'
                    else:
                        group = self.local_random.choice(['apical dendritic', 'distal apical dendritic',
                                                          'distal apical dendritic'])
                stim_inh_syns[group].append(syn)

        gauss_sigma = global_theta_cycle_duration * input_field_width / 3. / np.sqrt(2.)  # contains 99.7% gaussian area

        for group in stim_exc_syns.keys():
            if stim_exc_syns[group]:
                peak_locs[group] = np.arange(-0.75 * input_field_duration, (0.75 + track_length) * input_field_duration,
                                  (1.5 + track_length) * input_field_duration / int(len(stim_exc_syns[group])))
                peak_locs[group] = peak_locs[group][:len(stim_exc_syns[group])]
            self.local_random.shuffle(peak_locs[group])
            peak_locs[group] = list(peak_locs[group])

        for group in stim_exc_syns.keys():
            for syn in stim_exc_syns[group]:
                self.stim_exc_syn_list.append(syn)

        # modulate the weights of inputs that have peak_locs along this stretch of the track
        modulated_field_center = track_duration * 0.6
        gauss_mod_amp = {}

        for group in stim_exc_syns.keys():
            gauss_mod_amp[group] = 1.5 * np.exp(-((np.array(peak_locs[group]) - modulated_field_center) /
                                                  (gauss_sigma * 1.4)) ** 2.) + 1.
            for i, syn in enumerate(stim_exc_syns[group]):
                syn.netcon('AMPA_KIN').weight[0] = gauss_mod_amp[group][i]
        return True


def stim_single_exc_syn(index):
    """

    :param index: int
    """
    syn = local_container.stim_exc_syn_list[index]
    node_index = syn.node.index
    syn.source.play(h.Vector([equilibrate]))
    start_time = time.time()
    sim.run(v_init)
    syn.source.play(h.Vector())
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.parameters['spine_index'] = node_index
        sim.export_to_file(f, node_index)
    print 'Process: %i took %i s to stimulate synapse with index %i' % (os.getpid(), time.time() - start_time,
                                                                         node_index)
    return rec_filename


NMDA_type = 'NMDA_KIN5'

equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # (ms)
input_field_width = 20  # (theta cycles per 6 standard deviations)
excitatory_phase_extent = 450.  # (degrees)
# Geissler...Buzsaki, PNAS 2010
unit_theta_cycle_duration = global_theta_cycle_duration * input_field_width / (input_field_width +
                                                                               (excitatory_phase_extent / 360.))
input_field_duration = input_field_width * global_theta_cycle_duration
track_length = 2.5  # field widths
track_duration = track_length * input_field_duration
duration = equilibrate + 200.

v_init = -67.

syn_types = ['AMPA_KIN', NMDA_type]

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.set_terminal_branch_nas_gradient()
cell.insert_inhibitory_synapses_in_subset()

local_container = EngineContainer(cell)

trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
if trunk_bifurcation:
    trunk_branches = [branch for branch in trunk_bifurcation[0].children if branch.type == 'trunk']
    # get where the thickest trunk branch gives rise to the tuft
    trunk = max(trunk_branches, key=lambda node: node.sec(0.).diam)
    trunk = (node for node in cell.trunk if cell.node_in_subtree(trunk, node) and
             'tuft' in (child.type for child in node.children)).next()
else:
    trunk_bifurcation = [node for node in cell.trunk if 'tuft' in (child.type for child in node.children)]
    trunk = trunk_bifurcation[0]

sim = QuickSim(duration, verbose=0)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['global_theta_cycle_duration'] = global_theta_cycle_duration
sim.parameters['input_field_duration'] = input_field_duration
sim.parameters['track_length'] = track_length
sim.parameters['duration'] = duration
sim.append_rec(cell, cell.tree.root, description='soma', loc=0.5)
sim.append_rec(cell, trunk, description='distal_trunk', loc=0.)
sim.append_rec(cell, trunk_bifurcation[0], description='proximal_trunk', loc=1.)