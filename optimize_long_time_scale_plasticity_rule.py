__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import random
import sys

"""
In this version of the simulation, phase precession of CA3 inputs is implemented using the method from Chadwick et al.,
Elife, 2015, which uses a circular gaussian with a phase sensitivity factor that effectively compresses the range of
phases within each theta cycle that each input is active, which will reduce jitter across within-cycle input sequences.

"""
morph_filename = 'EB2-late-bifurcation.swc'

mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'

if len(sys.argv) > 1:
    synapses_seed = int(sys.argv[1])
else:
    synapses_seed = 0
if len(sys.argv) > 2:
    num_exc_syns = int(sys.argv[2])
else:
    num_exc_syns = 1600
if len(sys.argv) > 3:
    num_inh_syns = int(sys.argv[3])
else:
    num_inh_syns = 600
# allows parallel computation of multiple trials for the same spines with the same peak_locs, but with different
# input spike trains and stochastic synapses for each trial
if len(sys.argv) > 4:
    trial_seed = int(sys.argv[4])
else:
    trial_seed = 0
# location of first induced field, as a fraction of the track
if len(sys.argv) > 5:
    field1_loc = float(sys.argv[5])
else:
    field1_loc = 0.5
# constant run velocity
if len(sys.argv) > 6:
    run_vel = float(sys.argv[6])
else:
    run_vel = 30.  # cm/s


class History(object):
    def __init__(self):
        """

        """
        self.x = []
        self.Err = []

    def report_best(self):
        """

        :return: array
        """
        lowest_Err = min(self.Err)
        index = self.Err.index(lowest_Err)
        best_x = self.x[index]
        formatted_x = '[' + ', '.join(['%.2E' % xi for xi in best_x]) + ']'
        print 'best x: %s' % formatted_x
        print 'lowest Err: %.3E' % lowest_Err
        return best_x

    def export(self, hist_filename):
        """
        Save the history to .pkl
        """
        saved_history = {'x': self.x, 'Err': self.Err}
        write_to_pkl(data_dir+hist_filename+'.pkl', saved_history)


hist = History()


def get_dynamic_theta_phase_force(phase_ranges, peak_loc, input_field_duration, stim_t, dt):
    """
    Expects a list of tuples containing times and phases relative to peak_loc and the non-modulated phase preference
    (zero degrees). Returns a waveform of phase vs time.
    :param phase_ranges: list of tuple (ms, degrees)
    :param peak_loc:
    :param input_field_duration:
    :param stim_t:
    :param dt:
    :return: :class: 'np.array'
    """
    start_phase_val = phase_ranges[0][1] * 2. * np.pi / 360.  # convert degrees to radians
    end_phase_val = phase_ranges[-1][1] * 2. * np.pi / 360.  # convert degrees to radians
    phase_force = np.ones_like(stim_t) * start_phase_val
    phase_gradient = np.array([])
    for i in range(len(phase_ranges) - 1):
        t0 = phase_ranges[i][0]
        t1 = phase_ranges[i + 1][0]
        phase0 = phase_ranges[i][1] * 2. * np.pi / 360.  # convert degrees to radians
        phase1 = phase_ranges[i + 1][1] * 2. * np.pi / 360.
        del_t = t1 - t0
        del_phase = phase1 - phase0
        if abs(del_phase) > 0.:
            del_phase = del_phase / del_t * dt
            this_range_piece = np.arange(phase0, phase1, del_phase)
        else:
            this_range_piece = np.ones(int(del_t / dt)) * phase0
        phase_gradient = np.append(phase_gradient, this_range_piece)
    if stim_t[0] <= peak_loc - input_field_duration * 0.5 <= stim_t[-1]:
        phase_start = np.where(peak_loc - input_field_duration * 0.5 >= stim_t)[0]
        if np.any(phase_start):
            phase_start = phase_start[-1]
            phase_end = min(len(stim_t), phase_start + len(phase_gradient))
            phase_force[:phase_start] = start_phase_val
            phase_force[phase_start:phase_end] = phase_gradient[:phase_end - phase_start]
            phase_force[phase_end:] = end_phase_val
    elif stim_t[0] <= peak_loc + input_field_duration * 0.5 <= stim_t[-1]:
        phase_end = np.where(peak_loc + input_field_duration * 0.5 >= stim_t)[0]
        if np.any(phase_end):
            phase_end = phase_end[-1]
            phase_start = max(0, phase_end - len(phase_gradient))
            phase_force[:phase_start] = start_phase_val
            phase_force[phase_start:phase_end] = phase_gradient[-(phase_end - phase_start):]
            phase_force[phase_end:] = end_phase_val
    return phase_force


def run_trial(simiter, global_phase_offset=0.):
    """

    :param simiter: int
    :param global_phase_offset: float
    """
    local_random.seed(simiter)
    index = 0
    for group in stim_exc_syns:
        exc_stim_forces[group] = []
        for i, syn in enumerate(stim_exc_syns[group]):
            peak_loc = peak_locs[group][i]
            peak_time = peak_loc * dt / dx
            gauss_force = excitatory_peak_rate[group] * np.exp(-((stim_x - peak_loc) / gauss_sigma)**2.)
            before = np.array(gauss_force[:int(track_duration/dt)])
            after = np.array(gauss_force[int((2. * track_duration - track_equilibrate)/dt):int(3. * track_duration /
                                                                                               dt)])
            gauss_force[int((track_duration - track_equilibrate)/dt):int(2. * track_duration / dt)] += after
            gauss_force[int(track_duration / dt):int(2. * track_duration / dt)] += before
            gauss_force = np.array(gauss_force[int((track_duration - track_equilibrate)/dt):int(2. * track_duration /
                                                                                                dt)])
            if group in excitatory_precession_range:
                phase_force = get_dynamic_theta_phase_force(excitatory_precession_range[group], peak_time,
                                                            input_field_duration, extended_stim_t, dt)
                if peak_time - 1.2 * input_field_duration / 2. < 0.:
                    before = np.array(phase_force[int((track_duration - 1.2 * input_field_duration / 2. +
                                                       peak_time)/dt):int(track_duration / dt)])
                    phase_force[int((2. * track_duration - 1.2 * input_field_duration / 2. + peak_time)/dt):
                                int(2. * track_duration / dt)] = before
                elif peak_time + 1.2 * input_field_duration / 2. > track_duration:
                    after = np.array(phase_force[int((2. * track_duration - track_equilibrate)/dt):
                                                int((track_duration + 1.2 * input_field_duration / 2. + peak_time)/dt)])
                    phase_force[int((track_duration - track_equilibrate)/dt):
                                int((1.2 * input_field_duration / 2. + peak_time) / dt)] = after
                phase_force = np.array(phase_force[int((track_duration - track_equilibrate)/dt):
                                                    int(2. * track_duration / dt)])
                theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] * np.cos(phase_force +
                                        excitatory_theta_phase_offset[group] - 2. * np.pi * stim_t /
                                        global_theta_cycle_duration + global_phase_offset))
            else:
                theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] *
                                 np.cos(excitatory_theta_phase_offset[group] - 2. * np.pi * stim_t /
                                        global_theta_cycle_duration + global_phase_offset))
            theta_force -= np.min(theta_force)
            theta_force /= np.max(theta_force)
            theta_force *= excitatory_theta_modulation_depth[group]
            theta_force += 1. - excitatory_theta_modulation_depth[group]
            stim_force = np.multiply(gauss_force, theta_force)
            exc_stim_forces[group].append(stim_force)
            index += 1


NMDA_type = 'NMDA_KIN5'

dt = 0.02  # ms
dx = run_vel / 1000. * dt  # cm
equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = -0.5 * run_vel + 155. # 150.  # (ms)
input_field_width = 90.  # cm
input_field_duration = input_field_width * 1000. / run_vel  # ms
track_length = 180.  # cm
track_duration = track_length * 1000. / run_vel  # ms
track_equilibrate = 2. * global_theta_cycle_duration
duration = equilibrate + track_equilibrate + track_duration  # input_field_duration
excitatory_peak_rate = {'CA3': 40., 'ECIII': 40.}
excitatory_theta_modulation_depth = {'CA3': 0.7, 'ECIII': 0.7}
# From Chadwick et al., ELife 2015
excitatory_theta_phase_tuning_factor = {'CA3': 0.8, 'ECIII': 0.8}
excitatory_precession_range = {}
excitatory_precession_range['CA3'] = [(-input_field_duration*0.5, 180.), (-input_field_duration*0.35, 180.),
                                      (input_field_duration*0.35, -180.), (input_field_duration*0.5, -180.)]  # (ms, degrees)
excitatory_theta_phase_offset = {}
excitatory_theta_phase_offset['CA3'] = 165. / 360. * 2. * np.pi  # radians
excitatory_theta_phase_offset['ECIII'] = 0. / 360. * 2. * np.pi  # radians
excitatory_stochastic = 1
inhibitory_mean_rate = {'perisomatic': 25., 'axo-axonic': 25., 'apical dendritic': 25., 'distal apical dendritic': 25.,
                        'tuft feedforward': 25., 'tuft feedback': 25.}
inhibitory_theta_modulation_depth = {'perisomatic': 0.5, 'axo-axonic': 0.5, 'apical dendritic': 0.5,
                                     'distal apical dendritic': 0.5, 'tuft feedforward': 0.5, 'tuft feedback': 0.5}
inhibitory_theta_phase_tuning_factor = {'perisomatic': 0.6, 'axo-axonic': 0.6, 'apical dendritic': 0.6,
                                     'distal apical dendritic': 0.6, 'tuft feedforward': 0.6, 'tuft feedback': 0.6}
inhibitory_precession_range = {}
inhibitory_theta_phase_offset = {}
inhibitory_theta_phase_offset['perisomatic'] = 135. / 360. * 2. * np.pi  # Like PV+ Basket
inhibitory_theta_phase_offset['axo-axonic'] = 45. / 360. * 2. * np.pi  # Vargas et al., ELife, 2014
inhibitory_theta_phase_offset['apical dendritic'] = 200. / 360. * 2. * np.pi  # Like PYR-layer Bistratified
inhibitory_theta_phase_offset['distal apical dendritic'] = 180. / 360. * 2. * np.pi  # Like SR/SLM Border Cells
inhibitory_theta_phase_offset['tuft feedforward'] = 340. / 360. * 2. * np.pi  # Like Neurogliaform
inhibitory_theta_phase_offset['tuft feedback'] = 200. / 360. * 2. * np.pi  # Like SST+ O-LM

v_init = -67.

syn_types = ['AMPA_KIN', NMDA_type]

local_random = random.Random()

# choose a subset of synapses to stimulate with inhomogeneous poisson rates
local_random.seed(synapses_seed)

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.set_terminal_branch_na_gradient()
cell.insert_inhibitory_synapses_in_subset()

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

all_exc_syns = {sec_type: [] for sec_type in ['basal', 'trunk', 'apical', 'tuft']}
all_inh_syns = {sec_type: [] for sec_type in ['soma', 'ais', 'basal', 'trunk', 'apical', 'tuft']}
stim_exc_syns = {'CA3': [], 'ECIII': []}
stim_inh_syns = {'perisomatic': [], 'axo-axonic': [], 'apical dendritic': [], 'distal apical dendritic': [],
                 'tuft feedforward': [], 'tuft feedback': []}
stim_successes = []
peak_locs = {'CA3': [], 'ECIII': []}

# place synapses in trunk for inheritance of mechanisms (for testing)
if 'trunk' not in all_exc_syns:
    for node in cell.trunk:
        for spine in node.spines:
            syn = Synapse(cell, spine, syn_types, stochastic=excitatory_stochastic)

# place synapses in every spine
for sec_type in all_exc_syns:
    for node in cell.get_nodes_of_subtype(sec_type):
        for spine in node.spines:
            syn = Synapse(cell, spine, syn_types, stochastic=excitatory_stochastic)
            all_exc_syns[sec_type].append(syn)
cell.init_synaptic_mechanisms()

# collate inhibitory synapses
for sec_type in all_inh_syns:
    for node in cell.get_nodes_of_subtype(sec_type):
        for syn in node.synapses:
            if 'GABA_A_KIN' in syn._syn:
                all_inh_syns[sec_type].append(syn)

sim = QuickSim(duration, cvode=0, dt=0.01)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['track_equilibrate'] = track_equilibrate
sim.parameters['global_theta_cycle_duration'] = global_theta_cycle_duration
sim.parameters['input_field_duration'] = input_field_duration
sim.parameters['track_length'] = track_length
sim.parameters['duration'] = duration
sim.parameters['stim_dt'] = dt
sim.append_rec(cell, cell.tree.root, description='soma', loc=0.)
sim.append_rec(cell, trunk_bifurcation[0], description='proximal_trunk', loc=1.)
sim.append_rec(cell, trunk, description='distal_trunk', loc=1.)
spike_output_vec = h.Vector()
cell.spike_detector.record(spike_output_vec)

# get the fraction of total spines contained in each sec_type
total_exc_syns = {sec_type: len(all_exc_syns[sec_type]) for sec_type in ['basal', 'trunk', 'apical', 'tuft']}
fraction_exc_syns = {sec_type: float(total_exc_syns[sec_type]) / float(np.sum(total_exc_syns.values())) for sec_type in
                 ['basal', 'trunk', 'apical', 'tuft']}

for sec_type in all_exc_syns:
    for i in local_random.sample(range(len(all_exc_syns[sec_type])), int(num_exc_syns*fraction_exc_syns[sec_type])):
        syn = all_exc_syns[sec_type][i]
        if sec_type == 'tuft':
            stim_exc_syns['ECIII'].append(syn)
        else:
            stim_exc_syns['CA3'].append(syn)

# get the fraction of inhibitory synapses contained in each sec_type
total_inh_syns = {sec_type: len(all_inh_syns[sec_type]) for sec_type in ['soma', 'ais', 'basal', 'trunk', 'apical',
                                                                         'tuft']}
fraction_inh_syns = {sec_type: float(total_inh_syns[sec_type]) / float(np.sum(total_inh_syns.values())) for sec_type in
                 ['soma', 'ais', 'basal', 'trunk', 'apical', 'tuft']}
num_inh_syns = min(num_inh_syns, int(np.sum(total_inh_syns.values())))

for sec_type in all_inh_syns:
    for i in local_random.sample(range(len(all_inh_syns[sec_type])), int(num_inh_syns*fraction_inh_syns[sec_type])):
        syn = all_inh_syns[sec_type][i]
        if syn.node.type == 'tuft':
            if cell.is_terminal(syn.node):
                # GABAergic synapses on terminal tuft branches are about 25% feedforward
                group = local_random.choice(['tuft feedforward', 'tuft feedback', 'tuft feedback', 'tuft feedback'])
            else:
                # GABAergic synapses on intermediate tuft branches are about 50% feedforward
                group = local_random.choice(['tuft feedforward', 'tuft feedback'])
        elif syn.node.type == 'trunk':
            distance = cell.get_distance_to_node(cell.tree.root, syn.node, syn.loc)
            if distance <= 50.:
                group = 'perisomatic'
            elif distance <= 150.:
                group = local_random.choice(['apical dendritic', 'apical dendritic', 'distal apical dendritic'])
            else:
                group = local_random.choice(['apical dendritic', 'distal apical dendritic', 'distal apical dendritic'])
        elif syn.node.type == 'basal':
            distance = cell.get_distance_to_node(cell.tree.root, syn.node, syn.loc)
            group = 'perisomatic' if distance <= 50. and not cell.is_terminal(syn.node) else 'apical dendritic'
        elif syn.node.type == 'soma':
            group = 'perisomatic'
        elif syn.node.type == 'apical':
            distance = cell.get_distance_to_node(cell.tree.root, cell.get_dendrite_origin(syn.node), loc=1.)
            if distance <= 150.:
                group = local_random.choice(['apical dendritic', 'apical dendritic', 'distal apical dendritic'])
            else:
                group = local_random.choice(['apical dendritic', 'distal apical dendritic', 'distal apical dendritic'])
        elif syn.node.type == 'ais':
            group = 'axo-axonic'
        stim_inh_syns[group].append(syn)

stim_x = np.arange(-track_length, 2. * track_length, dx)
stim_t = np.arange(-track_equilibrate, track_duration, dt)
extended_stim_t = np.arange(-track_duration, 2. * track_duration, dt)

gauss_sigma = input_field_width / 3. / np.sqrt(2.)  # contains 99.7% gaussian area

rand_exc_seq_locs = {}
delta_peak_locs = {}
extended_peak_locs = {}
within_track = {}
for group in stim_exc_syns:
    rand_exc_seq_locs[group] = []
    delta_peak_locs[group] = track_length / len(stim_exc_syns[group])
    if stim_exc_syns[group]:
        extended_peak_locs[group] = np.arange(-track_length, 2. * track_length, delta_peak_locs[group])
        within_track[group] = np.where((extended_peak_locs[group] >= 0.) &
                                       (extended_peak_locs[group] <= track_length))[0]
        peak_locs[group] = extended_peak_locs[group][within_track[group]]

for group in stim_exc_syns:
    for syn in stim_exc_syns[group]:
        all_exc_syns[syn.node.parent.parent.type].remove(syn)

# modulate the weights of inputs with peak_locs along this stretch of the track
field_center1 = track_length * field1_loc
extended_cos_mod_weight = {}
cos_mod_weight = {}
peak_mod_weight = 2.5
tuning_amp = (peak_mod_weight - 1.)
start_loc = field_center1 - input_field_width * 1.2 / 2.
end_loc = field_center1 + input_field_width * 1.2 / 2.

for group in stim_exc_syns:
    extended_cos_mod_weight[group] = np.zeros_like(extended_peak_locs[group])
    within_field = np.where((extended_peak_locs[group] >= start_loc) & (extended_peak_locs[group] <= end_loc))[0]
    extended_cos_mod_weight[group][within_field] = tuning_amp / 2. * np.cos(2. * np.pi / (input_field_width * 1.2) *
                                                                            (extended_peak_locs[group][within_field] -
                                                                             field_center1)) + tuning_amp / 2.
    cos_mod_weight[group] = extended_cos_mod_weight[group][within_track[group]]
    before_track_indexes = np.where(extended_peak_locs[group] < 0.)[0]
    before_track = np.array(extended_cos_mod_weight[group][before_track_indexes[-min(len(before_track_indexes),
                                                                                     len(cos_mod_weight[group])):]])
    after_track_indexes = np.where(extended_peak_locs[group] > track_length)[0]
    after_track = np.array(extended_cos_mod_weight[group][after_track_indexes[:min(len(after_track_indexes),
                                                                                     len(cos_mod_weight[group])):]])
    cos_mod_weight[group][-len(before_track):] += before_track
    cos_mod_weight[group][:len(after_track)] += after_track
    cos_mod_weight[group] += 1.
exc_stim_forces = {}
inh_stim_forces = {}
global_phase_offset = 0.  # -2.*np.pi*165./360.
run_trial(trial_seed, global_phase_offset)

"""
weighted_forces_orig = {}
weighted_forces_orig_sum = np.zeros_like(stim_t)
for group in ['CA3', 'ECIII']:
    if group not in weighted_forces_orig:
        weighted_forces_orig[group] = []
    for i, stim_force in enumerate(exc_stim_forces[group]):
        weighted_forces_orig[group].append(stim_force * cos_mod_weight[group][i] / 1000.)
    weighted_forces_orig_sum = np.add(weighted_forces_orig_sum, np.sum(weighted_forces_orig[group], axis=0))
"""

ref_theta = 0.5 * np.cos(2. * np.pi * stim_t / global_theta_cycle_duration - global_phase_offset) + 1.75

rule_phase_offset = 45. / 360. * 2. * np.pi
induction_duration = 300.


def check_bounds(x, xmin, xmax):
    """
    Check that the current set of parameters are within the specified bounds.
    :param x: array
    :param xmin: array
    :param xmax: array
    :return: bool
    """
    for i in range(len(x)):
        if ((xmin[i] is not None and x[i] < xmin[i]) or
                (xmax[i] is not None and x[i] > xmax[i])):
            return False
    return True


def build_rule_waveform0(x, field_center, plot=False):
    """
    Construct a normalized rule waveform with the following components:
    0) phase-dependent rule during current injection
    :param x: array: [during_depth]
    :param field_center: float
    :param plot: bool
    :return: array
    """
    during_depth = x[0]
    pot_amp = (1. - during_depth) / 2.
    induction_loc = field_center

    induction_start = induction_loc * 1000. / run_vel
    rule_theta = pot_amp * np.cos(2. * np.pi * extended_stim_t / global_theta_cycle_duration -
                                               global_phase_offset - rule_phase_offset) + pot_amp + during_depth
    rule_waveform = np.zeros_like(extended_stim_t)
    left = int((track_duration + induction_start) / dt)
    right = left + int(induction_duration / dt)
    rule_waveform[left:right] = np.array(rule_theta[left:right])
    before_rule = np.array(rule_waveform[:int(track_duration / dt)])
    after_rule = np.array(rule_waveform[int(2. * track_duration / dt):])
    within_rule = np.array(rule_waveform[int(track_duration / dt):int(2. * track_duration / dt)])
    rule_waveform = np.zeros_like(stim_t)
    rule_waveform[int(track_equilibrate / dt):int((track_equilibrate + track_duration) / dt)] = before_rule + \
                                                                                                after_rule + \
                                                                                                within_rule
    if plot:
        left = int((track_equilibrate + induction_start) / dt)
        right = left + int(induction_duration / dt)
        fig, axes = plt.subplots(1)
        axes.plot(stim_t[left:right] - induction_start, rule_waveform[left:right], label='Plasticity rule')
        axes.plot(stim_t[left:right] - induction_start, ref_theta[left:right], label='LFP theta')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Time (ms)')
        axes.set_ylabel('Normalized change in synaptic weight')
        axes.set_ylim([0., 2.5])
        axes.set_xlim([0., induction_duration])
        clean_axes(axes)
        plt.show()
        plt.close()

    return rule_waveform


def build_rule_waveform1(x, field_center, plot=False):
    """
    Construct a normalized rule waveform with the following components:
    0) phase-dependent rule during current injection
    1) time-dependent potentiation before current injection
    2) time-dependent potentiation after current injection
    :param x: array: [during_depth, pre_pot_dur, post_pot_dur]
    :param field_center: float
    :param plot: bool
    :return: array
    """
    during_depth = x[0]
    pre_pot_dur = x[1]
    post_pot_dur = x[2]
    pot_amp = (1. - during_depth) / 2.
    induction_loc = field_center

    induction_start = induction_loc * 1000. / run_vel
    rule_theta = pot_amp * np.cos(2. * np.pi * extended_stim_t / global_theta_cycle_duration -
                                               global_phase_offset - rule_phase_offset) + pot_amp + during_depth
    rule_waveform = np.zeros_like(extended_stim_t)

    left = int((track_duration + induction_start - pre_pot_dur) / dt)
    right = left + int(pre_pot_dur / dt)
    pre_pot_phase_offset = 2. * np.pi * ((induction_start - pre_pot_dur) % (pre_pot_dur * 2.)) / (pre_pot_dur * 2.)
    pre_pot_waveform = during_depth * (0.5 * np.cos(2. * np.pi * extended_stim_t / (pre_pot_dur * 2.) - np.pi -
                                                    pre_pot_phase_offset) + 0.5)
    rule_waveform[left:right] = np.array(pre_pot_waveform[left:right])
    left = right
    right = left + int(induction_duration / dt)
    rule_waveform[left:right] = np.array(rule_theta[left:right])
    left = right
    right = left + int(post_pot_dur / dt)
    post_pot_phase_offset = 2. * np.pi * ((induction_start + induction_duration) % (post_pot_dur * 2.)) / \
                            (post_pot_dur * 2.)
    post_waveform1 = during_depth * (0.5 * np.cos(2. * np.pi * extended_stim_t / (post_pot_dur * 2.) -
                                                  post_pot_phase_offset) + 0.5)
    rule_waveform[left:right] = np.array(post_waveform1[left:right])
    before_rule = np.array(rule_waveform[:int(track_duration / dt)])
    after_rule = np.array(rule_waveform[int(2. * track_duration / dt):])
    within_rule = np.array(rule_waveform[int(track_duration / dt):int(2. * track_duration / dt)])
    rule_waveform = np.zeros_like(stim_t)
    rule_waveform[int(track_equilibrate / dt):int((track_equilibrate + track_duration) / dt)] = before_rule + \
                                                                                                after_rule + \
                                                                                                within_rule
    if plot:
        left = int((track_equilibrate + induction_start - pre_pot_dur) / dt)
        right = left + int((induction_duration + pre_pot_dur + post_pot_dur) / dt)
        fig, axes = plt.subplots(1)
        axes.plot(stim_t[left:right] - induction_start, rule_waveform[left:right], label='Plasticity rule')
        axes.plot(stim_t[left:right] - induction_start, ref_theta[left:right], label='LFP theta')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Time (ms)')
        axes.set_ylabel('Normalized change in synaptic weight')
        axes.set_ylim([0., 2.5])
        axes.set_xlim([-pre_pot_dur, induction_duration+post_pot_dur])
        clean_axes(axes)
        plt.show()
        plt.close()

    return rule_waveform


def calculate_delta_weights(rule_waveform, target_delta_weights, plot=False):
    """
    Given any rule waveform and target_delta_weights, compute the resulting delta_weights (without saturation).
    :param rule_waveform:
    :param target_delta_weights:
    :param plot: bool
    :return: delta_weights: array, rule_gain: float
    """
    delta_weights = np.zeros_like(target_delta_weights)
    group = 'CA3'
    for i, stim_force in enumerate(exc_stim_forces[group]):
        this_force = 0.001 * dt * np.multiply(stim_force, rule_waveform)
        this_area = np.trapz(this_force, dx=dt)
        delta_weights[i] = this_area
    rule_gain = 1.5 / np.max(delta_weights)
    delta_weights *= rule_gain
    if plot:
        fig, axes = plt.subplots(1)
        axes.scatter(peak_locs['CA3'], target_delta_weights+1., color='k', label='Target weights')
        axes.scatter(peak_locs['CA3'], delta_weights + 1., color='r', label='Weights from plasticity rule')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Normalized synaptic weight')
        axes.set_ylim([0.9, 3.5])
        axes.set_xlim([0., track_length])
        clean_axes(axes)
        plt.show()
        plt.close()
    return delta_weights, rule_gain


def get_mod_weights(x, field_center, existing_weights, rule_gain):
    """

    :param x: [pre_induction_pot_dur, post_induction_pot_dur, post_induction_depot_dur, pot_DC, depot_trough]
    :param field_center: float  # cm
    :param existing_weights: array
    :return: spdp_mod_weight, rule_waveform, new_rule_gain: array, array, float
    """
    pre_induction_pot_dur = x[0]
    post_induction_pot_dur = x[1]
    post_induction_depot_dur = x[2]
    pot_DC = x[3]
    depot_trough = x[4]
    pot_amp = (1. - pot_DC) / 2.
    induction_loc = field_center

    induction_start = induction_loc * 1000. / run_vel
    rule_theta = rule_gain * (pot_amp * np.cos(2. * np.pi * extended_stim_t / global_theta_cycle_duration -
                                               global_phase_offset - rule_phase_offset) + pot_amp + pot_DC)
    rule_waveform = np.zeros_like(extended_stim_t)
    left = int((track_duration + induction_start - pre_induction_pot_dur) / dt)
    right = left + int(pre_induction_pot_dur / dt)
    pre_phase_offset = 2. * np.pi * ((induction_start - pre_induction_pot_dur) % (pre_induction_pot_dur * 2.)) / \
                       (pre_induction_pot_dur * 2.)
    pre_waveform = rule_gain * pot_DC * (0.5 * np.cos(2. * np.pi * extended_stim_t / (pre_induction_pot_dur * 2.) -
                                                      np.pi - pre_phase_offset) + 0.5)
    rule_waveform[left:right] = np.array(pre_waveform[left:right])
    left = right
    right = left + int(induction_duration / dt)
    rule_waveform[left:right] = np.array(rule_theta[left:right])
    left = right
    right = left + int(post_induction_pot_dur / dt)
    post_phase_offset1 = 2. * np.pi * ((induction_start + induction_duration) % (post_induction_pot_dur * 2.)) / \
                         (post_induction_pot_dur * 2.)
    post_waveform1 = rule_gain * pot_DC * (0.5 * np.cos(2. * np.pi * extended_stim_t /
                                                        (post_induction_pot_dur * 2.) - post_phase_offset1) + 0.5)
    rule_waveform[left:right] = np.array(post_waveform1[left:right])
    left = right
    right = left + int(post_induction_depot_dur / dt)
    post_phase_offset2 = 2. * np.pi * ((induction_start + induction_duration + post_induction_pot_dur) %
                                       post_induction_depot_dur) / post_induction_depot_dur
    post_waveform2 = rule_gain * (depot_trough * 0.5 * np.cos(2. * np.pi * extended_stim_t /
                                                              post_induction_depot_dur - post_phase_offset2) -
                                  depot_trough * 0.5)
    rule_waveform[left:right] = np.array(post_waveform2[left:right])
    before_rule = np.array(rule_waveform[:int(track_duration / dt)])
    after_rule = np.array(rule_waveform[int(2. * track_duration / dt):])
    within_rule = np.array(rule_waveform[int(track_duration / dt):int(2. * track_duration / dt)])
    rule_waveform = np.zeros_like(stim_t)
    rule_waveform[int(track_equilibrate / dt):int((track_equilibrate + track_duration) / dt)] = before_rule + \
                                                                                                after_rule + \
                                                                                                within_rule
    spdp_mod_weight = copy.deepcopy(existing_weights)
    for group in ['CA3', 'ECIII']:
        for i, stim_force in enumerate(exc_stim_forces[group]):
            this_force = 0.001 * dt * np.multiply(stim_force, rule_waveform)
            this_area = np.trapz(this_force, dx=dt)
            this_weight = this_area + spdp_mod_weight[group][i]
            if this_weight < 1.:
                this_weight = 1.
            elif this_weight > 3.:
                this_weight = 3.
            spdp_mod_weight[group][i] = this_weight
    rule_gain *= 1.5 / (np.max(spdp_mod_weight['CA3']) - 1.)
    return spdp_mod_weight, rule_waveform, rule_gain


def delta_weights_error(x, xmin, xmax, rule_function, target_delta_weights, field_center, plot=False):
    """

    :param x: array [during_depth, pre_pot_dur, post_pot_dur, depot_dur, depot_depth]
    :param xmin: array
    :param xmax: array
    :param rule_function: callable
    :param target_delta_weights: array
    :param field_center: float
    :param plot: bool
    :return: float
    """
    if not check_bounds(x, xmin, xmax):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    rule_waveform = rule_function(x, field_center, plot)
    delta_weights, rule_gain = calculate_delta_weights(rule_waveform, target_delta_weights, plot)

    Err = 0.
    for i in range(len(target_delta_weights)):
        Err += abs(target_delta_weights[i] - delta_weights[i])
    Err **= 2.

    hist.x.append(x)
    hist.Err.append(Err)
    formatted_x = '[' + ', '.join(['%.2f' % xi for xi in x]) + ']'

    print 'x:', formatted_x, 'Err:', Err, 'Rule gain:', rule_gain
    return Err


def optimize_polish(x, xmin, xmax, error_function, rule_function, target_delta_weights, field_center, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param rule_function: callable
    :param target_delta_weights: array
    :param field_center: float
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 200

    result = optimize.minimize(error_function, x, method='Nelder-Mead', options={'ftol': 1e-3,
                                                    'xtol': 1e-3, 'disp': True, 'maxiter': maxfev},
                               args=(xmin, xmax, rule_function, target_delta_weights, field_center))
    formatted_x = '['+', '.join(['%.2E' % xi for xi in result.x])+']'
    print 'Process: %i completed optimize_polish after %i iterations with Error: %.4E and x: %s' % (os.getpid(),
                                                                            result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


def optimize_explore(x, xmin, xmax, error_function, rule_function, target_delta_weights, field_center, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param rule_function: callable
    :param target_delta_weights: array
    :param field_center: float
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 400

    take_step = Normalized_Step(x, xmin, xmax)
    minimizer_kwargs = dict(method=null_minimizer, args=(xmin, xmax, rule_function, target_delta_weights, field_center))
    result = optimize.basinhopping(error_function, x, niter=maxfev, niter_success=maxfev/2,
                                       disp=True, interval=20, minimizer_kwargs=minimizer_kwargs, take_step=take_step)
    formatted_x = '['+', '.join(['%.2E' % xi for xi in result.x])+']'
    print 'Process: %i completed optimize_explore after %i iterations with Error: %.4E and x: %s' % (os.getpid(),
                                                                            result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


# [pre_induction_pot_dur, post_induction_pot_dur, post_induction_depot_dur, pot_DC, depot_trough]
x0 = [1505.6, 1124.6, 564.5, 0.70, 1.94]  # rule gain: 4.15
xmin = [750., 750., 500., 0.2, 0.11]
xmax = [2250., 2250., 1500., 1., 3.]

# [during_depth, pre_pot_dur, post_pot_dur, depot_dur, depot_depth]

# [during_depth]
# x0 = [0.7]
x0 = [0.0594]
xmin0 = [0.]
xmax0 = [1.]

# [during_depth, pre_pot_dur, post_pot_dur]
# x1 = [0.7, 1505.6, 1124.6]
x1 = [1., 1477., 915.]
xmin1 = [0., 750., 750.]
xmax1 = [1., 2250., 2250.]

"""
# [during_depth]
result = optimize_explore(x0, xmin0, xmax0, delta_weights_error, build_rule_waveform0, cos_mod_weight['CA3']-1.,
                          field_center1)
polished_result = optimize_polish(result['x'], xmin0, xmax0, delta_weights_error, build_rule_waveform0,
                                  cos_mod_weight['CA3']-1., field_center1)
delta_weights_error(polished_result['x'], xmin0, xmax0, build_rule_waveform0, cos_mod_weight['CA3']-1., field_center1,
                    plot=True)

# [during_depth, pre_pot_dur, post_pot_dur]
result = optimize_explore(x1, xmin1, xmax1, delta_weights_error, build_rule_waveform1, cos_mod_weight['CA3']-1.,
                          field_center1)
polished_result = optimize_polish(result['x'], xmin1, xmax1, delta_weights_error, build_rule_waveform1,
                                  cos_mod_weight['CA3']-1., field_center1)
delta_weights_error(polished_result['x'], xmin1, xmax1, build_rule_waveform0, cos_mod_weight['CA3']-1., field_center1,
                    plot=True)
"""