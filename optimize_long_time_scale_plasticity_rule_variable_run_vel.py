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

synapses_seed = 0
num_exc_syns = 1600
num_inh_syns = 600
# allows parallel computation of multiple trials for the same spines with the same peak_locs, but with different
# input spike trains and stochastic synapses for each trial
trial_seed = 0
field1_loc = 0.5
# .pkl file contains traces for position vs. time, ramp vs. position, and difference between 1st and 2nd induction
if len(sys.argv) > 1:
    experimental_filename = str(sys.argv[1])
else:
    experimental_filename = '102716 katie newcell5 saved output'


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


def get_dynamic_theta_phase_force(phase_ranges, peak_loc, stim_x, x_interp, stim_t, dx):
    """
    Expects a list of tuples containing times and phases relative to peak_loc and the non-modulated phase preference
    (zero degrees). Returns a waveform of phase vs time.
    :param phase_ranges: list of tuple (ms, degrees)
    :param peak_loc: float (cm)
    :param stim_x: array (just one item from the outer list)
    :param x_interp: array (just one item from the outer list)
    :param stim_t: array (just one item from the outer list)
    :param dx: float (cm)
    :return: :class: 'np.array'
    """
    end_x_val = phase_ranges[-1][0] + peak_loc
    start_phase_val = phase_ranges[0][1] * 2. * np.pi / 360.  # convert degrees to radians
    end_phase_val = phase_ranges[-1][1] * 2. * np.pi / 360.  # convert degrees to radians
    phase_force = np.ones_like(stim_x) * start_phase_val
    phase_force[np.where(stim_x >= end_x_val)[0]] *= end_phase_val
    for i in range(len(phase_ranges) - 1):
        x0 = phase_ranges[i][0] + peak_loc
        x1 = phase_ranges[i + 1][0] + peak_loc
        phase0 = phase_ranges[i][1] * 2. * np.pi / 360.  # convert degrees to radians
        phase1 = phase_ranges[i + 1][1] * 2. * np.pi / 360.  # convert degrees to radians
        del_x = x1 - x0
        del_phase = phase1 - phase0
        if abs(del_x) > 0.:
            this_x_piece = np.arange(x0, x1, dx)
            this_indexes = np.where((stim_x >= x0) & (stim_x < x1))[0]
            this_x_interp_piece = stim_x[this_indexes]
            if abs(del_phase) > 0.:
                d_phase = del_phase / del_x * dx
                this_range_piece = np.arange(phase0, phase1, d_phase)
                this_range_interp_piece = np.interp(this_x_interp_piece, this_x_piece, this_range_piece)
            else:
                this_range_interp_piece = np.ones(len(this_x_interp_piece)) * phase0
            phase_force[this_indexes] = this_range_interp_piece
    before_start = peak_loc - track_length * 0.5
    after_end = peak_loc + track_length * 0.5
    after_end_index = np.where(stim_x >= after_end)[0][0]
    if before_start < 0.:
        before = np.array(phase_force[np.where((stim_x >= before_start) & (stim_x < 0.))[0]])
        phase_force[after_end_index:after_end_index+len(before)] = before
    elif after_end > track_length:
        after = np.array(phase_force[2*len(x_interp)-int(track_equilibrate/dt):after_end_index])
        phase_force[len(x_interp)-int(track_equilibrate/dt):len(x_interp)-int(track_equilibrate/dt)+len(after)] = after
    return phase_force[len(x_interp)-int(track_equilibrate/dt):len(x_interp)-int(track_equilibrate/dt)+len(stim_t)]


def generate_rate_maps(simiter, stim_x, x_interp, stim_t, global_phase_offset=0.):
    """

    :param simiter: int
    :param stim_x: array (just one item from the outer list)
    :param x_interp: array (just one item from the outer list)
    :param stim_t: array (just one item from the outer list)
    :param global_phase_offset: float
    """
    local_random.seed(simiter)
    rate_maps = []
    for group in ['CA3']:  # stim_exc_syns:
        for i, syn in enumerate(stim_exc_syns[group]):
            peak_loc = peak_locs[group][i]
            gauss_force = excitatory_peak_rate[group] * np.exp(-((stim_x - peak_loc) / gauss_sigma)**2.)
            before = np.array(gauss_force[:len(x_interp)])
            after = np.array(gauss_force[2*len(x_interp)-int(track_equilibrate/dt):])
            gauss_force[len(x_interp):len(x_interp) + len(before)] += before
            gauss_force[len(x_interp)-int(track_equilibrate/dt):len(x_interp)-int(track_equilibrate/dt)+len(after)] \
                += after
            gauss_force = np.array(gauss_force[len(x_interp)-int(track_equilibrate/dt):len(x_interp)-
                                                                                       int(track_equilibrate/dt)+
                                                                                       len(stim_t)])
            if group in excitatory_precession_range:
                phase_force = get_dynamic_theta_phase_force(excitatory_precession_range[group], peak_loc, stim_x,
                                                            x_interp, stim_t, dx)
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
            rate_maps.append(stim_force)
    return rate_maps


def subtract_baseline(waveform, baseline=None):
    """

    :param waveform: array
    :param baseline: float
    :return: array
    """
    if baseline is None:
        baseline = np.mean(waveform[np.where(np.array(waveform) <= np.percentile(waveform, 10.))[0]])
    waveform -= baseline
    return baseline


def get_expected_depolization(rate_maps, weights, x_interp, stim_t, track_duration, baseline=None):
    """
    Take pre-computed rate maps and weights for a set of inputs. Convolve with an EPSP kernel, sum, normalize, and
    downsample to 100 spatial bins.
    :param rate_maps: list of array
    :param weights: list of float
    :param x_interp: array (just one item from the outer list)
    :param stim_t: array (just one item from the outer list)
    :param track_duration: float (just one item from the outer list)
    :param baseline: float
    :return: array
    """
    filter_t = np.arange(0., 200., dt)
    epsp_filter = np.exp(-filter_t/20.) - np.exp(-filter_t/.5)
    epsp_filter /= np.sum(epsp_filter)
    weighted_rate_maps = []
    scaling_factor = 7.15e-4  # generates a predicted 6 mV depolarization from gaussian weights with peak = 2.5
    for i, rate_map in enumerate(rate_maps):
        this_weighted_rate_map = rate_map * weights[i] * scaling_factor
        this_weighted_rate_map = np.convolve(this_weighted_rate_map, epsp_filter)[:len(stim_t)]
        weighted_rate_maps.append(this_weighted_rate_map)
    expected_depolarization = np.sum(weighted_rate_maps, axis=0)
    expected_depolarization = low_pass_filter(expected_depolarization[int(track_equilibrate / dt):
                                                int(track_equilibrate / dt) + len(x_interp)], 2.,
                                                    track_equilibrate+track_duration, dt, 1.)
    expected_depolarization = np.interp(x_bins, x_interp, expected_depolarization)
    baseline = subtract_baseline(expected_depolarization, baseline)
    return expected_depolarization, baseline


NMDA_type = 'NMDA_KIN5'

dt = 1.  # ms
equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # -0.5 * run_vel + 155.  # (ms)
input_field_width = 90.  # cm
track_length = 190.  # cm
dx = track_length / 100.  # cm
x = np.arange(0., track_length+dx/2., dx)
x_bins = x[:100] + dx/2.
default_run_vel = 30.  # cm/s
if experimental_filename is not None:
    experimental_data = read_from_pkl(data_dir+experimental_filename+'.pkl')
    x_t = np.multiply(experimental_data['t'], 1000.)  # convert s to ms
    ramp = experimental_data['ramp']
    difference = experimental_data['difference']
    induction_locs = experimental_data['induction_locs']
else:
    x_t = [np.array(x / default_run_vel * 1000.)]
    ramp = None
    difference = None
    induction_locs = None
t = [np.arange(0., this_x_t[-1]+dt/2., dt) for this_x_t in x_t]  # ms
x_interp = [np.interp(t[i], x_t[i], x) for i in range(len(x_t))]

track_equilibrate = 2. * global_theta_cycle_duration
track_duration = [this_t[-1] for this_t in t]
duration = [equilibrate + track_equilibrate + this_track_duration for this_track_duration in track_duration]
induction_duration = 300.

excitatory_peak_rate = {'CA3': 40., 'ECIII': 40.}
excitatory_theta_modulation_depth = {'CA3': 0.7, 'ECIII': 0.7}
# From Chadwick et al., ELife 2015
excitatory_theta_phase_tuning_factor = {'CA3': 0.8, 'ECIII': 0.8}
excitatory_precession_range = {}  # # (ms, degrees)
excitatory_precession_range['CA3'] = [(-input_field_width*0.65, 0.), (-input_field_width*0.5, 180.),
                                      (-input_field_width*0.35, 180.), (input_field_width*0.35, -180.),
                                      (input_field_width*0.5, -180.), (input_field_width*0.65, 0.)]
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
sim.parameters['input_field_width'] = input_field_width
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

stim_x = []
stim_t = []
extended_stim_t = []
for i in range(len(x_interp)):
    this_stim_x = np.zeros(3*len(x_interp[i]))
    this_stim_x[:len(x_interp[i])] = x_interp[i] - track_length
    this_stim_x[len(x_interp[i]):2*len(x_interp[i])] = x_interp[i]
    this_stim_x[2*len(x_interp[i]):] = x_interp[i] + track_length
    stim_x.append(this_stim_x)
    stim_t.append(np.append(np.arange(-track_equilibrate, 0., dt), t[i]))
    this_extended_stim_t = np.zeros(3*len(t[i]))
    this_extended_stim_t[:len(t[i])] = t[i] - track_duration[i]
    this_extended_stim_t[len(t[i]):2*len(t[i])] = t[i]
    this_extended_stim_t[2*len(t[i]):] = t[i] + track_duration[i]
    extended_stim_t.append(this_extended_stim_t)

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
                                       (extended_peak_locs[group] < track_length))[0]
        peak_locs[group] = extended_peak_locs[group][within_track[group]]

for group in stim_exc_syns:
    for syn in stim_exc_syns[group]:
        all_exc_syns[syn.node.parent.parent.type].remove(syn)

global_phase_offset = 0.  # -2.*np.pi*165./360.

rate_maps = []
for i in range(len(stim_x)):
    rate_maps.append(generate_rate_maps(trial_seed, stim_x[i], x_interp[i], stim_t[i], global_phase_offset))

# modulate the weights of inputs with peak_locs along this stretch of the track
field_center1 = track_length * field1_loc

if ramp is None:
    i = 0
    extended_cos_mod_weight = {}
    cos_mod_weight = {}
    peak_mod_weight = 2.5
    tuning_amp = (peak_mod_weight - 1.)
    start_loc = field_center1 - input_field_width * 1.2 / 2.
    end_loc = field_center1 + input_field_width * 1.2 / 2.

    for group in ['CA3']:  # stim_exc_syns:
        extended_cos_mod_weight[group] = np.zeros_like(extended_peak_locs[group])
        within_field = np.where((extended_peak_locs[group] >= start_loc) & (extended_peak_locs[group] <= end_loc))[0]
        extended_cos_mod_weight[group][within_field] = tuning_amp / 2. * np.cos(2. * np.pi / (input_field_width * 1.2) *
                                                                                (extended_peak_locs[group][
                                                                                     within_field] -
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
    ramp = [get_expected_depolization(rate_maps[i], cos_mod_weight['CA3'], x_interp[i], stim_t[i], track_duration[i])]
else:
    baseline = None
    for this_ramp in ramp:
        baseline = subtract_baseline(this_ramp, baseline)

if induction_locs is None:
    induction_locs = field_center1 + 5.  # optimize for backward shift of peak of place field


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


def build_rule_waveform0(x, induction_loc, x_interp, stim_x, stim_t, extended_stim_t, plot=False):
    """
    Construct a normalized rule waveform with the following components:
    0) phase-dependent potentiation during current injection
    :param x: array: [during_depth, phase_offset]
    :param induction_loc: float (just one item from the outer list)
    :param x_interp: array (just one item from the outer list)
    :param stim_x: array (just one item from the outer list)
    :param stim_t: array (just one item from the outer list)
    :param extended_stim_t: array (just one item from the outer list)
    :param plot: bool
    :return: array
    """
    during_depth = x[0]
    during_pot_amp = (1. - during_depth) / 2.
    rule_phase_offset = x[1] * 2. * np.pi / 360.

    induction_start_index = np.where(stim_x >= induction_loc)[0][0]
    induction_start = extended_stim_t[induction_start_index]
    rule_theta = during_pot_amp * np.cos(2. * np.pi * extended_stim_t / global_theta_cycle_duration -
                                               global_phase_offset - rule_phase_offset) + during_pot_amp + during_depth
    rule_waveform = np.zeros_like(extended_stim_t)
    left = induction_start_index
    right = left + int(induction_duration / dt)
    rule_waveform[left:right] = np.array(rule_theta[left:right])

    before_rule = np.array(rule_waveform[:len(x_interp)])
    after_rule = np.array(rule_waveform[2*len(x_interp):])
    within_rule = np.array(rule_waveform[len(x_interp):2*len(x_interp)])
    rule_waveform = np.zeros_like(stim_t)
    rule_waveform[int(track_equilibrate / dt):int(track_equilibrate / dt)+len(x_interp)] = before_rule + after_rule + \
                                                                                           within_rule
    ref_theta = 0.5 * np.cos(2. * np.pi * stim_t / global_theta_cycle_duration - global_phase_offset) + 1.75

    if plot:
        left = int((track_equilibrate + induction_start) / dt)
        right = left + int(induction_duration / dt)
        fig, axes = plt.subplots(1)
        axes.plot(stim_t[left:right] - induction_start, rule_waveform[left:right], label='Plasticity rule')
        axes.plot(stim_t[left:right] - induction_start, ref_theta[left:right], label='LFP theta')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Time (ms)')
        axes.set_ylabel('Normalized change in synaptic weight')
        axes.set_ylim([math.floor(np.min(rule_waveform)) - 0.25, 2.5])
        axes.set_xlim([0., induction_duration])
        clean_axes(axes)
        plt.show()
        plt.close()

    return rule_waveform


def build_rule_waveform1(x, induction_loc, x_interp, stim_x, stim_t, extended_stim_t, plot=False):
    """
    Construct a normalized rule waveform with the following components:
    0) phase-dependent potentiation during current injection
    1) time-dependent potentiation before current injection
    2) time-dependent potentiation after current injection
    :param x: array: [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur]
    :param induction_loc: float (just one item from the outer list)
    :param x_interp: array (just one item from the outer list)
    :param stim_x: array (just one item from the outer list)
    :param stim_t: array (just one item from the outer list)
    :param extended_stim_t: array (just one item from the outer list)
    :param plot: bool
    :return: array
    """
    during_depth = x[0]
    during_pot_amp = (1. - during_depth) / 2.
    rule_phase_offset = x[1]
    pot_amp = x[2]
    pre_pot_dur = x[3]
    post_pot_dur = x[4]

    induction_start_index = np.where(stim_x >= induction_loc)[0][0]
    induction_start = extended_stim_t[induction_start_index]
    rule_theta = during_pot_amp * np.cos(2. * np.pi * extended_stim_t / global_theta_cycle_duration -
                                  global_phase_offset - rule_phase_offset) + during_pot_amp + during_depth
    rule_waveform = np.zeros_like(extended_stim_t)
    left = induction_start_index - int(pre_pot_dur / dt)
    right = induction_start_index
    pre_pot_phase_offset = 2. * np.pi * ((induction_start - pre_pot_dur) % (pre_pot_dur * 2.)) / (pre_pot_dur * 2.)
    pre_pot_waveform = pot_amp * (0.5 * np.cos(2. * np.pi * extended_stim_t / (pre_pot_dur * 2.) - np.pi -
                                                    pre_pot_phase_offset) + 0.5)
    rule_waveform[left:right] = np.array(pre_pot_waveform[left:right])
    left = right
    right = left + int(induction_duration / dt)
    rule_waveform[left:right] = np.array(rule_theta[left:right])
    left = right
    right = left + int(post_pot_dur / dt)
    post_pot_phase_offset = 2. * np.pi * ((induction_start + induction_duration) % (post_pot_dur * 2.)) / \
                            (post_pot_dur * 2.)
    post_waveform1 = pot_amp * (0.5 * np.cos(2. * np.pi * extended_stim_t / (post_pot_dur * 2.) -
                                                  post_pot_phase_offset) + 0.5)
    rule_waveform[left:right] = np.array(post_waveform1[left:right])

    ref_theta = 0.5 * np.cos(2. * np.pi * extended_stim_t / global_theta_cycle_duration - global_phase_offset) + 1.75
    if plot:
        left = induction_start_index - int(5000. / dt)
        right = induction_start_index + int(5000. / dt)
        fig, axes = plt.subplots(1)
        axes.plot(extended_stim_t[left:right] - induction_start, rule_waveform[left:right], label='Plasticity rule')
        axes.plot(extended_stim_t[left:right] - induction_start, ref_theta[left:right], label='LFP theta')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Time (ms)')
        axes.set_ylabel('Normalized change in synaptic weight')
        axes.set_ylim([math.floor(np.min(rule_waveform)) - 0.25, 2.5])
        axes.set_xlim([-5000., 5000.])
        clean_axes(axes)
        plt.show()
        plt.close()

    before_rule = np.array(rule_waveform[:len(x_interp)])
    after_rule = np.array(rule_waveform[2 * len(x_interp):])
    within_rule = np.array(rule_waveform[len(x_interp):2 * len(x_interp)])
    rule_waveform = np.zeros_like(stim_t)
    rule_waveform[int(track_equilibrate / dt):int(track_equilibrate / dt) + len(x_interp)] = before_rule + \
                                                                                             after_rule + \
                                                                                             within_rule
    return rule_waveform


def build_rule_waveform2(x, induction_loc, x_interp, stim_x, stim_t, extended_stim_t, plot=False):
    """
    Construct a normalized rule waveform with the following components:
    0) phase-dependent potentiation during current injection
    1) time-dependent potentiation before current injection
    2) time-dependent potentiation after current injection
    3) time-dependent depotentiation after current injection
    :param x: array: [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, post_depot_depth, post_depot_dur]
    :param induction_loc: float (just one item from the outer list)
    :param x_interp: array (just one item from the outer list)
    :param stim_x: array (just one item from the outer list)
    :param stim_t: array (just one item from the outer list)
    :param extended_stim_t: array (just one item from the outer list)
    :param plot: bool
    :return: array
    """
    during_depth = x[0]
    during_pot_amp = (1. - during_depth) / 2.
    rule_phase_offset = x[1]
    pot_amp = x[2]
    pre_pot_dur = x[3]
    post_pot_dur = x[4]
    post_depot_depth = x[5]
    post_depot_dur = x[6]

    induction_start_index = np.where(stim_x >= induction_loc)[0][0]
    induction_start = extended_stim_t[induction_start_index]
    rule_theta = during_pot_amp * np.cos(2. * np.pi * extended_stim_t / global_theta_cycle_duration -
                                  global_phase_offset - rule_phase_offset) + during_pot_amp + during_depth
    rule_waveform = np.zeros_like(extended_stim_t)
    left = induction_start_index - int(pre_pot_dur / dt)
    right = induction_start_index
    pre_pot_phase_offset = 2. * np.pi * ((induction_start - pre_pot_dur) % (pre_pot_dur * 2.)) / (pre_pot_dur * 2.)
    pre_pot_waveform = pot_amp * (0.5 * np.cos(2. * np.pi * extended_stim_t / (pre_pot_dur * 2.) - np.pi -
                                                    pre_pot_phase_offset) + 0.5)
    rule_waveform[left:right] = np.array(pre_pot_waveform[left:right])
    left = right
    right = left + int(induction_duration / dt)
    rule_waveform[left:right] = np.array(rule_theta[left:right])
    left = right
    right = left + int(post_pot_dur / dt)
    post_pot_phase_offset = 2. * np.pi * ((induction_start + induction_duration) % (post_pot_dur * 2.)) / \
                            (post_pot_dur * 2.)
    post_pot_waveform = pot_amp * (0.5 * np.cos(2. * np.pi * extended_stim_t / (post_pot_dur * 2.) -
                                                  post_pot_phase_offset) + 0.5)
    rule_waveform[left:right] = np.array(post_pot_waveform[left:right])
    left = right
    right = left + int(post_depot_dur / dt)
    post_depot_phase_offset = 2. * np.pi * ((induction_start + induction_duration + post_pot_dur) % post_depot_dur) / \
                              post_depot_dur
    post_depot_waveform = post_depot_depth * 0.5 * np.cos(2. * np.pi * extended_stim_t / post_depot_dur -
                                                          post_depot_phase_offset) - post_depot_depth * 0.5
    rule_waveform[left:right] = np.array(post_depot_waveform[left:right])

    ref_theta = 0.5 * np.cos(2. * np.pi * extended_stim_t / global_theta_cycle_duration - global_phase_offset) + 1.75
    if plot:
        left = induction_start_index - int(5000. / dt)
        right = induction_start_index + int(5000. / dt)
        fig, axes = plt.subplots(1)
        axes.plot(extended_stim_t[left:right] - induction_start, rule_waveform[left:right], label='Plasticity rule')
        axes.plot(extended_stim_t[left:right] - induction_start, ref_theta[left:right], label='LFP theta')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Time (ms)')
        axes.set_ylabel('Normalized change in synaptic weight')
        axes.set_ylim([math.floor(np.min(rule_waveform)) - 0.25, 2.5])
        axes.set_xlim([-5000., 5000.])
        clean_axes(axes)
        plt.show()
        plt.close()

    before_rule = np.array(rule_waveform[:len(x_interp)])
    after_rule = np.array(rule_waveform[2 * len(x_interp):])
    within_rule = np.array(rule_waveform[len(x_interp):2 * len(x_interp)])
    rule_waveform = np.zeros_like(stim_t)
    rule_waveform[int(track_equilibrate / dt):int(track_equilibrate / dt) + len(x_interp)] = before_rule + \
                                                                                             after_rule + \
                                                                                             within_rule
    return rule_waveform


def calculate_delta_weights(rule_waveform, i):
    """
    Given any rule waveform, compute the resulting delta_weights (without saturation).
    :param rule_waveform: array
    :param i: int: index of various location and time waves to use
    :param plot: bool
    :return: delta_weights: array, rule_gain: float
    """
    group = 'CA3'
    delta_weights = np.zeros_like(peak_locs[group])
    for i, stim_force in enumerate(rate_maps[i]):
        this_force = 0.001 * dt * np.multiply(stim_force, rule_waveform)
        this_area = np.trapz(this_force, dx=dt)
        delta_weights[i] = this_area
    rule_gain = 1. / np.max(delta_weights)
    delta_weights /= rule_gain
    return delta_weights, rule_gain


def ramp_error(x, xmin, xmax, rule_function, ramp, i, baseline=None, plot=False, full_output=False):
    """

    :param x: array [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, post_depot_depth, post_depot_dur]
    :param xmin: array
    :param xmax: array
    :param rule_function: callable
    :param ramp: array
    :param i: int: list index of various location and time waves to use
    :param baseline: float
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
    if not check_bounds(x, xmin, xmax):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    rule_waveform = rule_function(x, induction_locs[i], x_interp[i], stim_x[i], stim_t[i], extended_stim_t[i], plot)
    delta_weights, rule_gain = calculate_delta_weights(rule_waveform, i)
    new_ramp, new_baseline = get_expected_depolization(rate_maps[i], delta_weights + 1., x_interp[i], stim_t[i],
                                                       track_duration[i], baseline)
    ramp_peak = np.max(ramp)
    # ramp_peak_index = np.where(ramp == ramp_peak)[0][0]
    # ramp_gain = ramp_peak / new_ramp[ramp_peak_index]
    ramp_gain = ramp_peak / np.max(new_ramp)
    delta_weights *= ramp_gain
    rule_gain *= ramp_gain
    new_ramp *= ramp_gain
    Err = 0.
    for j in range(len(ramp)):
        Err += ((ramp[j] - new_ramp[j]) / 0.05) ** 2.

    if plot:
        fig, axes = plt.subplots(1)
        axes.scatter(peak_locs['CA3'], delta_weights + 1., color='r', label='Weights from plasticity rule')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Normalized synaptic weight')
        axes.set_ylim([math.floor(np.min(delta_weights + 1.)), math.ceil(np.max(delta_weights + 1.))])
        axes.set_xlim([0., track_length])
        clean_axes(axes)
        plt.show()
        plt.close()

        fig, axes = plt.subplots(1)
        axes.plot(x_bins, ramp, label='Experiment')
        axes.plot(x_bins, new_ramp, label='Model')
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Depolarization (mV)')
        axes.set_xlim([0., track_length])
        axes.set_ylim([math.floor(np.min(ramp)) - 0.5, math.ceil(np.max(ramp)) + 0.5])
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes)
        plt.show()
        plt.close()

    formatted_x = '[' + ', '.join(['%.2f' % xi for xi in x]) + ']'
    print 'x: %s, Err: %.4E, Rule gain: %.4E' % (formatted_x, Err, rule_gain)
    if full_output:
        return rule_waveform, delta_weights, new_ramp, new_baseline
    else:
        hist.x.append(x)
        hist.Err.append(Err)
        return Err


def optimize_polish(x, xmin, xmax, error_function, rule_function, ramp, i, baseline=None, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param rule_function: callable
    :param ramp: array
    :param i: int: index of various location and time waves to use
    :param baseline: float
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 200

    result = optimize.minimize(error_function, x, method='Nelder-Mead', options={'ftol': 1e-3,
                                                    'xtol': 1e-3, 'disp': True, 'maxiter': maxfev},
                               args=(xmin, xmax, rule_function, ramp, i, baseline))
    formatted_x = '['+', '.join(['%.2E' % xi for xi in result.x])+']'
    print 'Process: %i completed optimize_polish after %i iterations with Error: %.4E and x: %s' % (os.getpid(),
                                                                            result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


def optimize_explore(x, xmin, xmax, error_function, rule_function, ramp, i, baseline=None, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param rule_function: callable
    :param ramp: array
    :param i: int: index of various location and time waves to use
    :param baseline: float
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 400

    take_step = Normalized_Step(x, xmin, xmax)
    minimizer_kwargs = dict(method=null_minimizer, args=(xmin, xmax, rule_function, ramp, i, baseline))
    result = optimize.basinhopping(error_function, x, niter=maxfev, niter_success=maxfev/2,
                                   disp=True, interval=20, minimizer_kwargs=minimizer_kwargs,
                                   take_step=take_step)
    formatted_x = '['+', '.join(['%.2E' % xi for xi in result.x])+']'
    print 'Process: %i completed optimize_explore after %i iterations with Error: %.4E and x: %s' % (os.getpid(),
                                                                            result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


rule_waveforms = []
delta_weights = []
model_ramp = []


# [during_depth, phase_offset]
# x0 = [0.7, 45.]
x0 = [0., 45.]  # Err: 1.932E+05, Rule gain:
xmin0 = [0., 45.]
xmax0 = [1., 90.]

# [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur]
# x1 = [2.65E-01, 6.24E+01, 1.00, 2.5E+03, 5.E+02]
# x1 = [2.99E-04, 7.56E+01, 1.00E+00, 2.80E+03, 8.95E+02]  # Err: 3.649E+02, Rule gain: 1.8212E-04, normalized to
                                                         # experiment peak_loc
x1 = [0.00, 75.54, 1.00, 2826.48, 923.00]  # Err: 1.3359E+03, Rule gain: 1.7876E-04
xmin1 = [0., 45., 0., 1000., 500.]
xmax1 = [1., 90., 1., 3000., 2000.]

# [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, post_depot_depth, post_depot_dur]
# x2 = [0.00, 75.54, 1.00, 2826.48, 923.00, 2., 600.]
# x2 = [0.22, 81.98, 0.75, 2999.97, 500.00, 3.00, 861.61]
x2 = [0.01, 88.49, 0.97, 2999.73, 200.02, 2.17, 1492.99]
xmin2 = [0., 45., 0., 1000., 100., 0., 200.]
xmax2 = [1., 90., 1., 4000., 1500., 5., 2500.]

# [during_depth, phase_offset]
# induction_index = 0
# ramp_error(x0, xmin0, xmax0, build_rule_waveform0, ramp[induction_index], induction_index, plot=True)

# induction_index = 0
# [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur]
# ramp_error(x1, xmin1, xmax1, build_rule_waveform1, ramp[induction_index], induction_index, plot=True)

# induction_index = 1
# [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, post_depot_depth, post_depot_dur]
# ramp_error(x2, xmin2, xmax2, build_rule_waveform2, difference, induction_index, plot=True)

"""
induction_index = 0
result = optimize_explore(x0, xmin0, xmax0, ramp_error, build_rule_waveform0, ramp[induction_index], induction_index)
polished_result = optimize_polish(result['x'], xmin0, xmax0, ramp_error, build_rule_waveform0, ramp[induction_index],
                                  induction_index)
ramp_error(polished_result['x'], xmin0, xmax0, build_rule_waveform0, ramp[induction_index], induction_index, plot=True)
"""

"""
induction_index = 0
result = optimize_explore(x1, xmin1, xmax1, ramp_error, build_rule_waveform1, ramp[induction_index], induction_index,
                          maxfev=700)
polished_result = optimize_polish(result['x'], xmin1, xmax1, ramp_error, build_rule_waveform1, ramp[induction_index],
                                  induction_index)
ramp_error(polished_result['x'], xmin1, xmax1, build_rule_waveform1, ramp[induction_index], induction_index, plot=True)
"""

induction_index = 0
this_rule_waveform, this_delta_weights, this_model_ramp, baseline = ramp_error(x1, xmin1, xmax1, build_rule_waveform1,
                                                                             ramp[induction_index],
                                                                             induction_index, plot=True,
                                                                             full_output=True)
rule_waveforms.append(this_rule_waveform)
delta_weights.append(this_delta_weights)
model_ramp.append(this_model_ramp)

induction_index = 1

result = optimize_explore(x2, xmin2, xmax2, ramp_error, build_rule_waveform2, difference, induction_index, baseline,
                          maxfev=700)
polished_result = optimize_polish(result['x'], xmin2, xmax2, ramp_error, build_rule_waveform2, difference,
                                  induction_index, baseline)
this_rule_waveform, this_delta_weights, this_model_ramp, baseline = ramp_error(polished_result['x'], xmin2, xmax2,
                                                                               build_rule_waveform2, difference,
                                                                               induction_index, baseline=baseline,
                                                                               plot=True, full_output=True)
rule_waveforms.append(this_rule_waveform)
delta_weights.append(this_delta_weights)
model_ramp.append(this_model_ramp)
