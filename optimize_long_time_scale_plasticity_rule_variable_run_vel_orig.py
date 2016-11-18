__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import random
import sys
mpl.rcParams['font.size'] = 30.
mpl.rcParams['legend.fontsize'] = 30.
mpl.rcParams['figure.figsize'] = 6, 4.3

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
    cell_id = int(sys.argv[1])
else:
    cell_id = None

experimental_filename = {1: '110716 katie newcell1 saved output', 2: '110716 katie newcell2 saved output',
                         3: '110716 katie newcell3 saved output', 4: '110716 katie newcell4 saved output',
                         5: '110716 katie newcell5 saved output'}

rule_max_timescale = 9000.


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


def get_expected_depolization(rate_maps, weights, x_interp, stim_t, track_duration):
    """
    Take pre-computed rate maps and weights for a set of inputs. Convolve with an EPSP kernel, sum, normalize, and
    downsample to 100 spatial bins.
    :param rate_maps: list of array
    :param weights: list of float
    :param x_interp: array (just one item from the outer list)
    :param stim_t: array (just one item from the outer list)
    :param track_duration: float (just one item from the outer list)
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
    return expected_depolarization


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
if cell_id is not None:
    experimental_data = read_from_pkl(data_dir+experimental_filename[cell_id]+'.pkl')
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
    ramp_baseline_indexes = None
    for this_ramp in ramp:
        if baseline is None:
            ramp_baseline_indexes = np.where(np.array(this_ramp) <= np.percentile(this_ramp, 10.))[0]
            baseline = np.mean(this_ramp[ramp_baseline_indexes])
        ignore = subtract_baseline(this_ramp, baseline)

if induction_locs is None:
    induction_locs = [field_center1 + 5.]  # optimize for backward shift of peak of place field


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
    :return: array, array
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
    ref_theta = 0.5 * np.cos(2. * np.pi * extended_stim_t / global_theta_cycle_duration - global_phase_offset) + 1.75

    rule_waveform_offset = rule_waveform[induction_start_index-int(rule_max_timescale/dt):induction_start_index+
                                                                                          int(rule_max_timescale/dt)]

    if plot:
        left = induction_start_index
        right = left + int(induction_duration / dt)
        fig, axes = plt.subplots(1)
        axes.plot(extended_stim_t[left:right] - induction_start, rule_waveform[left:right], label='Plasticity rule')
        axes.plot(extended_stim_t[left:right] - induction_start, ref_theta[left:right], label='LFP theta')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Time (ms)')
        axes.set_ylabel('Relative change in synaptic weight per spike')
        axes.set_ylim([math.floor(np.min(rule_waveform)) - 0.25, 2.5])
        axes.set_xlim([0., induction_duration])
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

    return rule_waveform, rule_waveform_offset


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
    rule_waveform_offset = rule_waveform[
                           induction_start_index - int(rule_max_timescale / dt):induction_start_index +
                                                                                int(rule_max_timescale / dt)]

    if plot:
        left = induction_start_index - int(rule_max_timescale / dt)
        right = induction_start_index + int(rule_max_timescale / dt)
        fig, axes = plt.subplots(1)
        axes.plot(extended_stim_t[left:right] - induction_start, rule_waveform[left:right], label='Plasticity rule')
        axes.plot(extended_stim_t[left:right] - induction_start, ref_theta[left:right], label='LFP theta')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Time (s)')
        axes.set_ylabel('Relative change in synaptic weight per spike')
        axes.set_ylim([math.floor(np.min(rule_waveform)) - 0.25, 2.5])
        # axes.set_xlim([-rule_max_timescale, rule_max_timescale])
        axes.set_xlim(-4500., 2000.)
        axes.set_xticks([-4000., -3000., -2000., -1000., 0., 1000., 2000.])
        axes.set_xticklabels([-4, -3, -2, -1, 0, 1, 2])
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
    return rule_waveform, rule_waveform_offset


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
    """
    left = right
    right = left + int(post_depot_dur / dt)
    post_depot_phase_offset = 2. * np.pi * ((induction_start + induction_duration + post_pot_dur) % post_depot_dur) / \
                              post_depot_dur
    post_depot_waveform = post_depot_depth * 0.5 * np.cos(2. * np.pi * extended_stim_t / post_depot_dur -
                                                          post_depot_phase_offset) - post_depot_depth * 0.5
    rule_waveform[left:right] = np.array(post_depot_waveform[left:right])
    """
    right = left + int(post_depot_dur / dt)
    post_depot_phase_offset = 2. * np.pi * ((induction_start + induction_duration) % post_depot_dur) / \
                              post_depot_dur
    post_depot_waveform = post_depot_depth * 0.5 * np.cos(2. * np.pi * extended_stim_t / post_depot_dur -
                                                          post_depot_phase_offset) - post_depot_depth * 0.5
    rule_waveform[left:right] = np.add(rule_waveform[left:right], post_depot_waveform[left:right])
    ref_theta = 0.5 * np.cos(2. * np.pi * extended_stim_t / global_theta_cycle_duration - global_phase_offset) + 1.75
    rule_waveform_offset = rule_waveform[
                           induction_start_index - int(rule_max_timescale / dt):induction_start_index +
                                                                                int(rule_max_timescale / dt)]

    if plot:
        left = induction_start_index - int(rule_max_timescale / dt)
        right = induction_start_index + int(rule_max_timescale / dt)
        fig, axes = plt.subplots(1)
        axes.plot(extended_stim_t[left:right] - induction_start, rule_waveform[left:right], label='Plasticity rule')
        axes.plot(extended_stim_t[left:right] - induction_start, ref_theta[left:right], label='LFP theta')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Time (ms)')
        axes.set_ylabel('Relative change in synaptic weight per spike')
        axes.set_ylim([math.floor(np.min(rule_waveform)) - 0.25, 2.5])
        axes.set_xlim([-rule_max_timescale, rule_max_timescale])
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
    return rule_waveform, rule_waveform_offset


def calculate_plasticity_signal(rule_waveform, induction_index):
    """
    Given any rule waveform, compute the resulting plasticity_signal (without saturation).
    :param rule_waveform: array
    :param induction_index: int: index of various location and time waves to use
    :param plot: bool
    :return: plasticity_signal: array, rule_gain: float
    """
    group = 'CA3'
    plasticity_signal = np.zeros_like(peak_locs[group])
    for j, stim_force in enumerate(rate_maps[induction_index]):
        this_force = 0.001 * dt * np.multiply(stim_force, rule_waveform)
        this_area = np.trapz(this_force, dx=dt)
        plasticity_signal[j] = this_area
    return plasticity_signal


def calculate_delta_weights_orig(plasticity_signal, orig_weights=None, depot_depth=1.):
    """
    Implement a BCM-like transformation of plasticity signal into change in synaptic weight dependent on initial weight.
    :param plasticity_signal: array
    :param orig_weights: array
    :param depot_depth: float
    :return: array
    """
    if orig_weights is None:
        orig_weights = np.zeros_like(peak_locs['CA3'])
    this_filter = lambda this_plasticity_signal, this_w: \
        this_plasticity_signal - this_w if this_plasticity_signal >= this_w else \
            depot_depth * this_w / 2. * np.cos(2. * np.pi * this_plasticity_signal / this_w) - depot_depth * this_w / 2.
    this_filter = np.vectorize(this_filter)
    return this_filter(plasticity_signal, orig_weights) + orig_weights


def calculate_delta_weights(plasticity_signal, orig_weights=None, depot_depth=1., pot_peak=20.):
    """
    Implement a BCM-like transformation of plasticity signal into change in synaptic weight dependent on initial weight.
    :param plasticity_signal: array
    :param orig_weights: array
    :param depot_depth: float
    :param pot_peak: float
    :return: array
    """
    if orig_weights is None:
        orig_weights = np.zeros_like(peak_locs['CA3'])
    this_filter = lambda this_plasticity_signal, this_w: 0. if this_plasticity_signal == 0. \
        else 0.5 * depot_depth * this_w * np.cos(2. * np.pi * this_plasticity_signal / this_w) - \
             0.5 * depot_depth * this_w \
        if 0 < this_plasticity_signal <= this_w / 2. \
        else 0.5 * pot_peak * np.cos(np.pi * (this_plasticity_signal - this_w / 2.) / pot_peak - np.pi) + \
             0.5 * pot_peak - depot_depth * this_w \
        if this_plasticity_signal <= pot_peak + this_w / 2. else pot_peak - depot_depth * this_w
    this_filter = np.vectorize(this_filter)
    return this_filter(plasticity_signal, orig_weights) + orig_weights


def ramp_error_orig(x, xmin, xmax, rule_function, ramp, induction_index, orig_weights=None, baseline_indexes=None,
               baseline=None, plot=False, full_output=False):
    """

    :param x: array [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, post_depot_depth, post_depot_dur]
    :param xmin: array
    :param xmax: array
    :param rule_function: callable
    :param ramp: array
    :param induction_index: int: list index of various location and time waves to use
    :param orig_weights: array
    :param baseline_indexes: list
    :param baseline: float
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
    if not check_bounds(x, xmin, xmax):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    rule_waveform, rule_waveform_offset = rule_function(x, induction_locs[induction_index], x_interp[induction_index],
                                                        stim_x[induction_index], stim_t[induction_index],
                                                        extended_stim_t[induction_index], plot)
    delta_weights = calculate_delta_weights(rule_waveform, induction_index)
    if orig_weights is not None:
        this_weights = delta_weights + orig_weights
    else:
        this_weights = delta_weights
    new_ramp = get_expected_depolization(rate_maps[induction_index], this_weights + 1., x_interp[induction_index],
                                         stim_t[induction_index], track_duration[induction_index])
    if baseline is None:
        if baseline_indexes is None:
            new_baseline = subtract_baseline(new_ramp)
        else:
            new_baseline = np.mean(new_ramp[baseline_indexes])
            new_ramp -= new_baseline
    else:
        new_baseline = subtract_baseline(new_ramp, baseline)

    ramp_peak = np.max(ramp)
    rule_gain = ramp_peak / np.max(new_ramp)
    new_ramp *= rule_gain
    rule_waveform_offset *= rule_gain
    Err = 0.
    for j in range(len(ramp)):
        Err += ((ramp[j] - new_ramp[j]) / 0.05) ** 2.

    if plot:
        fig, axes = plt.subplots(1)
        axes.scatter(peak_locs['CA3'], this_weights * rule_gain + 1., color='r', label='Weights from plasticity rule')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Normalized synaptic weight')
        axes.set_ylim([math.floor(np.min(this_weights * rule_gain + 1.)),
                       math.ceil(np.max(this_weights * rule_gain + 1.))])
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
        return rule_waveform_offset, this_weights, new_ramp, new_baseline, rule_gain
    else:
        hist.x.append(x)
        hist.Err.append(Err)
        return Err


def ramp_error0(x, xmin, xmax, ramp, induction_index=None, orig_weights=None, baseline_indexes=None, baseline=None,
                plot=False, full_output=False):
    """
    Calculates a rule_waveform and set of weights to match the first place field induction. Uses baseline_indexes to
    subtract a baseline at the same spatial locations as the minimum of the experimental data.
    :param x: array [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, depot_depth]
    :param xmin: array
    :param xmax: array
    :param rule_function: callable
    :param ramp: array
    :param induction_index: int: list index of various location and time waves to use
    :param orig_weights: array
    :param baseline_indexes: list
    :param baseline: float
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
    if not check_bounds(x, xmin, xmax):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    if induction_index is None:
        induction_index = 0
    rule_waveform, rule_waveform_offset = build_rule_waveform0(x, induction_locs[induction_index],
                                                               x_interp[induction_index], stim_x[induction_index],
                                                               stim_t[induction_index],
                                                               extended_stim_t[induction_index], plot)
    plasticity_signal = calculate_plasticity_signal(rule_waveform, induction_index)
    this_weights = plasticity_signal
    if orig_weights is not None:
        this_weights = this_weights + orig_weights
    new_ramp = get_expected_depolization(rate_maps[induction_index], this_weights + 1., x_interp[induction_index],
                                         stim_t[induction_index], track_duration[induction_index])
    if baseline is None:
        if baseline_indexes is None:
            new_baseline = subtract_baseline(new_ramp)
        else:
            new_baseline = np.mean(new_ramp[baseline_indexes])
            new_ramp -= new_baseline
    else:
        new_baseline = subtract_baseline(new_ramp, baseline)

    ramp_peak = np.max(ramp)
    rule_gain = ramp_peak / np.max(new_ramp)
    new_ramp *= rule_gain
    rule_waveform_offset *= rule_gain
    Err = 0.
    for j in range(len(ramp)):
        Err += ((ramp[j] - new_ramp[j]) / 0.05) ** 2.

    if plot:
        fig, axes = plt.subplots(1)
        axes.scatter(peak_locs['CA3'], this_weights * rule_gain + 1., color='r', label='Weights from plasticity rule')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Normalized synaptic weight')
        axes.set_ylim([math.floor(np.min(this_weights * rule_gain + 1.)),
                       math.ceil(np.max(this_weights * rule_gain + 1.))])
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
        return rule_waveform_offset, this_weights, new_ramp, new_baseline, rule_gain
    else:
        hist.x.append(x)
        hist.Err.append(Err)
        return Err


def ramp_error1(x, xmin, xmax, ramp, induction_index=None, orig_weights=None, baseline_indexes=None, baseline=None,
                plot=False, full_output=False):
    """
    Calculates a rule_waveform and set of weights to match the first place field induction. Uses baseline_indexes to
    subtract a baseline at the same spatial locations as the minimum of the experimental data.
    :param x: array [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, depot_depth]
    :param xmin: array
    :param xmax: array
    :param rule_function: callable
    :param ramp: array
    :param induction_index: int: list index of various location and time waves to use
    :param orig_weights: array
    :param baseline_indexes: list
    :param baseline: float
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
    if not check_bounds(x, xmin, xmax):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    if induction_index is None:
        induction_index = 0
    rule_waveform, rule_waveform_offset = build_rule_waveform1(x, induction_locs[induction_index],
                                                               x_interp[induction_index], stim_x[induction_index],
                                                               stim_t[induction_index],
                                                               extended_stim_t[induction_index], plot)
    plasticity_signal = calculate_plasticity_signal(rule_waveform, induction_index)
    this_weights = plasticity_signal
    if orig_weights is not None:
        this_weights = this_weights + orig_weights
    new_ramp = get_expected_depolization(rate_maps[induction_index], this_weights + 1., x_interp[induction_index],
                                         stim_t[induction_index], track_duration[induction_index])
    if baseline is None:
        if baseline_indexes is None:
            new_baseline = subtract_baseline(new_ramp)
        else:
            new_baseline = np.mean(new_ramp[baseline_indexes])
            new_ramp -= new_baseline
    else:
        new_baseline = subtract_baseline(new_ramp, baseline)

    ramp_peak = np.max(ramp)
    rule_gain = ramp_peak / np.max(new_ramp)
    new_ramp *= rule_gain
    rule_waveform_offset *= rule_gain
    Err = 0.
    for j in range(len(ramp)):
        Err += ((ramp[j] - new_ramp[j]) / 0.05) ** 2.

    if plot:
        fig, axes = plt.subplots(1)
        axes.scatter(peak_locs['CA3'], this_weights * rule_gain + 1., color='r', label='Weights from plasticity rule')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Normalized synaptic weight')
        axes.set_ylim([math.floor(np.min(this_weights * rule_gain + 1.)),
                       math.ceil(np.max(this_weights * rule_gain + 1.))])
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
        return rule_waveform_offset, this_weights, new_ramp, new_baseline, rule_gain
    else:
        hist.x.append(x)
        hist.Err.append(Err)
        return Err


def ramp_error2(x, xmin, xmax, ramp, induction_index=None, orig_weights=None, baseline_indexes=None, baseline=None,
                rule_gain=None, plot=False, full_output=False):
    """
    Requires information from the result of optimizing ramp_error1, including orig_weights, baseline, and rule_gain.
    Calculates a separate rule_waveform and set of weights to match the second place field induction.
    :param x: array [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, depot_depth]
    :param xmin: array
    :param xmax: array
    :param rule_function: callable
    :param ramp: array
    :param induction_index: int: list index of various location and time waves to use
    :param orig_weights: array
    :param baseline_indexes: list
    :param baseline: float
    :param rule_gain: float
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
    if not check_bounds(x, xmin, xmax):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    if induction_index is None:
        induction_index = 1
    rule_waveform, rule_waveform_offset = build_rule_waveform1(x[:-1], induction_locs[induction_index],
                                                               x_interp[induction_index], stim_x[induction_index],
                                                               stim_t[induction_index],
                                                               extended_stim_t[induction_index], plot)
    plasticity_signal = calculate_plasticity_signal(rule_waveform, induction_index)
    this_weights = calculate_delta_weights(plasticity_signal, orig_weights, x[-1])
    new_ramp = get_expected_depolization(rate_maps[induction_index], this_weights + 1., x_interp[induction_index],
                                         stim_t[induction_index], track_duration[induction_index])
    if baseline is None:
        if baseline_indexes is None:
            new_baseline = subtract_baseline(new_ramp)
        else:
            new_baseline = np.mean(new_ramp[baseline_indexes])
            new_ramp -= new_baseline
    else:
        new_baseline = subtract_baseline(new_ramp, baseline)

    if rule_gain is None:
        ramp_peak = np.max(ramp)
        rule_gain = ramp_peak / np.max(new_ramp)
    new_ramp *= rule_gain
    rule_waveform_offset *= rule_gain

    Err = 0.
    for j in range(len(ramp)):
        Err += ((ramp[j] - new_ramp[j]) / 0.05) ** 2.

    if plot:
        fig, axes = plt.subplots(1)
        axes.scatter(peak_locs['CA3'], this_weights * rule_gain + 1., color='r', label='Weights from plasticity rule')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Normalized synaptic weight')
        axes.set_ylim([math.floor(np.min(this_weights * rule_gain + 1.)),
                       math.ceil(np.max(this_weights * rule_gain + 1.))])
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
        return rule_waveform_offset, this_weights, new_ramp, new_baseline, rule_gain
    else:
        hist.x.append(x)
        hist.Err.append(Err)
        return Err


def ramp_error3(x, xmin, xmax, ramp, induction_index=None, orig_weights=None, baseline_indexes=None, baseline=None,
                rule_gain=None, plot=False, full_output=False):
    """
    Attemps to find a single rule_waveform to account for both the first and second place field inductions.
    :param x: array [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, depot_depth, pot_peak]
    :param xmin: array
    :param xmax: array
    :param rule_function: callable
    :param ramp: list of array
    :param induction_index: int: list index of various location and time waves to use
    :param orig_weights: array
    :param baseline_indexes: list
    :param baseline: float
    :param rule_gain: float
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
    if not check_bounds(x, xmin, xmax):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    induction_index = 0
    rule_waveform0, rule_waveform_offset0 = build_rule_waveform1(x[:-2], induction_locs[induction_index],
                                                               x_interp[induction_index], stim_x[induction_index],
                                                               stim_t[induction_index],
                                                               extended_stim_t[induction_index], plot)
    plasticity_signal0 = calculate_plasticity_signal(rule_waveform0, induction_index)
    this_weights0 = calculate_delta_weights(plasticity_signal0, orig_weights, depot_depth=x[-2], pot_peak=x[-1])
    model_ramp0 = get_expected_depolization(rate_maps[induction_index], this_weights0 + 1., x_interp[induction_index],
                                         stim_t[induction_index], track_duration[induction_index])
    if baseline is None:
        if baseline_indexes is None:
            new_baseline = subtract_baseline(model_ramp0)
        else:
            new_baseline = np.mean(model_ramp0[baseline_indexes])
            model_ramp0 -= new_baseline
    else:
        new_baseline = subtract_baseline(model_ramp0, baseline)

    if rule_gain is None:
        ramp_peak = np.max(ramp[induction_index])
        rule_gain = ramp_peak / np.max(model_ramp0)
    model_ramp0 *= rule_gain
    rule_waveform_offset0 *= rule_gain

    Err = 0.
    for j in range(len(ramp[induction_index])):
        Err += ((ramp[induction_index][j] - model_ramp0[j]) / 0.05) ** 2.

    induction_index = 1
    rule_waveform1, rule_waveform_offset1 = build_rule_waveform1(x[:-2], induction_locs[induction_index],
                                                                 x_interp[induction_index], stim_x[induction_index],
                                                                 stim_t[induction_index],
                                                                 extended_stim_t[induction_index], plot)
    plasticity_signal1 = calculate_plasticity_signal(rule_waveform1, induction_index)
    this_weights1 = calculate_delta_weights(plasticity_signal1, this_weights0, depot_depth=x[-2], pot_peak=x[-1])
    model_ramp1 = get_expected_depolization(rate_maps[induction_index], this_weights1 + 1., x_interp[induction_index],
                                            stim_t[induction_index], track_duration[induction_index])
    model_ramp1 -= new_baseline
    model_ramp1 *= rule_gain
    rule_waveform_offset1 *= rule_gain

    for j in range(len(ramp[induction_index])):
        Err += ((ramp[induction_index][j] - model_ramp1[j]) / 0.05) ** 2.

    if plot:
        colors = ['r', 'k', 'c', 'g']
        fig, axes = plt.subplots(1)
        weight_min = min(np.min(this_weights0), np.min(this_weights1))
        weight_max = max(np.max(this_weights0), np.max(this_weights1))
        for j, this_weights in enumerate([this_weights0, this_weights1]):
            axes.scatter(peak_locs['CA3'], this_weights * rule_gain + 1., color=colors[j], label='Induction '+str(j+1))
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Relative synaptic weight')
        axes.set_ylim([math.floor(weight_min * rule_gain + 1.),
                       math.ceil(weight_max * rule_gain + 1.)])
        axes.set_xlim([0., track_length])
        axes.set_title('Model synaptic weights')
        clean_axes(axes)
        plt.show()
        plt.close()

        fig, axes = plt.subplots(1)
        ramp_min = min([np.min(ramp), np.min(model_ramp0), np.min(model_ramp1)])
        ramp_max = max([np.max(ramp), np.max(model_ramp0), np.max(model_ramp1)])
        for j, model_ramp in enumerate([model_ramp0, model_ramp1]):
            axes.plot(x_bins, ramp[j], label='Experiment'+str(j+1), color=colors[j+2])
            axes.plot(x_bins, model_ramp, label='Model'+str(j+1), color=colors[j])
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Depolarization (mV)')
        axes.set_xlim([0., track_length])
        axes.set_ylim([math.floor(ramp_min) - 0.5, math.ceil(ramp_max) + 0.5])
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes)
        plt.show()
        plt.close()

    formatted_x = '[' + ', '.join(['%.2f' % xi for xi in x]) + ']'
    print 'x: %s, Err: %.4E, Rule gain: %.4E' % (formatted_x, Err, rule_gain)
    if full_output:
        return [rule_waveform_offset0, rule_waveform_offset1], [this_weights0, this_weights1], \
               [model_ramp0, model_ramp1], new_baseline, rule_gain
    else:
        hist.x.append(x)
        hist.Err.append(Err)
        return Err


def optimize_polish(x, xmin, xmax, error_function, ramp, induction_index=None, orig_weights=None, baseline_indexes=None,
                    baseline=None, rule_gain=None, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param ramp: array
    :param induction_index: int: index of various location and time waves to use
    :param orig_weights: array
    :param baseline_indexes: list
    :param baseline: float
    :param rule_gain: float
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 400

    result = optimize.minimize(error_function, x, method='Nelder-Mead', options={'ftol': 1e-3,
                                                    'xtol': 1e-3, 'disp': True, 'maxiter': maxfev},
                               args=(xmin, xmax, ramp, induction_index, orig_weights, baseline_indexes,
                                     baseline, rule_gain))
    formatted_x = '['+', '.join(['%.2E' % xi for xi in result.x])+']'
    print 'Process: %i completed optimize_polish after %i iterations with Error: %.4E and x: %s' % (os.getpid(),
                                                                            result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


def optimize_explore(x, xmin, xmax, error_function, ramp, induction_index=None, orig_weights=None, baseline_indexes=None,
                     baseline=None, rule_gain=None, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param ramp: array
    :param induction_index: int: index of various location and time waves to use
    :param orig_weights: array
    :param baseline_indexes: list
    :param baseline: float
    :param rule_gain: float
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 700

    take_step = Normalized_Step(x, xmin, xmax)
    minimizer_kwargs = dict(method=null_minimizer, args=(xmin, xmax, ramp, induction_index, orig_weights,
                                                         baseline_indexes, baseline, rule_gain))
    result = optimize.basinhopping(error_function, x, niter=maxfev, niter_success=maxfev/2,
                                   disp=True, interval=min(20, int(maxfev/20)), minimizer_kwargs=minimizer_kwargs,
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
if cell_id == 1:
    x1 = [0.00, 80.72, 1.00, 4138.92, 1434.98]  # Err: 1.3314E+04, Rule gain: 1.3949E-01
elif cell_id == 2:
    x1 = [0.00, 64.48, 0.98, 6951.65, 100.04]  # Err: 1.8220E+05, Rule gain: 5.5148E-02
elif cell_id == 3:
    x1 = [0.00, 82.75, 0.25, 3562.13, 100.00]  # Err: 3.1113E+03, Rule gain: 1.9009E-01
elif cell_id == 4:
    x1 = [0.19, 68.51, 1.00, 6999.99, 6159.05]  # Err: 1.6502E+04, Rule gain: 1.3647E-02
elif cell_id == 5:
    x1 = [0.00, 75.58, 1.00, 2855.96, 932.52]  # Err: 1.1978E+03, Rule gain: 8.6899E-02
xmin1 = [0., 45., 0., 1000., 100.]
xmax1 = [1., 90., 1., 7000., 7000.]

# [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, depot_depth]
if cell_id == 1:
    # x2 = [0.00, 49.71, 1.00, 1749.80, 3998.71, 0.17]  # Err: 7.9991E+04, Rule gain: 8.3086E-02
    # x2 = [0.34, 71.20, 0.78, 1125.44, 3953.76, 0.73]
    x2 = [0.00, 81.86, 0.76, 1027.12, 7316.62, 0.57]  # Err: 4.6454E+04, Rule gain: 1.3950E-01
elif cell_id == 2:
    # x2 = [0.42, 81.10, 1.00, 2685.49, 2999.35, 0.5]
    x2 = [0.00, 78.83, 0.79, 1000.00, 2197.87, 0.23]
elif cell_id == 3:
    x2 = [0.00, 75.09, 0.06, 7000.00, 3996.26, 0.05]  # Err: 2.5270E+04, Rule gain: 1.3342E-01
elif cell_id == 4:
    x2 = [1.00, 84.44, 1.00, 4993.32, 4000.00, 0.47]  # Err: 7.5804E+04, Rule gain: 1.8681E-02
elif cell_id == 5:
    x2 = [1.00, 45.10, 0.39, 5999.77, 2571.92, 0.40]  # Err: 6.0807E+04, Rule gain: 9.1297E-02
xmin2 = [0., 45., 0., 100., 100., 0.]
xmax2 = [1., 90., 1., rule_max_timescale, rule_max_timescale, 1.]


# [during_depth, phase_offset, pot_amp, pre_pot_dur, post_pot_dur, depot_depth]
if cell_id == 1:
    # x2 = [0.00, 49.71, 1.00, 1749.80, 3998.71, 0.17]  # Err: 7.9991E+04, Rule gain: 8.3086E-02
    # x2 = [0.34, 71.20, 0.78, 1125.44, 3953.76, 0.73]
    # x3 = [0.00, 81.86, 0.76, 1027.12, 7316.62, 0.57, 20.]  # Err: 4.6454E+04, Rule gain: 1.3950E-01
    x3 = [0.40, 64.42, 0.56, 3430.22, 222.36, 0.01, 15.42]
elif cell_id == 2:
    # x2 = [0.42, 81.10, 1.00, 2685.49, 2999.35, 0.5]
    x3 = [0.00, 78.83, 0.79, 1000.00, 2197.87, 0.23, 20.]
elif cell_id == 3:
    x3 = [0.00, 75.09, 0.06, 7000.00, 3996.26, 0.05, 20.]  # Err: 2.5270E+04, Rule gain: 1.3342E-01
elif cell_id == 4:
    x3 = [1.00, 84.44, 1.00, 4993.32, 4000.00, 0.47, 20.]  # Err: 7.5804E+04, Rule gain: 1.8681E-02
elif cell_id == 5:
    x3 = [1.00, 45.10, 0.39, 5999.77, 2571.92, 0.40, 20.]  # Err: 6.0807E+04, Rule gain: 9.1297E-02
xmin3 = [0., 45., 0., 100., 100., 0., 10.]
xmax3 = [1., 90., 1., rule_max_timescale, rule_max_timescale, 2., 30.]

"""
induction_index = 0
result = optimize_explore(x0, xmin0, xmax0, ramp_error0, ramp[induction_index], induction_index,
                          baseline_indexes=ramp_baseline_indexes, maxfev=700)
polished_result = optimize_polish(result['x'], xmin0, xmax0, ramp_error0, ramp[induction_index], induction_index,
                                  baseline_indexes=ramp_baseline_indexes)
ramp_error0(polished_result['x'], xmin0, xmax0, ramp[induction_index], induction_index,
            baseline_indexes=ramp_baseline_indexes, plot=True)

induction_index = 0
this_rule_waveform0, this_raw_delta_weights0, this_model_ramp0, this_baseline0, this_rule_gain0 = \
    ramp_error0(x0, xmin0, xmax0, ramp[induction_index], induction_index, baseline_indexes=ramp_baseline_indexes,
                plot=False, full_output=True)
rule_waveforms.append(this_rule_waveform0)
delta_weights.append(this_raw_delta_weights0 * this_rule_gain0)
model_ramp.append(this_model_ramp0)


induction_index = 0
result = optimize_explore(x1, xmin1, xmax1, ramp_error1, ramp[induction_index], induction_index,
                          baseline_indexes=ramp_baseline_indexes, maxfev=700)
polished_result = optimize_polish(result['x'], xmin1, xmax1, ramp_error1, ramp[induction_index], induction_index,
                                  baseline_indexes=ramp_baseline_indexes)

ramp_error1(polished_result['x'], xmin1, xmax1, ramp[induction_index], induction_index,
            baseline_indexes=ramp_baseline_indexes, plot=True)


induction_index = 0
this_rule_waveform0, this_raw_delta_weights0, this_model_ramp0, this_baseline0, this_rule_gain0 = \
    ramp_error1(x1, xmin1, xmax1, ramp[induction_index], induction_index, baseline_indexes=ramp_baseline_indexes,
                plot=False, full_output=True)
rule_waveforms.append(this_rule_waveform0)
delta_weights.append(this_raw_delta_weights0 * this_rule_gain0)
model_ramp.append(this_model_ramp0)


induction_index = 1

result = optimize_explore(x2, xmin2, xmax2, ramp_error2,  ramp[induction_index], induction_index,
                          orig_weights=this_raw_delta_weights0, baseline=this_baseline0, rule_gain=this_rule_gain0,
                          maxfev=700)
polished_result = optimize_polish(result['x'], xmin2, xmax2, ramp_error2, ramp[induction_index], induction_index,
                                  orig_weights=this_raw_delta_weights0, baseline=this_baseline0,
                                  rule_gain=this_rule_gain0)
ramp_error2(polished_result['x'], xmin2, xmax2, ramp[induction_index], induction_index,
                orig_weights=this_raw_delta_weights0, baseline=this_baseline0, rule_gain=this_rule_gain0, plot=True)

this_rule_waveform1, this_raw_delta_weights1, this_model_ramp1, this_baseline1, this_rule_gain1 = \
    ramp_error2(x2, xmin2, xmax2, ramp[induction_index], induction_index, orig_weights=this_raw_delta_weights0,
                baseline=this_baseline0, rule_gain=this_rule_gain0, plot=False, full_output=True)

rule_waveforms.append(this_rule_waveform1)
delta_weights.append(this_raw_delta_weights1 * this_rule_gain1)
model_ramp.append(this_model_ramp1)
"""

# result = optimize_explore(x3, xmin3, xmax3, ramp_error3, ramp, baseline_indexes=ramp_baseline_indexes, maxfev=1000)
"""
polished_result = optimize_polish(result['x'], xmin3, xmax3, ramp_error3, ramp, baseline_indexes=ramp_baseline_indexes)
ramp_error3(polished_result['x'], xmin3, xmax3, ramp, baseline_indexes=ramp_baseline_indexes, plot=True)

rule_waveforms, delta_weights, model_ramp, new_baseline, rule_gain = ramp_error3(x3, xmin3, xmax3, ramp, plot=False,
                                                                                 full_output=True)

# colors = ['r', 'k', 'c', 'g']
colors = ['k', 'r', 'g', 'c']
x_start = [induction_loc/track_length for induction_loc in induction_locs]

ylim = max(np.max(ramp[0]), np.max(model_ramp))
#ylim = max(np.max(ramp), np.max(model_ramp))
fig, axes = plt.subplots(1)
for i in range(len(model_ramp)):
    axes.plot(x_bins, ramp[i], color=colors[i+2], label='Experiment'+str(i+1))
    axes.plot(x_bins, model_ramp[i], color=colors[i], label='Model'+str(i+1))
    axes.axhline(y=ylim+0.3, xmin=x_start[i], xmax=x_start[i]+0.02, c=colors[i], linewidth=3., zorder=0)
axes.set_xlabel('Location (cm)')
axes.set_xlim(0., track_length)
axes.set_ylabel('Depolarization (mV)')
axes.legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes)
plt.show()
plt.close()

fig, axes = plt.subplots(1)
t_waveform = np.arange(-rule_max_timescale, rule_max_timescale, dt)
for i in range(len(model_ramp)):
    axes.plot(t_waveform, rule_waveforms[i], label='Induction kernel '+str(i+1), color=colors[i])
axes.set_xlabel('Time (ms)')
axes.set_xlim(-rule_max_timescale, rule_max_timescale)
axes.set_ylabel('Relative change in synaptic weight per spike')
axes.legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes)
plt.show()
plt.close()

ylim = np.max(delta_weights) + 1.
fig, axes = plt.subplots(1)
for i in range(len(model_ramp)):
    axes.scatter(peak_locs['CA3'], delta_weights[i] + 1., label='Weights'+str(i+1), color=colors[i])
    axes.axhline(y=ylim + 0.1, xmin=x_start[i], xmax=x_start[i] + 0.01, c=colors[i], linewidth=3., zorder=0)
axes.set_xlabel('Location (cm)')
axes.set_xlim(0., track_length)
axes.set_ylabel('Relative synaptic weight')
axes.legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes)
plt.show()
plt.close()

if len(delta_weights) > 1:
    fig, axes = plt.subplots(1)
    axes.scatter(delta_weights[0]+1., delta_weights[1]+1.)
    clean_axes(axes)
    plt.show()
    plt.close()
"""