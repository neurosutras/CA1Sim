__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import random
import sys
import mkl

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
weights_filename = '121316 induced field long rule weights'

if len(sys.argv) > 1:
    induction_run_vel = float(sys.argv[1])
else:
    induction_run_vel = 30.
# allows parallel computation of multiple trials for the same spines with the same peak_locs, but with different
# input spike trains and stochastic synapses for each trial
if len(sys.argv) > 2:
    induction_seed = int(sys.argv[2])
else:
    induction_seed = 0

mkl.set_num_threads(8)


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
        formatted_x = '[' + ', '.join(['%.3E' % xi for xi in best_x]) + ']'
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


def calculate_ramp_features(ramp, induction_loc):
    """

    :param x: array
    :param ramp: array
    :param induction_loc: float
    """
    extended_binned_x = np.concatenate([binned_x - track_length, binned_x, binned_x + track_length])
    extended_binned_ramp = np.concatenate([ramp for i in range(3)])
    extended_interp_x = np.concatenate([default_interp_x - track_length, default_interp_x,
                                        default_interp_x + track_length])
    dx = extended_interp_x[1] - extended_interp_x[0]
    extended_ramp = np.interp(extended_interp_x, extended_binned_x, extended_binned_ramp)
    interp_ramp = extended_ramp[len(default_interp_x):2*len(default_interp_x)]
    baseline_indexes = np.where(interp_ramp <= np.percentile(interp_ramp, 10.))[0]
    baseline = np.mean(interp_ramp[baseline_indexes])
    interp_ramp -= baseline
    extended_ramp -= baseline
    peak_index = np.where(interp_ramp == np.max(interp_ramp))[0][0] + len(interp_ramp)
    # use center of mass in 10 spatial bins instead of literal peak for determining peak_shift
    before_peak_index = peak_index-int(track_length/10./2./dx)
    after_peak_index = peak_index + int(track_length/10./2./dx)
    area_around_peak = np.trapz(extended_ramp[before_peak_index:after_peak_index], dx=dx)
    for i in range(before_peak_index+1, after_peak_index):
        this_area = np.trapz(extended_ramp[before_peak_index:i], dx=dx)
        if this_area/area_around_peak >= 0.5:
            center_of_mass_index = i
            break
    center_of_mass_val = np.mean(extended_ramp[before_peak_index:after_peak_index])

    if extended_interp_x[center_of_mass_index] > induction_loc + 30.:
        center_of_mass_index -= len(interp_ramp)
    center_of_mass_x = extended_interp_x[center_of_mass_index]
    start_index = np.where(extended_ramp[:center_of_mass_index] <= 0.15*center_of_mass_val)[0][-1]
    end_index = center_of_mass_index + np.where(extended_ramp[center_of_mass_index:] <= 0.15*center_of_mass_val)[0][0]
    peak_shift = center_of_mass_x - induction_loc
    ramp_width = extended_interp_x[end_index] - extended_interp_x[start_index]
    before_width = induction_loc - extended_interp_x[start_index]
    after_width = extended_interp_x[end_index] - induction_loc
    ratio = before_width / after_width
    return center_of_mass_val, ramp_width, peak_shift, ratio


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


def wrap_around_and_compress(waveform, interp_x):
    """

    :param waveform: array of len(3 * interp_x)
    :param interp_x: array
    :return: array of len(interp_x)
    """
    before = np.array(waveform[:len(interp_x)])
    after = np.array(waveform[2 * len(interp_x):])
    within = np.array(waveform[len(interp_x):2 * len(interp_x)])
    waveform = within[:len(interp_x)] + before[:len(interp_x)] + after[:len(interp_x)]
    
    return waveform


def get_dynamic_theta_phase_force(phase_ranges, peak_loc, extended_x, interp_x, dx):
    """
    Expects a list of tuples containing times and phases relative to peak_loc and the non-modulated phase preference
    (zero degrees). Returns a waveform of phase vs time.
    :param phase_ranges: list of tuple (ms, degrees)
    :param peak_loc: float (cm)
    :param extended_x: array (just one item from the outer list)
    :param interp_x: array (just one item from the outer list)
    :param dx: float (cm)
    :return: :class: 'np.array'
    """
    end_x_val = phase_ranges[-1][0] + peak_loc
    start_phase_val = phase_ranges[0][1] * 2. * np.pi / 360.  # convert degrees to radians
    end_phase_val = phase_ranges[-1][1] * 2. * np.pi / 360.  # convert degrees to radians
    phase_force = np.ones_like(extended_x) * start_phase_val
    phase_force[np.where(extended_x >= end_x_val)[0]] = end_phase_val
    for i in range(len(phase_ranges) - 1):
        x0 = phase_ranges[i][0] + peak_loc
        x1 = phase_ranges[i + 1][0] + peak_loc
        phase0 = phase_ranges[i][1] * 2. * np.pi / 360.  # convert degrees to radians
        phase1 = phase_ranges[i + 1][1] * 2. * np.pi / 360.  # convert degrees to radians
        del_x = x1 - x0
        del_phase = phase1 - phase0
        if abs(del_x) > 0.:
            this_x_piece = np.arange(x0, x1 + dx / 2., dx)
            this_indexes = np.where((extended_x >= x0) & (extended_x < x1))[0]
            this_interp_x_piece = extended_x[this_indexes]
            if abs(del_phase) > 0.:
                d_phase = del_phase / del_x * dx
                this_range_piece = np.arange(phase0, phase1 + d_phase / 2., d_phase)
                this_range_interp_piece = np.interp(this_interp_x_piece, this_x_piece, this_range_piece)
            else:
                this_range_interp_piece = np.ones(len(this_interp_x_piece)) * phase0
            phase_force[this_indexes] = this_range_interp_piece
    phase_force = wrap_around_and_compress(phase_force, interp_x)
    return phase_force


def generate_generic_rate_and_phase_maps():
    """

    """
    spatial_rate_maps = {'CA3': [], 'ECIII': []}
    phase_maps = {'CA3': []}
    for group in stim_exc_syns:
        for i, syn in enumerate(stim_exc_syns[group]):
            peak_loc = peak_locs[group][i]
            gauss_force = excitatory_peak_rate[group] * np.exp(-((extended_x - peak_loc) / gauss_sigma)**2.)
            gauss_force = wrap_around_and_compress(gauss_force, generic_x)
            spatial_rate_maps[group].append(gauss_force)
            if group in excitatory_precession_range:
                phase_force = get_dynamic_theta_phase_force(excitatory_precession_range[group], peak_loc, extended_x,
                                                            generic_x, generic_dx)
                phase_maps[group].append(phase_force)
            
    return spatial_rate_maps, phase_maps


def generate_complete_rate_maps(induction, spatial_rate_maps, phase_maps, global_phase_offset=None):
    """

    :param induction: int
    :param spatial_rate_maps: array
    :param phase_maps: array
    :param global_phase_offset: float
    """
    if global_phase_offset is None:
        global_phase_offset = local_random.uniform(-np.pi, np.pi)
    complete_rate_maps = {'CA3': [], 'ECIII': []}
    running_delta_t = 0.
    for i in range(len(interp_x[induction])):
        this_interp_t = interp_t[induction][i]
        this_interp_x = interp_x[induction][i]
        this_run_vel_gate = run_vel_gate[induction][i]
        if i == 0:
            complete_t = this_interp_t
            complete_x = this_interp_x - track_length
            complete_run_vel_gate = np.array(this_run_vel_gate)
            running_delta_t += len(this_interp_t) * dt
        complete_t = np.append(complete_t, this_interp_t + running_delta_t)
        complete_x = np.append(complete_x, this_interp_x + i * track_length)
        complete_run_vel_gate = np.append(complete_run_vel_gate, this_run_vel_gate)
        running_delta_t += len(this_interp_t) * dt
        if i == len(interp_x[induction]) - 1:
            complete_t = np.append(complete_t, this_interp_t + running_delta_t)
            complete_x = np.append(complete_x, this_interp_x + (i+1) * track_length)
            complete_run_vel_gate = np.append(complete_run_vel_gate, this_run_vel_gate)
    for group in peak_locs:
        for j in range(len(spatial_rate_maps[group])):
            for i in range(len(interp_x[induction])):
                this_interp_x = interp_x[induction][i]
                this_rate_map = np.interp(this_interp_x, generic_x, spatial_rate_maps[group][j])
                if i == 0:
                    this_complete_rate_map = np.array(this_rate_map)
                this_complete_rate_map = np.append(this_complete_rate_map, this_rate_map)
                if i == len(interp_x[induction]) - 1:
                    this_complete_rate_map = np.append(this_complete_rate_map, this_rate_map)
                if group in excitatory_precession_range:
                    this_phase_map = np.interp(this_interp_x, generic_x, phase_maps[group][j])
                    if i == 0:
                        this_complete_phase_map = np.array(this_phase_map)
                    this_complete_phase_map = np.append(this_complete_phase_map, this_phase_map)
                    if i == len(interp_x[induction]) - 1:
                        this_complete_phase_map = np.append(this_complete_phase_map, this_phase_map)
            if group in excitatory_precession_range:
                theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] *
                                     np.cos(this_complete_phase_map + excitatory_theta_phase_offset[group] -
                                            2. * np.pi * complete_t / global_theta_cycle_duration +
                                            global_phase_offset))
            else:
                theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] *
                                     np.cos(excitatory_theta_phase_offset[group] - 2. * np.pi * complete_t /
                                            global_theta_cycle_duration + global_phase_offset))
            theta_force -= np.min(theta_force)
            theta_force /= np.max(theta_force)
            theta_force *= excitatory_theta_modulation_depth[group]
            theta_force += 1. - excitatory_theta_modulation_depth[group]
            this_complete_rate_map = np.multiply(this_complete_rate_map, theta_force)
            this_complete_rate_map = np.multiply(this_complete_rate_map, complete_run_vel_gate)
            complete_rate_maps[group].append(this_complete_rate_map)

    return complete_t, complete_x, complete_rate_maps


def generate_default_rate_maps(spatial_rate_maps, phase_maps, global_phase_offset=None):
    """

    :param spatial_rate_maps: array
    :param phase_maps: array
    :param global_phase_offset: float
    """
    if global_phase_offset is None:
        global_phase_offset = local_random.uniform(-np.pi, np.pi)
    default_rate_maps = {'CA3': [], 'ECIII': []}
    for group in peak_locs:
        for j in range(len(spatial_rate_maps[group])):
            this_rate_map = np.interp(default_interp_x, generic_x, spatial_rate_maps[group][j])
            if group in excitatory_precession_range:
                this_phase_map = np.interp(default_interp_x, generic_x, phase_maps[group][j])
                theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] *
                                     np.cos(this_phase_map + excitatory_theta_phase_offset[group] -
                                            2. * np.pi * default_interp_t / global_theta_cycle_duration +
                                            global_phase_offset))
            else:
                theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] *
                                     np.cos(excitatory_theta_phase_offset[group] -
                                            2. * np.pi * default_interp_t / global_theta_cycle_duration +
                                            global_phase_offset))
            theta_force -= np.min(theta_force)
            theta_force /= np.max(theta_force)
            theta_force *= excitatory_theta_modulation_depth[group]
            theta_force += 1. - excitatory_theta_modulation_depth[group]
            this_rate_map = np.multiply(this_rate_map, theta_force)
            default_rate_maps[group].append(this_rate_map)

    return default_rate_maps


def generate_complete_induction_gate(induction):
    """

    :param induction: int
    :return:
    """
    for i in range(len(interp_x[induction])):
        this_interp_t = interp_t[induction][i]
        this_interp_x = interp_x[induction][i]
        if i == 0:
            complete_induction_gate = np.zeros_like(this_interp_t)
        this_induction_loc = induction_locs[induction][i]
        this_induction_dur = induction_durs[induction][i]
        this_induction_gate = np.zeros_like(this_interp_t)
        start_index = np.where(this_interp_x >= this_induction_loc)[0][0]
        end_index = start_index + int(this_induction_dur / dt)
        this_induction_gate[start_index:end_index] = 1.
        complete_induction_gate = np.append(complete_induction_gate, this_induction_gate)
        if i == len(interp_x[induction]) - 1:
            complete_induction_gate = np.append(complete_induction_gate, np.zeros_like(this_interp_t))
    return complete_induction_gate


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


def get_expected_depolarization(rate_maps, weights, interp_x):
    """
    Take pre-computed rate maps and weights for a set of inputs. Convolve with an EPSP kernel, sum, normalize, and
    downsample to 100 spatial bins.
    :param rate_maps: list of array
    :param weights: list of float
    :param interp_x: array (just one item from the outer list)
    :return: array
    """
    filter_t = np.arange(0., 200., dt)
    epsp_filter = np.exp(-filter_t/20.) - np.exp(-filter_t/.5)
    epsp_filter /= np.sum(epsp_filter)
    weighted_rate_maps = []
    scaling_factor = 7.15e-4  # generates a predicted 6 mV depolarization from gaussian weights with peak = 2.5
    for group in peak_locs:
        for i, rate_map in enumerate(rate_maps[group]):
            this_weighted_rate_map = rate_map * weights[group][i] * scaling_factor
            this_weighted_rate_map = np.concatenate([this_weighted_rate_map for i in range(3)])
            this_weighted_rate_map = np.convolve(this_weighted_rate_map, epsp_filter)[:3*len(interp_x)]
            weighted_rate_maps.append(this_weighted_rate_map)
    expected_depolarization = np.sum(weighted_rate_maps, axis=0)
    expected_depolarization = low_pass_filter(expected_depolarization, 2.,
                                              3*len(interp_x)*dt, dt, 1.)[len(interp_x): 2 * len(interp_x)]
    expected_depolarization = np.interp(binned_x, interp_x, expected_depolarization)
    return expected_depolarization


def generate_spike_trains(rate_maps, t):
    """

    :param rate_maps: list of array
    :return: list of array
    """
    spike_trains = {}
    for group in rate_maps:
        spike_trains[group] = []
        for rate_map in rate_maps[group]:
            this_spike_train = get_inhom_poisson_spike_times_by_thinning(rate_map, t, dt=0.02, generator=local_random)
            spike_trains[group].append(this_spike_train)
    return spike_trains


def filter_spike_trains(spike_trains):
    """

    :param spike_trains: list of array
    :return: list of array
    """
    dynamics = [0.2, 1.769, 67.351, 0.878, 92.918]
    successes = {}
    for group in spike_trains:
        successes[group] = []
        for spike_train in spike_trains[group]:
            this_Pr = Pr(*dynamics)
            this_success_train = []
            for spike_time in spike_train:
                P = this_Pr.stim(spike_time)
                if local_random.random() < P:
                    this_success_train.append(spike_time)
            successes[group].append(this_success_train)
    return successes



NMDA_type = 'NMDA_KIN5'

dt = 1.  # ms
down_dt = 10.  # ms, to speed up optimization
equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # -0.5 * run_vel + 155.  # (ms)
input_field_width = 90.  # cm
track_length = 187.  # cm
field1_loc = 0.5

binned_dx = track_length / 100.  # cm
binned_x = np.arange(0., track_length+binned_dx/2., binned_dx)[:100] + binned_dx/2.
generic_dx = binned_dx / 100.  # cm
generic_x = np.arange(0., track_length, generic_dx)

default_run_vel = 30.  # cm/s
generic_position_dt = generic_dx / induction_run_vel * 1000.
generic_t = np.arange(0., len(generic_x)*generic_position_dt, generic_position_dt)[:len(generic_x)]
generic_track_duration = len(generic_t) * generic_position_dt

default_interp_t = np.arange(0., generic_t[-1], dt)
default_interp_x = np.interp(default_interp_t, generic_t, generic_x)

extended_x = np.concatenate([generic_x - track_length, generic_x, generic_x + track_length])

induction_locs = {}
induction_durs = {}
interp_t = {}
interp_x = {}
run_vel = {}
run_vel_gate = {}

vel_window_dur = 10.  # ms
vel_window_bins = int(vel_window_dur/dt)/2

pause_duration = 3000.

for induction in [1]:
    interp_x[induction] = [np.append(np.zeros_like(np.arange(0., pause_duration, dt)), default_interp_x) for i in range(5)]
    induction_locs[induction] = [field1_loc*track_length for i in range(5)]
    induction_durs[induction] = [300. for i in range(5)]
    interp_t[induction] = [np.append(np.arange(0., pause_duration, dt), default_interp_t + pause_duration) for i in range(5)]
    run_vel[induction] = []
    run_vel_gate[induction] = []
    for i in range(5):
        this_interp_x = interp_x[induction][i]
        this_interp_t = interp_t[induction][i]
        padded_t = np.insert(this_interp_t, 0, this_interp_t[-vel_window_bins:]-generic_track_duration-pause_duration)
        padded_t = np.append(padded_t, this_interp_t[:vel_window_bins + 1] + generic_track_duration+pause_duration)
        padded_x = np.insert(this_interp_x, 0, this_interp_x[-vel_window_bins:]-track_length)
        padded_x = np.append(padded_x, this_interp_x[:vel_window_bins+1]+track_length)
        this_run_vel = []
        for j in range(vel_window_bins, vel_window_bins+len(this_interp_x)):
            this_run_vel.append(np.sum(np.diff(padded_x[j-vel_window_bins:j+vel_window_bins+1]))/
                                np.sum(np.diff(padded_t[j - vel_window_bins:j + vel_window_bins + 1]))*1000.)
        this_run_vel = np.array(this_run_vel)
        this_run_vel_gate = np.zeros_like(this_run_vel)
        this_run_vel_gate[np.where(this_run_vel>5.)[0]] = 1.
        run_vel[induction].append(this_run_vel)
        run_vel_gate[induction].append(this_run_vel_gate)


track_equilibrate = 2. * global_theta_cycle_duration
track_duration = {induction:
                      [interp_t[induction][i][-1] for i in range(len(interp_t[induction]))] for induction in interp_t}

excitatory_peak_rate = {'CA3': 40., 'ECIII': 40.}
excitatory_theta_modulation_depth = {'CA3': 0.7, 'ECIII': 0.7}
# From Chadwick et al., ELife 2015
excitatory_theta_phase_tuning_factor = {'CA3': 0.8, 'ECIII': 0.8}
excitatory_precession_range = {}  # (ms, degrees)
excitatory_precession_range['CA3'] = [(-input_field_width*0.7, 0.), (-input_field_width*0.6, 180.),
                                      (-input_field_width*0.35, 180.), (input_field_width*0.35, -180.),
                                      (input_field_width*0.6, -180.), (input_field_width*0.7, 0.)]
excitatory_theta_phase_offset = {}
excitatory_theta_phase_offset['CA3'] = 165. / 360. * 2. * np.pi  # radians
excitatory_theta_phase_offset['ECIII'] = 0. / 360. * 2. * np.pi  # radians
excitatory_stochastic = 1

v_init = -67.

syn_types = ['AMPA_KIN', NMDA_type]

local_random = random.Random()

# choose a subset of synapses to stimulate with inhomogeneous poisson rates
local_random.seed(synapses_seed)

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

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
stim_exc_syns = {'CA3': [], 'ECIII': []}
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
# cell.init_synaptic_mechanisms()

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

gauss_sigma = input_field_width / 3. / np.sqrt(2.)  # contains 99.7% gaussian area

delta_peak_locs = {}
extended_peak_locs = {}
within_track = {}
for group in stim_exc_syns:
    delta_peak_locs[group] = track_length / len(stim_exc_syns[group])
    if stim_exc_syns[group]:
        extended_peak_locs[group] = np.arange(-track_length, 2. * track_length, delta_peak_locs[group])
        within_track[group] = np.where((extended_peak_locs[group] >= 0.) &
                                       (extended_peak_locs[group] < track_length))[0]
        peak_locs[group] = extended_peak_locs[group][within_track[group]]

for group in stim_exc_syns:
    for syn in stim_exc_syns[group]:
        all_exc_syns[syn.node.parent.parent.type].remove(syn)

default_global_phase_offset = 0.

spatial_rate_maps, phase_maps = generate_generic_rate_and_phase_maps()

local_random.seed(induction_seed)

complete_t, complete_x, complete_rate_maps, complete_induction_gates = {}, {}, {}, {}
for induction in [1]:
    complete_t[induction], complete_x[induction], complete_rate_maps[induction] = \
        generate_complete_rate_maps(induction, spatial_rate_maps, phase_maps)
    complete_induction_gates[induction] = generate_complete_induction_gate(induction)

default_rate_maps = generate_default_rate_maps(spatial_rate_maps, phase_maps, default_global_phase_offset)

induction = 1
spike_trains = generate_spike_trains(complete_rate_maps[induction], complete_t[induction])
successes = filter_spike_trains(spike_trains)


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


def build_kernels(x, plot=False):
    """
    Construct two kernels with exponential rise and decay:
    1) Local kernel that generates a plasticity signal at each spine
    2) Global kernal that generates a plasticity signal during dendritic calcium spiking
    :param x: array: [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio, global_scale]
    :param plot: bool
    :return: array, array
    """
    local_rise_tau = x[0]
    local_decay_tau = x[1]
    global_rise_tau = x[2]
    global_decay_tau = x[3]
    filter_ratio = x[4]

    max_time_scale = np.max([local_rise_tau+local_decay_tau, global_rise_tau+global_decay_tau])
    filter_t = np.arange(0., 6.*max_time_scale, dt)
    local_filter = np.exp(-filter_t/local_decay_tau) - np.exp(-filter_t/local_rise_tau)
    peak_index = np.where(local_filter == np.max(local_filter))[0][0]
    decay_indexes = np.where(local_filter[peak_index:] < 0.005*np.max(local_filter))[0]
    if np.any(decay_indexes):
        local_filter = local_filter[:peak_index+decay_indexes[0]]
    local_filter /= np.sum(local_filter)

    global_filter = np.exp(-filter_t / global_decay_tau) - np.exp(-filter_t / global_rise_tau)
    peak_index = np.where(global_filter == np.max(global_filter))[0][0]
    decay_indexes = np.where(global_filter[peak_index:] < 0.005 * np.max(global_filter))[0]
    if np.any(decay_indexes):
        global_filter = global_filter[:peak_index + decay_indexes[0]]
    global_filter /= np.sum(global_filter)

    if plot:
        fig, axes = plt.subplots(1)
        axes.plot(filter_t[:len(local_filter)], local_filter / np.max(local_filter), color='g',
                  label='Local plasticity kernel')
        axes.plot(filter_t[:len(global_filter)], global_filter / np.max(global_filter) * filter_ratio, color='k',
                  label='Global plasticity kernel')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Time (ms)')
        axes.set_ylabel('Relative kernel amplitude (a.u.)')
        axes.set_xlim(-500., max(5000., max(filter_t[:len(local_filter)][-1], filter_t[:len(global_filter)][-1])))
        clean_axes(axes)
        plt.show()
        plt.close()

    return local_filter, global_filter


def calculate_plasticity_signal_discrete(x, local_kernel, global_kernel, induction, plot=False,
                                         use_global_saturation_factor=False):
    """
    Given the local and global kernels, convolve each input spike train (successes) with the local kernel, and convolve
    the current injection with the global kernel. The weight change for each input is proportional to the area under the
    product of the two signals. Incremental weight changes accrue across multiple induction trials.
    :param x: array: [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio, kernel_scale]
    :param local_kernel: array
    :param global_kernel: array
    :param induction: int: key for dicts of arrays
    :param plot: bool
    :param use_global_saturation_factor: bool  : calculate locally during optimization. use one global value for all
                                                running speeds for simulations (use 30. cm/s and peak_weights=3.
                                                as a reference)
    :return: plasticity_signal: array
    """
    if use_global_saturation_factor:
        saturation_factor = global_saturation_factor
        attempt_range = range(1,2)
    else:
        saturation_factor = 0.02
        attempt_range = range(2)
    filter_ratio = x[4]
    # kernel_scale = x[5]
    global_signal = np.convolve(complete_induction_gates[induction], global_kernel)[:len(complete_t[induction])] * \
                    kernel_scale
    max_global_signal = np.max(global_signal)
    filter_t = np.arange(0., len(local_kernel) * dt, dt)
    plasticity_signal = {}
    for attempt in attempt_range:
        local_signal_array = []
        start_time = time.time()
        for group in peak_locs:
            plasticity_signal[group] = np.zeros_like(peak_locs[group])
            for j, train in enumerate(successes[group]):
                indexes = (np.array(train)/dt).astype(int)
                this_stim_force = np.zeros_like(complete_t[induction])
                this_stim_force[indexes] = 1.
                local_signal = np.convolve(this_stim_force, local_kernel)[:len(complete_t[induction])] / \
                               saturation_factor * kernel_scale / filter_ratio
                if group == 'CA3':
                    local_signal_array.append(np.max(local_signal))
                this_signal = np.minimum(local_signal, global_signal)
                this_area = np.trapz(this_signal, dx=dt)
                plasticity_signal[group][j] += this_area
                if plot and j == int(len(successes[group])/2) and group == 'CA3' and attempt == 1:
                    ylim = max(np.max(local_signal), max_global_signal)
                    start_index = np.where(interp_x[induction][0] >= induction_locs[induction][0])[0][0]
                    this_induction_start = interp_t[induction][0][start_index]
                    this_induction_dur = induction_durs[induction][0]
                    start_time = -5000.
                    end_time = interp_t[induction][0][-1] + 5000.
                    this_duration = end_time - start_time
                    x_start = (5000. + this_induction_start) / this_duration
                    x_end = (5000. + this_induction_start + this_induction_dur) / this_duration
                    fig, axes = plt.subplots(1)
                    axes_t = np.array((complete_t[induction] - len(interp_t[induction][0])*dt)/1000.)
                    axes.plot(axes_t, local_signal, label='Local signal', color='g')
                    axes.plot(axes_t, global_signal, label='Global signal', color='k')
                    axes.fill_between(axes_t, 0., this_signal, label='Overlap', facecolor='r', alpha=0.5)
                    axes.axhline(y=ylim*1.05, xmin=x_start, xmax=x_end, linewidth=3, c='k')
                    axes.legend(loc='best', frameon=False, framealpha=0.5)
                    axes.set_xlabel('Time (s)')
                    axes.set_ylabel('Signal amplitude (a.u.)')
                    axes.set_xlim(-5., interp_t[induction][0][-1]/1000. + 5.)
                    axes.set_ylim(-0.05*ylim, ylim*1.1)
                    axes.set_title('Induced plasticity signal')
                    clean_axes(axes)
                    plt.show()
                    plt.close()
        if attempt == 0:
            local_signal_array = np.array(local_signal_array)
            max_local_signal = np.mean(local_signal_array[np.where(local_signal_array >=
                                                                   np.percentile(local_signal_array, 90.))[0]])
            saturation_factor *= filter_ratio * max_local_signal / max_global_signal
        # print 'Computed weights in %i s' % (time.time() - start_time)
        print 'saturation_factor after attempt %i: %.3E' % (attempt, saturation_factor)

    if plot:
        x_start = np.mean(induction_locs[induction]) / track_length
        ylim = np.max(plasticity_signal['CA3']) + 1.
        fig, axes = plt.subplots(1)
        axes.scatter(peak_locs['CA3'], plasticity_signal['CA3']+1., color='k')
        axes.axhline(y=ylim + 0.1, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='k')
        axes.set_ylabel('Relative synaptic weight')
        axes.set_xlabel('Location (cm)')
        axes.set_xlim(0., track_length)
        axes.set_title('Induced synaptic weight distribution')
        clean_axes(axes)
        plt.show()
        plt.close()

    return plasticity_signal


def calculate_weights_discrete(x, induction=None, plot=False, full_output=False, use_global_saturation_factor=False):
    """
    Calculates a rule_waveform and set of weights to match the first place field induction from trains of spike times
    (successes).
    :param x: array [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio, kernel_scale]
    :param induction: int: key for dicts of arrays
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :param use_global_saturation_factor: bool
    :return: tuple of array
    """
    start_time = time.time()
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    print 'Trying x: %s' % formatted_x
    if induction is None:
        induction = 1
    local_kernel, global_kernel = build_kernels(x, plot)
    this_weights = calculate_plasticity_signal_discrete(x, local_kernel, global_kernel, induction, plot,
                                                        use_global_saturation_factor)
    model_ramp = get_expected_depolarization(default_rate_maps,
                                             {group: this_weights[group] + 1. for group in this_weights},
                                             default_interp_x)
    baseline_indexes = np.where(np.array(model_ramp) <= np.percentile(model_ramp, 10.))[0]
    model_ramp_baseline = np.mean(model_ramp[baseline_indexes])
    model_ramp -= model_ramp_baseline

    if plot:
        x_start = np.mean(induction_locs[induction])/track_length
        ylim = np.max(model_ramp)
        ymin = np.min(model_ramp)
        fig, axes = plt.subplots(1)
        axes.plot(binned_x, model_ramp, label='Model', color='k')
        axes.axhline(y=ylim + 0.1, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='r')
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Depolarization (mV)')
        axes.set_xlim([0., track_length])
        axes.set_ylim([math.floor(ymin) - 0.5, math.ceil(ylim) + 0.5])
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_title('Induced Vm depolarization')
        clean_axes(axes)
        plt.show()
        plt.close()

    if full_output:
        return local_kernel, global_kernel, this_weights, model_ramp


def ramp_error_discrete(x, xmin, xmax, induction=None, plot=False, full_output=False,
                        use_global_saturation_factor=False):
    """
    Calculates a rule_waveform and set of weights to match the first place field induction.
    :param x: array [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio, kernel_scale]
    :param xmin: array
    :param xmax: array
    :param induction: int: key for dicts of arrays
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :param use_global_saturation_factor: bool
    :return: float
    """
    if not check_bounds(x, xmin, xmax):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    elif (x[1] <= x[0]) or (x[3] <= x[2]):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    start_time = time.time()
    if induction is None:
        induction = 1
    local_kernel, global_kernel, this_weights, ramp = calculate_weights_discrete(x, induction, plot,
                                                                                 full_output=True,
                                                                                 use_global_saturation_factor=
                                                                                 use_global_saturation_factor)
    this_induction_loc = np.mean(induction_locs[induction])
    amp, width, shift, ratio = calculate_ramp_features(ramp, this_induction_loc)
    result = {'ramp_width': width,
              'peak_shift': shift,
              'ratio': ratio,
              'min_weights': np.min(this_weights['CA3'])}

    Err = 0.
    for feature in target_vals:
        Err_piece = ((target_vals[feature] - result[feature]) / target_var[feature]) ** 2.
        Err += Err_piece

    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    print 'x: %s, Err: %.4E took %i s' % (formatted_x, Err, time.time() - start_time)
    # print 'ramp_width: %.1f, peak_weights: %.1f' % (result['ramp_width'], result['peak_weights'])
    print 'ramp_width: %.1f, peak_shift: %.1f, ratio: %.1f' % (result['ramp_width'], result['peak_shift'],
                                                               result['ratio'])
    if full_output:
        return local_kernel, global_kernel, this_weights, ramp
    else:
        hist.x.append(x)
        hist.Err.append(Err)
        return Err


def optimize_polish(x, xmin, xmax, error_function, induction=None, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param induction: int: key for dicts of arrays
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 400

    result = optimize.minimize(error_function, x, method='Nelder-Mead', options={'ftol': 1e-3,
                                                                                 'xtol': 1e-3, 'disp': True,
                                                                                 'maxiter': maxfev},
                               args=(xmin, xmax, induction))
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in result.x]) + ']'
    print 'Process: %i completed optimize_polish after %i iterations with Error: %.4E and x: %s' % (os.getpid(),
                                                                                                    result.nit,
                                                                                                    result.fun,
                                                                                                    formatted_x)
    return {'x': result.x, 'Err': result.fun}


def optimize_explore(x, xmin, xmax, error_function, induction=None, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param induction: int: key for dicts of arrays
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 700

    take_step = Normalized_Step(x, xmin, xmax)
    minimizer_kwargs = dict(method=null_minimizer, args=(xmin, xmax, induction))
    result = optimize.basinhopping(error_function, x, niter=maxfev, niter_success=maxfev / 2,
                                   disp=True, interval=min(20, int(maxfev / 20)), minimizer_kwargs=minimizer_kwargs,
                                   take_step=take_step)
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in result.x]) + ']'
    print 'Process: %i completed optimize_explore after %i iterations with Error: %.4E and x: %s' % (os.getpid(),
                                                                                                     result.nit,
                                                                                                     result.fun,
                                                                                                     formatted_x)
    return {'x': result.x, 'Err': result.fun}


# x1 = [2.795E+02, 1.236E+03, 2.268E+01, 4.522E+02, 1.397E+00, 2.763E-03]  # just induced
# kernel_scale corrected to produce closer to 2.5 mean weight in center 10 spatial bins 6 mV with 30 cm/s during
# induction laps
# x1 = [2.614E+02, 1.103E+03, 2.686E+01, 4.574E+02, 1.388E+00, 4.533E-03*0.4619]  # induced + spontaneous
# x1 = [2.066E+02, 1.033E+03, 2.861E+01, 4.317E+02, 1.205E+00, 3.708E-03]  # induced + spontaneous 121316
# x1 = [4.136E+02, 1.870E+03, 4.407E+01, 3.376E+02, 8.294E-01, 5.298E-04]
# x1 = [3.016E+02, 1.166E+03, 4.640E+01, 4.363E+02, 9.224E-01]
x1 = [1.002E+01, 1.996E+03, 3.104E+01, 6.577E+02, 8.169E-01]  # Error: 8.2169E+01

xmin = [10., 300., 10., 25., 0.7]  # , 5.e-4]
xmax = [500., 5000., 50., 2000., 1.5]  # , 2.e-2]

kernel_scale = 2.5e-3 * 0.9804 # 121416 target weights: 3.
global_saturation_factor = 2.747E-02

target_vals = {'ramp_width': 140., 'peak_shift': -10., 'ratio': 1.7, 'min_weights': 0.}
target_var = {'ramp_width': 0.1, 'peak_shift': 0.05, 'ratio': 0.05, 'min_weights': 0.05}

# to avoid saturation and reduce variability of time courses across cells, constrain the relative amplitude
# of global and local kernels:
# [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio, kernel_scale]

induction = 1

"""
result = optimize_polish(x1, xmin, xmax, ramp_error_discrete, induction)
hist.report_best()
hist.export('121416_magee_data_optimization_long_discrete')
"""

local_kernel, global_kernel, this_weights, ramp = ramp_error_discrete(x1, xmin, xmax, induction, plot=False,
                                                                                   full_output=True,
                                                                                   use_global_saturation_factor=True)
peak_weight_index = np.where(this_weights['CA3'] == np.max(this_weights['CA3']))[0][0]
center_indexes = np.where((peak_locs['CA3'] >= peak_locs['CA3'][peak_weight_index] - 5.) &
                          (peak_locs['CA3'] < peak_locs['CA3'][peak_weight_index] + 5.))[0]
print 'peak_weights: %.2f' %  np.mean(this_weights['CA3'][center_indexes])


with h5py.File(data_dir+weights_filename+'.hdf5', 'a') as f:
    run_vel_key = str(int(induction_run_vel))
    if run_vel_key not in f:
        f.create_group(run_vel_key)
    f[run_vel_key].create_group(str(induction_seed))
    for group in peak_locs:
        f[run_vel_key][str(induction_seed)].create_group(group)
        f[run_vel_key][str(induction_seed)][group].create_dataset('weights', compression='gzip', compression_opts=9,
                                                                  data=this_weights[group]+1.)
        f[run_vel_key][str(induction_seed)][group].create_dataset('peak_locs', compression='gzip', compression_opts=9,
                                                                  data=peak_locs[group])
