__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import random
import sys
import scipy.signal as signal
import mkl
# import matplotlib.gridspec as gridspec

"""
These methods determine the shape of the function that translates plasticity signal into changes in synaptic weight
induced by plateaus.

Assumptions:
1) Synaptic weights are all = 1 prior to field induction 1. w(t0) = 1
2) The transfer function for field induction 1 is a linear function of plasticity signal, with some gain.

"""
morph_filename = 'EB2-late-bifurcation.swc'

mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'

synapses_seed = 0
num_exc_syns = 1600
num_inh_syns = 600
# allows parallel computation of multiple trials for the same spines with the same peak_locs, but with different
# input spike trains and stochastic synapses for each trial
trial_seed = 0
# .hdf5 file contains traces for position vs. time for each induction trial, and ramp vs. position for each induction
if len(sys.argv) > 1:
    cell_id = str(sys.argv[1])
else:
    cell_id = None

experimental_file_dir = data_dir
experimental_filename = '121216 magee lab first induction'

mkl.set_num_threads(4)


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


def calculate_ramp_features(ramp, induction_loc, offset=False):
    """

    :param ramp: array
    :param induction_loc: float
    :param offset: bool
    """
    extended_binned_x = np.concatenate([binned_x - track_length, binned_x, binned_x + track_length])
    smoothed_ramp = signal.savgol_filter(ramp, 19, 3, mode='wrap')
    extended_binned_ramp = np.concatenate([smoothed_ramp for i in range(3)])
    extended_interp_x = np.concatenate([default_interp_x - track_length, default_interp_x,
                                        default_interp_x + track_length])
    extended_ramp = np.interp(extended_interp_x, extended_binned_x, extended_binned_ramp)
    interp_ramp = extended_ramp[len(default_interp_x):2*len(default_interp_x)]
    baseline_indexes = np.where(interp_ramp <= np.percentile(interp_ramp, 10.))[0]
    baseline = np.mean(interp_ramp[baseline_indexes])
    if offset:
        interp_ramp -= baseline
        extended_ramp -= baseline
    peak_index = np.where(interp_ramp == np.max(interp_ramp))[0][0] + len(interp_ramp)
    peak_val = extended_ramp[peak_index]
    peak_x = extended_interp_x[peak_index]
    start_index = np.where(extended_ramp[:peak_index] <=
                           0.15*(peak_val - baseline) + baseline)[0][-1]
    end_index = peak_index + np.where(extended_ramp[peak_index:] <= 0.15*
                                                (peak_val - baseline) + baseline)[0][0]
    start_loc = float(start_index % len(default_interp_x)) / float(len(default_interp_x)) * track_length
    end_loc = float(end_index % len(default_interp_x)) / float(len(default_interp_x)) * track_length
    peak_shift = peak_x - induction_loc
    if peak_shift > track_length / 2.:
        peak_shift = -(track_length - peak_shift)
    elif peak_shift < -track_length / 2.:
        peak_shift += track_length
    ramp_width = extended_interp_x[end_index] - extended_interp_x[start_index]
    before_width = induction_loc - start_loc
    if induction_loc < start_loc:
        before_width += track_length
    after_width = end_loc - induction_loc
    if induction_loc > end_loc:
        after_width += track_length
    ratio = before_width / after_width
    return peak_val, ramp_width, peak_shift, ratio, start_loc, end_loc


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


def generate_theta_phase_maps():
    """
    Given a set of place field peak locations, return theta phase vs. location computed at a resolution of
    track_length/10000 bins.
    :return: array
    """
    phase_maps = []
    for group in ['CA3']:
        for i, peak_loc in enumerate(peak_locs[group]):
            if group in excitatory_precession_range:
                phase_force = get_dynamic_theta_phase_force(excitatory_precession_range[group], peak_loc, extended_x,
                                                            generic_x, generic_dx)
                phase_maps.append(phase_force)
    return phase_maps


def generate_spatial_rate_maps():
    """
    Given a set of place field peak locations, return firing rate vs. location computed at a resolution of
    track_length/10000 bins.
    :return: array
    """
    spatial_rate_maps = []
    for group in ['CA3']:
        for i, peak_loc in enumerate(peak_locs[group]):
            gauss_force = excitatory_peak_rate[group] * np.exp(-((extended_x - peak_loc) / gauss_sigma) ** 2.)
            gauss_force = wrap_around_and_compress(gauss_force, generic_x)
            spatial_rate_maps.append(gauss_force)
    return spatial_rate_maps


def generate_complete_rate_maps_theta(induction, spatial_rate_maps, phase_maps, global_phase_offset=None):
    """
    Use spatial maps for firing rate and theta phase, time vs. position maps for each induction trial, and a binary
    function of running velocity, to compute a set of complete spatial and temporal rate maps for the entire induction
    period. If no trajectory data was loaded for the laps before and after induction, this method duplicates the first
    and last lap.
    :param induction: int
    :param spatial_rate_maps: array
    :param phase_maps: array
    :param global_phase_offset: float
    """
    if global_phase_offset is None:
        global_phase_offset = local_random.uniform(-np.pi, np.pi)
    complete_rate_maps = []
    running_delta_t = 0.
    for i in range(len(position[induction])):
        this_interp_t = interp_t[induction][i]
        this_interp_x = interp_x[induction][i]
        this_run_vel_gate = run_vel_gate[induction][i]
        if i == 0:
            complete_t = this_interp_t - len(this_interp_t) * dt
            complete_x = this_interp_x - track_length
            complete_run_vel_gate = np.array(this_run_vel_gate)
        complete_t = np.append(complete_t, this_interp_t + running_delta_t)
        complete_x = np.append(complete_x, this_interp_x + i * track_length)
        complete_run_vel_gate = np.append(complete_run_vel_gate, this_run_vel_gate)
        running_delta_t += len(this_interp_t) * dt
        if i == len(position[induction]) - 1:
            complete_t = np.append(complete_t, this_interp_t + running_delta_t)
            complete_x = np.append(complete_x, this_interp_x + (i+1) * track_length)
            complete_run_vel_gate = np.append(complete_run_vel_gate, this_run_vel_gate)
    group = 'CA3'
    for j in range(len(spatial_rate_maps)):
        for i in range(len(position[induction])):
            this_interp_x = interp_x[induction][i]
            this_rate_map = np.interp(this_interp_x, generic_x, spatial_rate_maps[j])
            this_phase_map = np.interp(this_interp_x, generic_x, phase_maps[j])
            if i == 0:
                this_complete_rate_map = np.array(this_rate_map)
                this_complete_phase_map = np.array(this_phase_map)
            this_complete_rate_map = np.append(this_complete_rate_map, this_rate_map)
            this_complete_phase_map = np.append(this_complete_phase_map, this_phase_map)
            if i == len(position[induction]) - 1:
                this_complete_rate_map = np.append(this_complete_rate_map, this_rate_map)
                this_complete_phase_map = np.append(this_complete_phase_map, this_phase_map)
        theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] * np.cos(this_complete_phase_map +
                                                                                  excitatory_theta_phase_offset[group] -
                                                                                       2. * np.pi * complete_t /
                                                                                  global_theta_cycle_duration +
                                                                                       global_phase_offset))
        theta_force -= np.min(theta_force)
        theta_force /= np.max(theta_force)
        theta_force *= excitatory_theta_modulation_depth[group]
        theta_force += 1. - excitatory_theta_modulation_depth[group]
        this_complete_rate_map = np.multiply(this_complete_rate_map, theta_force)
        this_complete_rate_map = np.multiply(this_complete_rate_map, complete_run_vel_gate)
        complete_rate_maps.append(this_complete_rate_map)

    return complete_t, complete_x, complete_rate_maps


def generate_complete_rate_maps_no_theta(induction, spatial_rate_maps):
    """
    Use spatial maps for firing rate, time vs. position maps for each induction trial, and a binary function of
    running velocity, to compute a set of complete spatial and temporal rate maps for the entire induction period. If
    no trajectory data was loaded for the laps before and after induction, this method duplicates the first and last
    lap.
    :param induction: int
    :param spatial_rate_maps: array
    """
    complete_rate_maps = []
    running_delta_t = 0.
    for i in range(len(position[induction])):
        this_interp_t = interp_t[induction][i]
        this_interp_x = interp_x[induction][i]
        this_run_vel_gate = run_vel_gate[induction][i]
        if i == 0:
            complete_t = this_interp_t - len(this_interp_t) * dt
            complete_x = this_interp_x - track_length
            complete_run_vel_gate = np.array(this_run_vel_gate)
        complete_t = np.append(complete_t, this_interp_t + running_delta_t)
        complete_x = np.append(complete_x, this_interp_x + i * track_length)
        complete_run_vel_gate = np.append(complete_run_vel_gate, this_run_vel_gate)
        running_delta_t += len(this_interp_t) * dt
        if i == len(position[induction]) - 1:
            complete_t = np.append(complete_t, this_interp_t + running_delta_t)
            complete_x = np.append(complete_x, this_interp_x + (i+1) * track_length)
            complete_run_vel_gate = np.append(complete_run_vel_gate, this_run_vel_gate)
    group = 'CA3'
    for j in range(len(spatial_rate_maps)):
        for i in range(len(position[induction])):
            this_interp_x = interp_x[induction][i]
            this_rate_map = np.interp(this_interp_x, generic_x, spatial_rate_maps[j])
            if i == 0:
                this_complete_rate_map = np.array(this_rate_map)
            this_complete_rate_map = np.append(this_complete_rate_map, this_rate_map)
            if i == len(position[induction]) - 1:
                this_complete_rate_map = np.append(this_complete_rate_map, this_rate_map)
        this_complete_rate_map = np.multiply(this_complete_rate_map, complete_run_vel_gate)
        if len(this_complete_rate_map) != len(complete_run_vel_gate):
            print 'generate_complete_rate_maps_no_theta: Mismatched array length'
        complete_rate_maps.append(this_complete_rate_map)
    return complete_t, complete_x, complete_rate_maps


def generate_complete_induction_gate(induction):
    """

    :param induction: int
    :return:
    """
    for i in range(len(position[induction])):
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
        if i == len(position[induction]) - 1:
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


def compute_EPSP_matrix(rate_maps, this_interp_x):
    """

    :param rate_maps: list of array
    :param this_interp_x: array
    :return: array
    """
    filter_t = np.arange(0., 200., dt)
    epsp_filter = np.exp(-filter_t / 20.) - np.exp(-filter_t / .5)
    epsp_filter /= np.sum(epsp_filter)
    EPSP_maps = []
    scaling_factor = 4.150E-04  # generates a predicted 6 mV depolarization from gaussian weights with peak = 2.5
    for i, rate_map in enumerate(rate_maps):
        this_EPSP_map = np.interp(default_interp_x, this_interp_x, rate_map) * scaling_factor
        this_EPSP_map = np.concatenate([this_EPSP_map for i in range(3)])
        this_EPSP_map = np.convolve(this_EPSP_map, epsp_filter)[:3 * len(default_interp_x)]
        EPSP_maps.append(this_EPSP_map[len(default_interp_x):2 * len(default_interp_x)])
    return np.array(EPSP_maps)


def generate_spike_trains(rate_maps, t):
    """

    :param rate_maps: list of array
    :return: list of array
    """
    spike_trains = []
    for rate_map in rate_maps:
        this_spike_train = get_inhom_poisson_spike_times_by_thinning(rate_map, t, dt=0.02, generator=local_random)
        spike_trains.append(this_spike_train)
    return spike_trains


def filter_spike_trains(spike_trains):
    """

    :param spike_trains: list of array
    :return: list of array
    """
    dynamics = [0.2, 1.769, 67.351, 0.878, 92.918]
    successes = []
    for spike_train in spike_trains:
        this_Pr = Pr(*dynamics)
        this_success_train = []
        for spike_time in spike_train:
            P = this_Pr.stim(spike_time)
            if local_random.random() < P:
                this_success_train.append(spike_time)
        successes.append(this_success_train)
    return successes


NMDA_type = 'NMDA_KIN5'

dt = 1.  # ms
down_dt = 10.  # ms, to speed up optimization
equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # ms
input_field_width = 90.  # cm
track_length = 187.  # cm

binned_dx = track_length / 100.  # cm
binned_x = np.arange(0., track_length+binned_dx/2., binned_dx)[:100] + binned_dx/2.
generic_dx = binned_dx / 100.  # cm
generic_x = np.arange(0., track_length, generic_dx)

default_run_vel = 30.  # cm/s
generic_position_dt = generic_dx / default_run_vel * 1000.
generic_t = np.arange(0., len(generic_x)*generic_position_dt, generic_position_dt)[:len(generic_x)]
generic_track_duration = len(generic_t) * generic_position_dt

default_interp_t = np.arange(0., generic_t[-1], dt)
default_interp_x = np.interp(default_interp_t, generic_t, generic_x)

extended_x = np.concatenate([generic_x - track_length, generic_x, generic_x + track_length])

position = {}
ramp = {}
induction_locs = {}
induction_durs = {}
t = {}
interp_t = {}
interp_x = {}
run_vel = {}
run_vel_gate = {}


if cell_id is not None:
    with h5py.File(experimental_file_dir+experimental_filename+'.hdf5', 'r') as f:
        sampling_rate = f[cell_id].attrs['sampling_rate']
        track_length = f[cell_id].attrs['track_length']
        position_dt = 1000./sampling_rate  # convert s to ms
        for induction in f[cell_id]:
            position[int(induction)] = []
            induction_locs[int(induction)] = []
            induction_durs[int(induction)] = []
            t[int(induction)] = []
            for trial in f[cell_id][induction]['position'].itervalues():
                this_position = np.array(trial) / np.max(trial) * track_length
                position[int(induction)].append(this_position)
                t[int(induction)].append(
                    np.arange(0., len(this_position)*position_dt, position_dt)[:len(this_position)])
                induction_locs[int(induction)].append(trial.attrs['induction_loc'])
                induction_durs[int(induction)].append(trial.attrs['induction_dur'])
            ramp[int(induction)] = signal.savgol_filter(f[cell_id][induction]['ramp'][:], 19, 4, mode='wrap')
else:
    raise Exception('No data imported')

vel_window_dur = 10.  # ms
vel_window_bins = int(vel_window_dur/dt)/2

for induction in position:
    interp_t[induction] = []
    interp_x[induction] = []
    run_vel[induction] = []
    run_vel_gate[induction] = []
    for i in range(len(position[induction])):
        this_interp_t = np.arange(0., t[induction][i][-1] + dt / 2., dt)
        interp_t[induction].append(this_interp_t)
        this_interp_x = np.interp(interp_t[induction][i], t[induction][i], position[induction][i])
        interp_x[induction].append(this_interp_x)
        padded_t = np.insert(this_interp_t, 0, this_interp_t[-vel_window_bins:] - this_interp_t[-1] - dt)
        padded_t = np.append(padded_t, this_interp_t[:vel_window_bins + 1] + this_interp_t[-1] + dt)
        padded_x = np.insert(this_interp_x, 0, this_interp_x[-vel_window_bins:] - track_length)
        padded_x = np.append(padded_x, this_interp_x[:vel_window_bins + 1] + track_length)
        this_run_vel = []
        for j in range(vel_window_bins, vel_window_bins+len(this_interp_x)):
            this_run_vel.append(np.sum(np.diff(padded_x[j - vel_window_bins:j + vel_window_bins + 1])) /
                                np.sum(np.diff(padded_t[j - vel_window_bins:j + vel_window_bins + 1])) * 1000.)
        this_run_vel = np.array(this_run_vel)
        this_run_vel_gate = np.zeros_like(this_run_vel)
        this_run_vel_gate[np.where(this_run_vel>5.)[0]] = 1.
        run_vel[induction].append(this_run_vel)
        run_vel_gate[induction].append(this_run_vel_gate)


track_equilibrate = 2. * global_theta_cycle_duration
track_duration = {induction:
                      [interp_t[induction][i][-1] for i in range(len(position[induction]))] for induction in position}

excitatory_peak_rate = {'CA3': 40., 'ECIII': 40.}
excitatory_theta_modulation_depth = {'CA3': 0.7, 'ECIII': 0.7}
# From Chadwick et al., ELife 2015
excitatory_theta_phase_tuning_factor = {'CA3': 0.8, 'ECIII': 0.8}
excitatory_precession_range = {}  # # (ms, degrees)
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

spatial_rate_maps = generate_spatial_rate_maps()  # x=generic_x
theta_phase_maps = generate_theta_phase_maps()

local_random.seed(10000+int(cell_id))

complete_t, complete_x, complete_rate_maps, complete_induction_gates, spike_trains, successes = {}, {}, {}, {}, {}, {}
for induction in [1]:  # position:
    complete_t[induction], complete_x[induction], complete_rate_maps[induction] = \
        generate_complete_rate_maps_theta(induction, spatial_rate_maps, theta_phase_maps)
    complete_induction_gates[induction] = generate_complete_induction_gate(induction)

spike_trains[induction] = generate_spike_trains(complete_rate_maps[induction], complete_t[induction])
successes[induction] = filter_spike_trains(spike_trains[induction])


input_matrix = compute_EPSP_matrix(spatial_rate_maps, generic_x)  # x=default_interp_x

baseline = None
for induction in position:
    if baseline is None:
        ramp_baseline_indexes = np.where(np.array(ramp[induction]) <= np.percentile(ramp[induction], 10.))[0]
        baseline = np.mean(ramp[induction][ramp_baseline_indexes])
    ignore = subtract_baseline(ramp[induction], baseline)


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
    :param x: array: [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio]
    :param plot: bool
    :return: array, array
    """
    local_rise_tau = x[0]
    local_decay_tau = x[1]
    global_rise_tau = x[2]
    global_decay_tau = x[3]

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
        axes.plot(filter_t[:len(local_filter)]/1000., local_filter / np.max(local_filter), color='k',
                  label='Local plasticity kernel')
        axes.plot(filter_t[:len(global_filter)]/1000., global_filter / np.max(global_filter), color='b',
                  label='Global plasticity kernel')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_xlabel('Time (s)')
        axes.set_ylabel('Kernel ampl. (norm.)')
        axes.set_xlim(-0.5, max(5000., max(filter_t[:len(local_filter)][-1], filter_t[:len(global_filter)][-1]))/1000.)
        clean_axes(axes)
        fig.tight_layout()
        plt.show()
        plt.close()

    return local_filter, global_filter


def calculate_plasticity_signal(x, local_kernel, global_kernel, spike_trains, induction, plot=False):
    """
    Given the local and global kernels, convolve each input rate_map with the local kernel, and convolve the
    current injection with the global kernel. The weight change for each input is proportional to the area under the
    overlap of the two signals. Incremental weight changes accrue across multiple induction trials.
    :param x: array: [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio]
    :param local_kernel: array
    :param global_kernel: array
    :param spike_trains: dict of array
    :param induction: int: key for dicts of arrays
    :param plot: bool
    :return: plasticity_signal: array
    """
    filter_ratio = x[4]
    this_kernel_scale = kernel_scale['mean']
    group = 'CA3'
    local_signal_peaks = []
    local_signal_array = []
    plasticity_signal = np.zeros_like(peak_locs[group])
    global_signal = np.convolve(complete_induction_gates[induction], global_kernel)[:len(complete_t[induction])]
    max_global_signal = np.max(global_signal)
    for j, train in enumerate(spike_trains[induction]):
        indexes = (np.array(train-complete_t[induction][0]) / dt).astype(int)
        this_stim_force = np.zeros_like(complete_t[induction])
        this_stim_force[indexes] = 1.
        this_local_signal = np.convolve(this_stim_force, local_kernel)[:len(complete_t[induction])]
        local_signal_array.append(this_local_signal)
        local_signal_peaks.append(np.max(this_local_signal))
    max_local_signal = np.mean(np.array(local_signal_peaks)[np.where(local_signal_peaks >=
                                                           np.percentile(local_signal_peaks, 90.))[0]])
    saturation_factor = filter_ratio * max_local_signal / max_global_signal
    # print 'saturation factor: %.3E' % saturation_factor
    for j, train in enumerate(spike_trains[induction]):
        this_local_signal = local_signal_array[j] / saturation_factor
        this_signal = np.minimum(this_local_signal, global_signal)
        this_area = np.trapz(this_signal, dx=dt)
        plasticity_signal[j] += this_area
        if plot and j == int(len(spike_trains[induction])/2):
            # buffer = 5000.
            buffer = 0.
            # orig_font_size = mpl.rcParams['font.size']
            # orig_fig_size = mpl.rcParams['figure.figsize']
            # mpl.rcParams['font.size'] = 8.
            # mpl.rcParams['figure.figsize'] = 7.34, 3.25
            # fig1 = plt.figure()
            # gs1 = gridspec.GridSpec(2, 2)
            # axes = plt.subplot(gs1[0, 0])
            fig1, axes = plt.subplots(1)
            ylim = max(np.max(this_local_signal), max_global_signal)
            ylim *= this_kernel_scale
            this_global_signal = np.multiply(global_signal, this_kernel_scale)
            this_local_signal *= this_kernel_scale
            this_signal *= this_kernel_scale
            start_index = np.where(interp_x[induction][0] >= induction_locs[induction][0])[0][0]
            this_induction_start = interp_t[induction][0][start_index]
            this_induction_dur = induction_durs[induction][0]
            start_time = -buffer
            end_time = interp_t[induction][0][-1] + buffer
            this_duration = end_time - start_time
            x_start = (buffer + this_induction_start) / this_duration
            x_end = (buffer + this_induction_start + this_induction_dur) / this_duration
            axes.plot(complete_t[induction] / 1000., this_global_signal, label='Global signal', color='b')
            axes.plot(complete_t[induction]/1000., this_local_signal, label='Local signal', color='k')
            axes.fill_between(complete_t[induction]/1000., 0., this_signal, label='Overlap', facecolor='grey',
                              alpha=0.5)
            axes.axhline(y=ylim*1.05, xmin=x_start, xmax=x_end, linewidth=1, c='k')
            axes.legend(loc='best', frameon=False, framealpha=0.5)
            axes.set_xlabel('Time (s)')
            axes.set_ylabel('Signal amplitude (a.u.)')
            axes.set_xlim(-buffer/1000., interp_t[induction][0][-1]/1000. + buffer/1000.)
            axes.set_ylim(-0.05*ylim, ylim*1.1)
            axes.set_title('Plasticity signal')
            clean_axes(axes)
            # gs1.tight_layout(fig1)
            fig1.tight_layout()
            plt.show()
            plt.close()
            # mpl.rcParams['font.size'] = orig_font_size
            # mpl.rcParams['figure.figsize'] = orig_fig_size

    if plot:
        x_start = np.mean(induction_locs[induction]) / track_length
        ylim = np.max(plasticity_signal) * this_kernel_scale
        fig, axes = plt.subplots(1)
        axes.plot(peak_locs['CA3'], plasticity_signal * this_kernel_scale, color='r')
        axes.axhline(y=ylim + 0.1, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='k')
        axes.set_ylabel('Plasticity signal (a.u.)')
        axes.set_xlabel('Location (cm)')
        axes.set_xlim(0., track_length)
        axes.set_title('Induction %i: Plasticity signal' % induction)
        clean_axes(axes)
        fig.tight_layout()
        plt.show()
        plt.close()

    return plasticity_signal, saturation_factor


def ramp_error_parametric(x, xmin, xmax, input_matrix, spike_trains, ramp, induction=None, transform=None,
                          baseline=None, plot=False, full_output=False):
    """
    Given time courses of rise and decay for local and global plasticity kernels, and run velocities during field
    induction, this method calculates a set of synaptic weights to match the first place field induction.
    :param x: array [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio]
    :param xmin: array
    :param xmax: array
    :param input_matrix: array
    :param spike_trains: dict of list of array
    :param ramp: dict of array
    :param induction: int: key for dicts of arrays
    :param transform: callable function
    :param baseline: float
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    if induction is None:
        induction = 1
    this_induction_loc = np.mean([induction_loc for induction_loc in induction_locs[induction] if
                                  induction_loc is not None])
    print 'Trying x: %s for cell %s, induction_loc: %.1f' % (formatted_x, cell_id, this_induction_loc)
    if not check_bounds(x, xmin, xmax) or x[3] <= x[2]:
        print 'Aborting: Invalid parameter values.'
        hist.x.append(x)
        hist.Err.append(1e9)
        return 1e9
    exp_ramp = np.array(ramp[induction])
    start_time = time.time()
    local_kernel, global_kernel = build_kernels(x, plot)
    delta_weights, saturation_factor = calculate_plasticity_signal(x, local_kernel, global_kernel, spike_trains,
                                                                   induction, plot)
    amp, width, peak_shift, ratio, start_loc, end_loc = {}, {}, {}, {}, {}, {}
    amp['exp'], width['exp'], peak_shift['exp'], ratio['exp'], start_loc['exp'], end_loc['exp'] = \
        calculate_ramp_features(exp_ramp, this_induction_loc)
    plasticity_signal = np.multiply(delta_weights, kernel_scale['mean'])
    this_kernel_scale = 1.
    for attempt in range(2):
        weights = np.add(np.multiply(delta_weights, this_kernel_scale), 1.)
        model_ramp = weights.dot(input_matrix)  # x=default_interp_x
        if baseline is None:
            model_baseline = subtract_baseline(model_ramp)
        else:
            model_baseline = baseline
            model_ramp -= model_baseline
        if attempt == 0:
            this_kernel_scale = amp['exp'] / np.max(model_ramp)
    model_ramp = np.interp(binned_x, default_interp_x, model_ramp)
    amp['model'], width['model'], peak_shift['model'], ratio['model'], start_loc['model'], end_loc['model'] = \
        calculate_ramp_features(model_ramp, this_induction_loc)
    Err = 0.
    for feature, sigma in zip((width, peak_shift, ratio), (0.1, 0.05, 0.05)):
        Err += ((feature['exp'] - feature['model']) / sigma) ** 2.

    delta_start = abs(start_loc['exp'] - start_loc['model'])
    if delta_start > track_length / 2.:
        delta_start = track_length - delta_start
    Err += (delta_start / 0.1) ** 2.

    delta_end = abs(end_loc['exp'] - end_loc['model'])
    if delta_end > track_length / 2.:
        delta_end = track_length - delta_end
    Err += (delta_end / 0.1) ** 2.

    for j in range(len(exp_ramp)):
        Err += ((exp_ramp[j] - model_ramp[j]) / 0.1) ** 2.

    # penalize DC drifts in minimum weight
    Err += ((np.min(weights) - 1.)/0.005) ** 2.

    if plot:
        x_start = this_induction_loc/track_length
        ylim = max(np.max(exp_ramp), np.max(model_ramp))
        ymin = min(np.min(exp_ramp), np.min(model_ramp))
        fig, axes = plt.subplots(1)
        axes.plot(binned_x, exp_ramp, label='Exp. data', color='b')
        axes.plot(binned_x, model_ramp, label='Long signal integration model', color='k')
        axes.axhline(y=ylim + 0.2, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='k')
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Vm (mV)')
        axes.set_xlim([0., track_length])
        axes.set_ylim([math.floor(ymin), max(math.ceil(ylim), ylim + 0.4)])
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_title('Induction %i: Vm ramp' % induction)
        clean_axes(axes)
        fig.tight_layout()
        ylim = np.max(weights)
        ymin = np.min(weights)
        fig1, axes1 = plt.subplots(1)
        axes1.plot(peak_locs['CA3'], weights, c='r')
        axes1.axhline(y=ylim + 0.2, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='k')
        axes1.set_xlabel('Location (cm)')
        axes1.set_ylabel('Candidate synaptic weights')
        axes1.set_title('Induction %i: Synaptic weights' % induction)
        axes1.set_xlim([0., track_length])
        axes1.set_ylim([math.floor(ymin), max(math.ceil(ylim), ylim + 0.4)])
        clean_axes(axes1)
        fig1.tight_layout()
        plt.show()
        plt.close()

    print 'cell: %s, x: %s, kernel_scale: %.3E, Err: %.4E took %i s' % (cell_id, formatted_x, this_kernel_scale, Err,
                                                                        time.time()-start_time)
    print 'exp: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, end_loc: %.1f' % \
          (amp['exp'], width['exp'], peak_shift['exp'], ratio['exp'], start_loc['exp'], end_loc['exp'])
    print 'model: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, end_loc: %.1f' % \
          (amp['model'], width['model'], peak_shift['model'], ratio['model'], start_loc['model'], end_loc['model'])
    sys.stdout.flush()
    hist.x.append(x)
    hist.Err.append(Err)
    if full_output:
        return local_kernel, global_kernel, plasticity_signal, weights, model_ramp, model_baseline, this_kernel_scale, \
               saturation_factor, Err
    else:
        return Err


def estimate_weights_nonparametric(ramp, input_matrix, induction=None, baseline=None, plot=False, full_output=False):
    """
    Uses singular value decomposition to estimate a set of weights to match any arbitrary place field ramp, agnostic
    about underlying kernel, induction velocity, etc.
    :param ramp: dict of array
    :param input_matrix: array; x=default_interp_x
    :param induction: int: key for dicts of arrays
    :param baseline: float
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
    alpha = 2.  # regularization parameter
    if induction is None:
        induction = 1
    exp_ramp = np.array(ramp[induction])
    this_induction_loc = np.mean([induction_loc for induction_loc in induction_locs[induction] if
                                  induction_loc is not None])
    start_time = time.time()
    upsampled_ramp = np.interp(default_interp_x, binned_x, exp_ramp)
    """
    # synthetic ramp with amp = 6 mV for peak weight = 2.5
    modulated_field_center = track_length * 0.4
    peak_mod_weight = 2.5
    tuning_amp = (peak_mod_weight - 1.) / 2.
    tuning_offset = tuning_amp + 1.
    group = 'CA3'
    forced_weights = tuning_amp * np.cos(2. * np.pi / (input_field_width * 1.2) *
                                              (peak_locs[group] - modulated_field_center)) + tuning_offset
    left = np.where(peak_locs[group] >= modulated_field_center - input_field_width * 1.2 / 2.)[0][0]
    right = np.where(peak_locs[group] > modulated_field_center + input_field_width * 1.2 / 2.)[0][0]
    forced_weights[:left] = 1.
    forced_weights[right:] = 1.
    upsampled_ramp = forced_weights.dot(input_matrix)  # x=default_interp_x
    upsampled_ramp -= np.min(upsampled_ramp)
    exp_ramp = np.interp(binned_x, default_interp_x, upsampled_ramp)
    subtract_baseline(exp_ramp)
    """
    amp, width, peak_shift, ratio, start_loc, end_loc = {}, {}, {}, {}, {}, {}
    amp['exp'], width['exp'], peak_shift['exp'], ratio['exp'], start_loc['exp'], end_loc['exp'] = \
        calculate_ramp_features(exp_ramp, this_induction_loc)
    [U, s, Vh] = np.linalg.svd(input_matrix)
    V = Vh.T
    D = np.zeros_like(input_matrix)

    D[np.where(np.eye(*D.shape))] = s / (s ** 2. + alpha ** 2.)
    input_matrix_inv = V.dot(D.conj().T).dot(U.conj().T)

    delta_weights = upsampled_ramp.dot(input_matrix_inv)
    amp_factor = 1.
    for attempt in range(2):
        weights = np.add(np.multiply(delta_weights, amp_factor), 1.)
        model_ramp = weights.dot(input_matrix)  # x=default_interp_x
        if baseline is None:
            model_baseline = subtract_baseline(model_ramp)
        else:
            model_baseline = baseline
            model_ramp -= model_baseline
        if attempt == 0:
            amp_factor = amp['exp'] / np.max(model_ramp)
    model_ramp = np.interp(binned_x, default_interp_x, model_ramp)
    amp['model'], width['model'], peak_shift['model'], ratio['model'], start_loc['model'], end_loc['model'] = \
        calculate_ramp_features(model_ramp, this_induction_loc)
    asymmetry = np.zeros_like(peak_locs['CA3'])
    start_index = int(start_loc['exp']/track_length*len(peak_locs['CA3']))
    induction_index = int(this_induction_loc/track_length*len(peak_locs['CA3']))
    end_index = int(end_loc['exp']/track_length*len(peak_locs['CA3']))
    if start_loc['exp'] > this_induction_loc:
        asymmetry[start_index:] = 1.
        asymmetry[:induction_index] = 1.
        asymmetry[induction_index:end_index] = 2.
    else:
        asymmetry[start_index:induction_index] = 1.
        if end_loc['exp'] > this_induction_loc:
            asymmetry[induction_index:end_index] = 2.
        else:
            asymmetry[:end_index] = 2.
            asymmetry[induction_index:] = 2.
    Err = 0.
    Err += ((amp['exp'] - amp['model']) / 0.01) ** 2.

    for j in range(len(exp_ramp)):
        Err += ((exp_ramp[j] - model_ramp[j]) / 0.01) ** 2.

    # penalize DC drifts in minimum weight
    if induction == 1:
        Err += ((np.min(weights)-1.)/0.005) ** 2.
    elif induction == 2:
        Err += ((np.min(exp_ramp) - np.min(model_ramp)) / 0.005) ** 2.
    if plot:
        x_start = this_induction_loc/track_length
        ylim = max(np.max(exp_ramp), np.max(model_ramp))
        ymin = min(np.min(exp_ramp), np.min(model_ramp))
        fig, axes = plt.subplots(1)
        axes.plot(binned_x, exp_ramp, label='Experiment', color='k')
        axes.plot(binned_x, model_ramp, label='Model (SVD)', color='c')
        axes.axhline(y=ylim + 0.2, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='k')
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Depolarization (mV)')
        axes.set_xlim([0., track_length])
        axes.set_ylim([math.floor(ymin), max(math.ceil(ylim), ylim + 0.4)])
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_title('Induced Vm ramp')
        clean_axes(axes)
        fig.tight_layout()
        ylim = np.max(weights)
        ymin = np.min(weights)
        fig1, axes1 = plt.subplots(1)
        axes1.plot(peak_locs['CA3'], weights, c='r')
        axes1.axhline(y=ylim + 0.2, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='k')
        axes1.set_xlabel('Location (cm)')
        axes1.set_ylabel('Candidate synaptic weights')
        axes1.set_xlim([0., track_length])
        axes1.set_ylim([math.floor(ymin), max(math.ceil(ylim), ylim + 0.4)])
        clean_axes(axes1)
        fig1.tight_layout()
        plt.show()
        plt.close()

    print 'SVD took %i s, Err: %.4E' % (time.time()-start_time, Err)
    print 'exp: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, end_loc: %.1f' % \
          (amp['exp'], width['exp'], peak_shift['exp'], ratio['exp'], start_loc['exp'], end_loc['exp'])
    print 'model: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, end_loc: %.1f' % \
          (amp['model'], width['model'], peak_shift['model'], ratio['model'], start_loc['model'], end_loc['model'])
    sys.stdout.flush()
    if full_output:
        return weights, model_ramp, model_baseline, asymmetry
    else:
        return Err


def optimize_polish(x, xmin, xmax, error_function, input_matrix, spike_trains, ramp, induction=None,
                    transform=None, baseline=None, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param input_matrix: array
    :param spike_trains: dict of list of array
    :param ramp: dict of array
    :param induction: int: key for dicts of arrays
    :param transform: callable function
    :param baseline: float
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 600

    result = optimize.minimize(error_function, x, method='Nelder-Mead', options={'fatol': 1e-3, 'xatol': 1e-3,
                                                                                 'disp': True, 'maxiter': maxfev,
                                                                                 'maxfev': maxfev},
                               args=(xmin, xmax, input_matrix, spike_trains, ramp, induction, transform,
                                     baseline))
    formatted_x = '['+', '.join(['%.3E' % xi for xi in result.x])+']'
    print 'Process: %i completed optimize_polish for cell %s after %i iterations with Error: %.4E and x: %s' % \
          (os.getpid(), cell_id, result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


def optimize_explore(x, xmin, xmax, error_function, input_matrix, spike_trains, ramp, induction=None,
                     transform=None, baseline=None, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param input_matrix: array
    :param spike_trains: dict of list of array
    :param ramp: dict of array
    :param induction: int: key for dicts of arrays
    :param transform: callable function
    :param baseline: float
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 700

    take_step = Normalized_Step(x, xmin, xmax)
    minimizer_kwargs = dict(method=null_minimizer, args=(xmin, xmax, input_matrix, spike_trains, ramp,
                                                         induction, transform, baseline))
    result = optimize.basinhopping(error_function, x, niter=maxfev, niter_success=maxfev/2,
                                   disp=True, interval=min(20, int(maxfev/20)), minimizer_kwargs=minimizer_kwargs,
                                   take_step=take_step)
    formatted_x = '['+', '.join(['%.3E' % xi for xi in result.x])+']'
    print 'Process: %i completed optimize_explore for cell %s after %i iterations with Error: %.4E and x: %s' % \
          (os.getpid(), cell_id, result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


local_kernel = {}
global_kernel = {}
plasticity_signal = {}
weights_parametric = {}
model_ramp_parametric = {}
kernel_scale = {}
weights_SVD = {}
model_ramp_SVD = {}
model_baseline = {}
asymmetry = {}

x0 = {}

# to avoid saturation, ensure that the peak amplitude of the local signal is lower than the global signal:
# [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio]

x0['1'] = [1.045E+01, 1.627E+03, 1.357E+02, 2.508E+02, 7.501E-01]  # Error: 9.410E+04
x0['2'] = [3.338E+01, 5.001E+02, 1.299E+02, 1.348E+02, 1.500E+00]  # Error: 1.198E+05
x0['4'] = [1.124E+01, 3.170E+03, 2.943E+02, 9.815E+02, 1.477E+00]  # Error: 4.901E+04
x0['5'] = [4.946E+02, 6.667E+02, 1.003E+01, 5.271E+02, 7.514E-01]  # Error: 2.002E+03
x0['6'] = [1.415E+02, 5.148E+02, 1.448E+01, 2.891E+02, 1.435E+00]  # Error: 2.164E+04
x0['7'] = [3.627E+02, 1.732E+03, 2.149E+01, 1.256E+03, 7.623E-01]  # Error: 2.024E+05
# Jeff and Sandro excluded cell8
x0['8'] = [3.330E+02, 5.510E+02, 7.597E+01, 4.142E+02, 8.297E-01]  # Error: 4.575E+04
x0['9'] = [5.000E+02, 1.464E+03, 1.558E+01, 1.928E+03, 8.979E-01]  # Error: 1.343E+05
x0['10'] = [3.603E+02, 1.687E+03, 2.102E+01, 1.009E+02, 1.478E+00]  # Error: 1.382E+05
x0['11'] = [1.000E+01, 5.000E+02, 1.402E+01, 1.009E+02, 1.470E+00]  # Error: 4.512E+04
x0['12'] = [2.994E+01, 5.954E+02, 1.040E+02, 1.217E+02, 1.102E+00]  # Error: 7.386E+03
x0['13'] = [4.924E+02, 6.762E+02, 1.013E+01, 1.099E+02, 7.500E-01]  # Error: 6.205E+03
x0['14'] = [9.616E+01, 5.156E+02, 2.034E+01, 3.181E+02, 1.462E+00]  # Error: 1.894E+04
x0['15'] = [4.775E+02, 5.000E+02, 1.038E+01, 3.549E+02, 8.845E-01]  # Error: 2.872E+04
x0['17'] = [4.899E+02, 9.503E+02, 1.431E+01, 7.245E+02, 7.512E-01]  # Error: 2.188E+05
x0['18'] = [4.983E+02, 3.651E+03, 1.022E+01, 7.529E+02, 7.545E-01]  # Error: 8.871E+04
x0['19'] = [1.002E+01, 1.484E+03, 1.720E+02, 3.228E+02, 7.586E-01]  # Error: 1.568E+04
x0['20'] = [1.011E+01, 3.149E+03, 2.964E+02, 3.669E+02, 1.381E+00]  # Error: 7.976E+04
x0['21'] = [1.080E+01, 1.735E+03, 1.444E+01, 8.001E+02, 7.528E-01]  # Error: 4.522E+04
x0['22'] = [2.660E+01, 5.004E+02, 2.135E+02, 2.259E+02, 1.435E+00]  # Error: 1.230E+04
x0['23'] = [4.932E+02, 7.026E+02, 2.059E+02, 6.568E+02, 1.496E+00]  # Error: 2.252E+04

x0['mean'] = [2.129E+02, 1.196E+03, 1.015E+02, 4.548E+02, 1.116E+00]  # Induced + Spontaneous 032117

kernel_scale['1'] = 4.525E-03
kernel_scale['2'] = 2.591E-03
kernel_scale['4'] = 8.166E-04
kernel_scale['5'] = 1.539E-03
kernel_scale['6'] = 2.779E-03
kernel_scale['7'] = 1.813E-03
kernel_scale['8'] = 1.980E-03
kernel_scale['9'] = 2.110E-03
kernel_scale['10'] = 3.692E-03
kernel_scale['11'] = 1.310E-03
kernel_scale['12'] = 4.227E-03
kernel_scale['13'] = 3.022E-03
kernel_scale['14'] = 1.407E-03
kernel_scale['15'] = 1.842E-03
kernel_scale['17'] = 2.071E-03
kernel_scale['18'] = 1.219E-03
kernel_scale['19'] = 3.809E-03
kernel_scale['20'] = 1.628E-03
kernel_scale['21'] = 2.447E-03
kernel_scale['22'] = 3.898E-03
kernel_scale['23'] = 2.510E-03

kernel_scale['mean'] = 2.535E-03


if cell_id in x0:
    x1 = x0[cell_id]
else:
    x1 = x0['1']

# x1 = x0['mean']

xmin1 = [10., 500., 10., 100., 0.75]
xmax1 = [500., 5000., 300., 2000., 1.5]

for i in range(len(x1)):
    if x1[i] < xmin1[i]:
        x1[i] = xmin1[i]
    elif x1[i] > xmax1[i]:
        x1[i] = xmax1[i]

"""
induction = 1
polished_result = optimize_polish(x1, xmin1, xmax1, ramp_error_parametric, input_matrix, successes, ramp,
                                  induction, maxfev=600)
x1 = polished_result['x']
"""

for induction in position:
    if induction == 2 and 1 in position:
        this_model_baseline = model_baseline[1]
    else:
        this_model_baseline = None
    local_kernel[induction], global_kernel[induction], plasticity_signal[induction], weights_parametric[induction], \
        model_ramp_parametric[induction], model_baseline[induction], this_kernel_scale, this_saturation_factor, \
        Err = ramp_error_parametric(x1, xmin1, xmax1, input_matrix, successes, ramp, induction,
                                    baseline=this_model_baseline, plot=False, full_output=True)

if 1 not in plasticity_signal:
    plasticity_signal[1] = np.zeros_like(peak_locs['CA3'])
    weights_parametric[1] = np.ones_like(peak_locs['CA3'])

mean_induction_loc, mean_induction_dur = {}, {}

for induction in ramp:
    if induction == 2:
        this_model_baseline = model_baseline[1]
    else:
        this_model_baseline = None
    weights_SVD[induction], model_ramp_SVD[induction], model_baseline[induction], asymmetry[induction] = \
        estimate_weights_nonparametric(ramp, input_matrix, induction, plot=False, full_output=True)
    mean_induction_loc[induction] = np.mean([induction_loc for induction_loc in induction_locs[induction] if
                                             induction_loc is not None])
    mean_induction_dur[induction] = np.mean([induction_dur for induction_dur in induction_durs[induction] if
                                             induction_dur is not None])
"""
fig1, axes1 = plt.subplots(1)
fig2, axes2 = plt.subplots(1)

from matplotlib import cm
this_cm = cm.get_cmap()
colors = [this_cm(1.*i/2) for i in range(3)]

ylim1 = max(np.max(ramp.values()), np.max(model_ramp_SVD.values()), np.max(model_ramp_parametric.values()))
ylim2 = max(np.max(weights_SVD.values()), np.max(weights_parametric.values()))

for induction in ramp:
    start_index = np.where(interp_x[induction][1] >= mean_induction_loc[induction])[0][0]
    end_index = start_index + int(mean_induction_dur[induction] / dt)
    x_start = mean_induction_loc[induction] / track_length
    x_end = interp_x[induction][1][end_index] / track_length
    axes1.plot(binned_x, ramp[induction], label='Exp. data', c=colors[0])
    axes1.plot(binned_x, model_ramp_parametric[induction], label='Model (Long signal integration)', c=colors[1])
    axes1.plot(binned_x, model_ramp_SVD[induction], label='Model (SVD)', c=colors[2])
    axes1.axhline(y=ylim1 + 0.25, xmin=x_start, xmax=x_end, linewidth=2, c='k')
    axes2.plot(peak_locs['CA3'], weights_parametric[induction], label='Model (Long signal integration)', c=colors[1])
    axes2.plot(peak_locs['CA3'], weights_SVD[induction], label='Model (SVD)', c=colors[2])
    axes2.axhline(y=ylim2 + 0.25, xmin=x_start, xmax=x_end, linewidth=2, c='k')
axes1.set_xlabel('Location (cm)')
axes1.set_ylabel('Ramp depolarization (mV)')
axes1.set_xlim([0., track_length])
axes1.legend(loc='best', frameon=False, framealpha=0.5)
axes1.set_title('Induction 1: Vm ramp')
clean_axes(axes1)
fig1.tight_layout()
axes2.set_xlabel('Location (cm)')
axes2.set_ylabel('Candidate synaptic weights (a.u.)')
axes2.set_title('Induction 1: Synaptic weights')
axes2.legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes2)
fig2.tight_layout()

label_handles = []
label_handles.append(mlines.Line2D([], [], color=colors[0], label='Out of field'))
label_handles.append(mlines.Line2D([], [], color=colors[1], label='Before induction loc'))
label_handles.append(mlines.Line2D([], [], color=colors[2], label='After induction loc'))

delta_weights = {1: weights_SVD[1] - 1.}
if 2 in weights_SVD:
    delta_weights[2] = np.subtract(weights_SVD[2], weights_SVD[1])

for induction in ramp:
    if induction == 1:
        fig3, axes3 = plt.subplots(1)
        axes3.scatter(plasticity_signal[induction], delta_weights[induction], c=asymmetry[induction], linewidth=0)
        axes3.set_xlabel('Plasticity signal (a.u.)')
        axes3.set_ylabel('Change in synaptic weight (a.u.)')
        axes3.set_title('Induction %i: Plasticity transfer function' % induction)
        axes3.legend(loc='best', handles=label_handles, framealpha=0.5, frameon=False)
        clean_axes(axes3)
        fig3.tight_layout()
    elif induction == 2:
        fig4 = plt.figure()
        axes4 = fig4.add_subplot(111, projection='3d')
        axes4.scatter(plasticity_signal[induction], weights_SVD[1], delta_weights[induction], c=asymmetry[induction],
                      linewidth=0)
        axes4.set_xlabel('Plasticity signal (a.u.)')
        axes4.set_ylabel('Initial synaptic weight (a.u.)')
        axes4.set_zlabel('Change in synaptic weight (a.u.)')
        axes4.set_title('Induction %i: Plasticity transfer function' % induction)
        axes4.legend(loc='center left', handles=label_handles, framealpha=0.5, frameon=False, bbox_to_anchor=(1, 0.5))
        fig4.tight_layout()

plt.show()
plt.close()
"""


induction = 1
output_filename = '032617 discrete plasticity summary'
with h5py.File(data_dir+output_filename+'.hdf5', 'a') as f:
    if 'long' not in f:
        f.create_group('long')
    if 'position' not in f:
        f.create_dataset('position', compression='gzip', compression_opts=9, data=binned_x)
    if 'peak_locs' not in f:
        f.create_dataset('peak_locs', compression='gzip', compression_opts=9, data=peak_locs['CA3'])
    f['long'].create_group(cell_id)
    f['long'][cell_id].attrs['track_length'] = track_length
    f['long'][cell_id].attrs['induction_loc'] = mean_induction_loc[induction]
    f['long'][cell_id].attrs['induction_dur'] = mean_induction_dur[induction]
    f['long'][cell_id].attrs['parameters'] = x1
    f['long'][cell_id].attrs['kernel_scale'] = this_kernel_scale
    f['long'][cell_id].attrs['saturation_factor'] = this_saturation_factor
    f['long'][cell_id].attrs['error'] = Err
    f['long'][cell_id].create_dataset('local_kernel', compression='gzip', compression_opts=9,
                                      data=local_kernel[induction])
    f['long'][cell_id].create_dataset('global_kernel', compression='gzip', compression_opts=9,
                                      data=global_kernel[induction])
    f['long'][cell_id].attrs['dt'] = dt
    f['long'][cell_id].create_dataset('exp_ramp', compression='gzip', compression_opts=9, data=ramp[induction])
    f['long'][cell_id].create_dataset('model_ramp_parametric', compression='gzip', compression_opts=9,
                                      data=model_ramp_parametric[induction])
    f['long'][cell_id].create_dataset('model_ramp_SVD', compression='gzip', compression_opts=9,
                                      data=model_ramp_SVD[induction])
    f['long'][cell_id].create_dataset('plasticity_signal', compression='gzip', compression_opts=9,
                                      data=plasticity_signal[induction])
    f['long'][cell_id].create_dataset('weights_parametric', compression='gzip', compression_opts=9,
                                      data=weights_parametric[induction])
    f['long'][cell_id].create_dataset('weights_SVD', compression='gzip', compression_opts=9,
                                      data=weights_SVD[induction])
print 'Exported data for cell %s' % cell_id
