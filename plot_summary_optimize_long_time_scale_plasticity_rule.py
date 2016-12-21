__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import random
import sys
import scipy.signal as signal
import matplotlib.gridspec as gridspec

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
# .hdf5 file contains traces for position vs. time for each induction trial, and ramp vs. position for each induction
if len(sys.argv) > 1:
    cell_id = str(sys.argv[1])
else:
    cell_id = None
if len(sys.argv) > 2:
    input_index = int(sys.argv[2])
else:
    input_index = None


# experimental_filename = '112116 magee lab first induction'
# experimental_filename = '112516 magee lab first induction'
experimental_filename = '121216 magee lab first induction'

rule_max_timescale = 9000.


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
    spatial_rate_maps = []
    phase_maps = []
    for group in ['CA3']:  # stim_exc_syns:
        for i, syn in enumerate(stim_exc_syns[group]):
            peak_loc = peak_locs[group][i]
            gauss_force = excitatory_peak_rate[group] * np.exp(-((extended_x - peak_loc) / gauss_sigma)**2.)
            gauss_force = wrap_around_and_compress(gauss_force, generic_x)
            spatial_rate_maps.append(gauss_force)
            if group in excitatory_precession_range:
                phase_force = get_dynamic_theta_phase_force(excitatory_precession_range[group], peak_loc, extended_x,
                                                            generic_x, generic_dx)
                phase_maps.append(phase_force)
            
    return spatial_rate_maps, phase_maps


def generate_complete_rate_maps(simiter, induction, spatial_rate_maps, phase_maps, global_phase_offset=None):
    """

    :param simiter: int
    :param induction: int
    :param spatial_rate_maps: array
    :param phase_maps: array
    :param global_phase_offset: float
    """
    local_random.seed(simiter)
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


def generate_default_rate_maps(simiter, spatial_rate_maps, phase_maps, global_phase_offset=None):
    """

    :param simiter: int
    :param spatial_rate_maps: array
    :param phase_maps: array
    :param global_phase_offset: float
    """
    local_random.seed(simiter)
    if global_phase_offset is None:
        global_phase_offset = local_random.uniform(-np.pi, np.pi)
    default_rate_maps = []
    group = 'CA3'
    for j in range(len(spatial_rate_maps)):
        this_rate_map = np.interp(default_interp_x, generic_x, spatial_rate_maps[j])
        this_phase_map = np.interp(default_interp_x, generic_x, phase_maps[j])
        theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] * np.cos(this_phase_map +
                                                                                  excitatory_theta_phase_offset[group] -
                                                                                       2. * np.pi * default_interp_t /
                                                                                  global_theta_cycle_duration +
                                                                                       global_phase_offset))
        theta_force -= np.min(theta_force)
        theta_force /= np.max(theta_force)
        theta_force *= excitatory_theta_modulation_depth[group]
        theta_force += 1. - excitatory_theta_modulation_depth[group]
        this_rate_map = np.multiply(this_rate_map, theta_force)
        default_rate_maps.append(this_rate_map)

    return default_rate_maps


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
    for i, rate_map in enumerate(rate_maps):
        this_weighted_rate_map = rate_map * weights[i] * scaling_factor
        this_weighted_rate_map = np.concatenate([this_weighted_rate_map for i in range(3)])
        this_weighted_rate_map = np.convolve(this_weighted_rate_map, epsp_filter)[:3*len(interp_x)]
        weighted_rate_maps.append(this_weighted_rate_map)
    expected_depolarization = np.sum(weighted_rate_maps, axis=0)
    expected_depolarization = low_pass_filter(expected_depolarization, 2.,
                                              3*len(interp_x)*dt, dt, 1.)[len(interp_x): 2 * len(interp_x)]
    expected_depolarization = np.interp(binned_x, interp_x, expected_depolarization)
    return expected_depolarization


NMDA_type = 'NMDA_KIN5'

dt = 1.  # ms
down_dt = 10.  # ms, to speed up optimization
equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # -0.5 * run_vel + 155.  # (ms)
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
    with h5py.File(data_dir+experimental_filename+'.hdf5', 'r') as f:
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

spatial_rate_maps, phase_maps = generate_generic_rate_and_phase_maps()

complete_t, complete_x, complete_rate_maps, complete_induction_gates = {}, {}, {}, {}
for induction in [1]:  # position:
    complete_t[induction], complete_x[induction], complete_rate_maps[induction] = \
        generate_complete_rate_maps(trial_seed, induction, spatial_rate_maps, phase_maps, default_global_phase_offset)
    complete_induction_gates[induction] = generate_complete_induction_gate(induction)

default_rate_maps = generate_default_rate_maps(trial_seed, spatial_rate_maps, phase_maps, default_global_phase_offset)

baseline = None
for induction in position:
    if baseline is None:
        ramp_baseline_indexes = np.where(np.array(ramp[induction]) <= np.percentile(ramp[induction], 10.))[0]
        baseline = np.mean(ramp[induction][ramp_baseline_indexes])
    ignore = subtract_baseline(ramp[induction], baseline)


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


def plot_plasticity_signal(x, local_kernel, global_kernel, induction, rule, saturation_factor, input_index=None):
    """
    Given the local and global kernels, convolve each input rate_map with the local kernel, and convolve the
    current injection with the global kernel. The weight change for each input is proportional to the area under the
    product of the two signals. Incremental weight changes accrue across multiple induction trials.
    :param x: array: [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio, kernel_scale]
    :param local_kernel: array
    :param global_kernel: array
    :param induction: int: key for dicts of arrays
    :param rule: str
    :param saturation_factor: float
    :param input_index: int
    """
    if input_index is None:
        t_index = len(interp_t[induction][0]) - int(1500. / dt)
        input_index = np.where(peak_locs['CA3'] >= interp_x[induction][0][t_index])[0][0]
    filter_ratio = x[4]
    kernel_scale = x[5]
    global_signal = np.convolve(complete_induction_gates[induction], global_kernel)[:len(complete_t[induction])] * \
                    kernel_scale
    this_stim_force = np.array(complete_rate_maps[induction][input_index])
    local_signal = np.convolve(0.001 * dt * this_stim_force, local_kernel)[:len(this_stim_force)] / \
                   saturation_factor * kernel_scale / filter_ratio
    this_signal = np.minimum(local_signal, global_signal)

    ylim = max(np.max(local_signal), np.max(global_signal))
    start_index = np.where(interp_x[induction][1] >= induction_locs[induction][1])[0][0]
    this_induction_start = interp_t[induction][1][start_index]
    this_induction_dur = induction_durs[induction][1]
    end_index = np.where(interp_t[induction][1] >= this_induction_start+this_induction_dur)[0][0]
    start_index += 2*len(interp_t[induction][0])
    end_index += 2 * len(interp_t[induction][0])
    this_duration = 6000.
    back_buffer = int(4000. / dt)
    forward_buffer = int(2000. / dt)
    x_start = 4000. / this_duration
    x_end = (4000. + this_induction_dur) / this_duration

    fig = plt.figure()
    gs = gridspec.GridSpec(3, 3)
    ax1 = plt.subplot(gs[0, :-1])
    ax1.plot(complete_t[induction][start_index - back_buffer:start_index + forward_buffer],
             local_signal[start_index - back_buffer:start_index + forward_buffer],
             label='Local plasticity signal (%s)' % rule, color='b')
    ax1.plot(complete_t[induction][start_index - back_buffer:start_index + forward_buffer],
             global_signal[start_index - back_buffer:start_index + forward_buffer],
             label='Global plasticity signal (%s)' % rule, color='k')
    ax1.fill_between(complete_t[induction][start_index - back_buffer:start_index + forward_buffer], 0.,
                     this_signal[start_index - back_buffer:start_index + forward_buffer], label='Overlap',
                     facecolor='b', alpha=0.5)
    ax1.axhline(y=ylim * 0.95, xmin=x_start, xmax=x_end, linewidth=3, c='k')
    ax1.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Signal amplitude (a.u.)')
    ax1.set_xlim(complete_t[induction][start_index - back_buffer], complete_t[induction][start_index + forward_buffer])
    ax1.set_ylim(0., ylim * 1.05)
    # axes.set_title('Induced plasticity signal')
    clean_axes(ax1)
    gs.tight_layout(fig)
    plt.show()
    plt.close()


x0 = {}

# x0['1'] = [1.972E+01, 1.216E+03, 4.765E+01, 2.681E+02, 7.022E-01, 2.635E-03]  # Error: 7.9850E+04
x0['1'] = [2.031E+01, 1.211E+03, 4.852E+01, 2.627E+02, 7.004E-01, 2.546E-03]  # Error: 7.9104E+04
# x0['2'] = [2.097E+01, 4.662E+02, 1.180E+01, 2.500E+01, 1.500E+00, 4.678E-03]  # Error: 1.2717E+05
x0['2'] = [2.868E+01, 4.608E+02, 1.377E+01, 2.500E+01, 1.500E+00, 5.052E-03]  # Error: 8.6251E+04
# Don't use cell3, it's the same as cell15
# x0['4'] = [8.155E+01, 1.624E+03, 3.712E+01, 1.967E+03, 1.499E+00, 8.198E-04]  # Error: 1.9124E+04
x0['4'] = [1.260E+02, 1.109E+03, 4.557E+01, 1.682E+03, 1.402E+00, 7.828E-04]  # Error: 1.3152E+04
# x0['5'] = [4.299E+02, 6.960E+02, 1.449E+01, 1.755E+02, 1.488E+00, 1.560E-03]  # Error: 1.4439E+04
x0['5'] = [4.820E+02, 6.828E+02, 1.000E+01, 1.565E+02, 1.418E+00, 1.474E-03]  #Error: 8.5719E+03
# x0['6'] = [2.017E+02, 3.000E+02, 1.231E+01, 2.615E+02, 1.445E+00, 2.682E-03]  # Error: 9.5538E+03
x0['6'] = [1.303E+02, 4.048E+02, 1.000E+01, 1.829E+02, 1.498E+00, 2.844E-03]  # Error: 1.0126E+04
# x0['7'] = [3.492E+02, 1.937E+03, 2.225E+01, 8.297E+02, 1.457E+00, 1.925E-03]  # Error: 9.6109E+04
x0['7'] = [3.006E+02, 2.072E+03, 2.326E+01, 7.880E+02, 1.500E+00, 1.959E-03]  # Error: 8.3039E+04
# x0['8'] = [3.342E+02, 4.438E+02, 3.464E+01, 2.306E+02, 1.232E+00, 1.795E-03]  # Error: 2.2846E+04
x0['8'] = [3.868E+02, 3.869E+02, 2.828E+01, 2.294E+02, 1.381E+00, 1.886E-03]  # Error: 1.7586E+04
# x0['9'] = [8.459E+01, 1.684E+03, 4.721E+01, 1.042E+03, 7.038E-01, 1.638E-03]  # Error: 4.8672E+04
x0['9'] = [8.912E+01, 1.651E+03, 1.137E+01, 1.220E+03, 7.057E-01, 1.663E-03]  # Error: 4.2240E+04
# x0['10'] = [5.000E+02, 4.975E+03, 1.393E+01, 2.500E+01, 1.381E+00, 4.110E-03]  # Error: 5.9070E+05
x0['10'] = [4.976E+02, 1.343E+03, 1.430E+01, 2.500E+01, 1.500E+00, 4.060E-03]  # Error: 7.0348E+04
# x0['11'] = [1.141E+01, 3.000E+02, 1.381E+01, 2.500E+01, 1.491E+00, 1.574E-03]  # Error: 2.0528E+04
x0['11'] = [1.146E+01, 3.000E+02, 1.372E+01, 2.509E+01, 1.500E+00, 1.637E-03]  # Error: 1.9207E+04
# x0['12'] = [6.790E+01, 5.297E+02, 2.998E+01, 4.299E+01, 1.469E+00, 4.953E-03]  # Error: 7.8383E+03
x0['12'] = [6.343E+01, 5.478E+02, 3.133E+01, 4.611E+01, 1.473E+00, 4.857E-03]  # Error: 6.6872E+03
# x0['13'] = [4.799E+02, 6.360E+02, 1.053E+01, 3.315E+01, 7.000E-01, 2.668E-03]  # Error: 2.5155E+03
x0['13'] = [4.975E+02, 6.337E+02, 1.077E+01, 3.313E+01, 7.004E-01, 2.635E-03]  # Error: 2.3783E+03
# x0['14'] = [1.302E+02, 3.000E+02, 1.544E+01, 1.630E+02, 1.253E+00, 1.256E-03]  # Error: 1.7495E+04
x0['14'] = [5.500E+01, 3.008E+02, 2.217E+01, 2.515E+02, 1.333E+00, 1.243E-03]  # Error: 5.9102E+03
# x0['15'] = [5.000E+02, 5.747E+02, 1.048E+01, 1.370E+02, 1.276E+00, 1.715E-03]  # Error: 3.6931E+04
x0['15'] = [4.348E+02, 6.156E+02, 1.000E+01, 3.663E+01, 1.225E+00, 2.099E-03]  # Error: 2.1445E+04
# Don't use cell16, it's the same as cell8
# x0['17'] = [5.000E+02, 1.805E+03, 2.395E+01, 4.724E+02, 1.371E+00, 1.762E-03]  # Error: 1.1026E+05
x0['17'] = [4.439E+02, 1.916E+03, 2.904E+01, 3.248E+02, 1.469E+00, 2.093E-03]  # Error: 1.0647E+05
# x0['18'] = [4.320E+02, 3.374E+03, 2.297E+01, 6.273E+02, 7.003E-01, 1.054E-03]  # Error: 8.3713E+04
x0['18'] = [4.920E+02, 4.974E+03, 2.950E+01, 1.475E+02, 7.015E-01, 1.200E-03]  # Error: 4.7170E+04
x0['19'] = [4.373E+01, 7.450E+02, 4.811E+01, 7.224E+02, 7.681E-01, 2.999E-03]  # Error: 2.1140E+04
x0['20'] = [6.382E+01, 2.153E+03, 4.881E+01, 1.092E+03, 1.129E+00, 1.524E-03]  # Error: 3.7379E+04
x0['21'] = [1.041E+01, 2.095E+03, 4.619E+01, 5.977E+02, 9.226E-01, 2.153E-03]  # Error: 1.7537E+04
x0['22'] = [6.949E+01, 3.162E+02, 4.587E+01, 4.365E+02, 8.259E-01, 2.938E-03]  # Error: 1.3404E+04
x0['23'] = [3.842E+02, 4.472E+02, 4.399E+01, 8.584E+02, 7.044E-01, 1.940E-03]  # Error: 3.2947E+04

# x0['mean'] = [2.795E+02, 1.236E+03, 2.268E+01, 4.522E+02, 1.397E+00, 2.763E-03]  # induced + spontaneous 120416
# x0['mean'] = [2.369E+02, 1.114E+03, 2.555E+01, 3.873E+02, 1.242E+00, 3.948E-03]  # induced + spontaneous 121116
x0['mean'] = [2.066E+02, 1.033E+03, 2.861E+01, 4.317E+02, 1.205E+00, 3.708E-03]  # induced + spontaneous 121316


# to avoid saturation and reduce variability of time courses across cells, constrain the relative amplitude
# of global and local kernels:
# [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio, kernel_scale]
if cell_id in x0:
    x1 = x0[cell_id]
else:
    x1 = x0['1']
xmin1 = [10., 300., 10., 25., 0.7, 5.e-4]
xmax1 = [500., 5000., 50., 2000., 1.5, 2.e-2]

global_saturation_factor = 2.056E-02  # cell 12

induction = 1

local_kernel, global_kernel = build_kernels(x1, False)

plot_plasticity_signal(x1, local_kernel, global_kernel, induction, 'long', global_saturation_factor)


"""
fig, axes = plt.subplots(2, 2)
fig.set_size_inches(5.2, 3.9)
mean_induction_loc = np.mean(induction_locs[induction])
mean_induction_dur = np.mean(induction_durs[induction])
start_index = np.where(interp_x[induction][0] >= mean_induction_loc)[0][0]
end_index = start_index + int(mean_induction_dur / dt)
x_start = mean_induction_loc/track_length
x_end = interp_x[induction][0][end_index] / track_length
ylim = max(np.max(ramp[induction]), np.max(model_ramp))
print 'ylim: ', ylim
ymin = min(np.min(ramp[induction]), np.min(model_ramp))
axes[0][0].plot(binned_x, ramp[induction], label='Experiment', color='k')
axes[0][0].plot(binned_x, model_ramp, label='Long model', color='b')
axes[0][0].axhline(y=ylim + 0.25, xmin=x_start, xmax=x_end, linewidth=2, c='k')
axes[0][0].set_ylabel('Depolarization (mV)')
axes[0][0].set_xlabel('Location (cm)')
axes[0][0].set_xlim(0., track_length)
axes[0][0].set_ylim(-0.5, ylim + 0.5)
axes[0][0].set_title('Induced Vm ramp', fontsize=12.)
axes[0][0].legend(loc='best', frameon=False, framealpha=0.5, fontsize=12.)
clean_axes(axes[0])
plt.tight_layout()
plt.show()
plt.close()
"""