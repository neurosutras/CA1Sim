__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import random
import sys
import scipy.signal as signal

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
# .hdf5 file contains traces for position vs. time for each induction trial, and ramp vs. position for each induction
if len(sys.argv) > 1:
    cell_id = str(sys.argv[1])
else:
    cell_id = None
if len(sys.argv) > 2:
    rule_type = int(sys.argv[2])
else:
    rule_type = 1


# experimental_filename = '112116 magee lab first induction'
experimental_filename = '112516 magee lab first induction'

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


def wrap_around(waveform, interp_x):
    """

    :param waveform: array
    :param interp_x: array
    :return: array
    """
    before = np.array(waveform[:len(interp_x)])
    after = np.array(waveform[2 * len(interp_x):])
    within = np.array(waveform[len(interp_x):2 * len(interp_x)])
    waveform[len(interp_x):2 * len(interp_x)] += before + after
    waveform[2 * len(interp_x):] += before + within
    waveform[:len(interp_x)] += within + after
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
    phase_force[np.where(extended_x >= end_x_val)[0]] *= end_phase_val
    for i in range(len(phase_ranges) - 1):
        x0 = phase_ranges[i][0] + peak_loc
        x1 = phase_ranges[i + 1][0] + peak_loc
        phase0 = phase_ranges[i][1] * 2. * np.pi / 360.  # convert degrees to radians
        phase1 = phase_ranges[i + 1][1] * 2. * np.pi / 360.  # convert degrees to radians
        del_x = x1 - x0
        del_phase = phase1 - phase0
        if abs(del_x) > 0.:
            this_x_piece = np.arange(x0, x1, dx)
            this_indexes = np.where((extended_x >= x0) & (extended_x < x1))[0]
            this_interp_x_piece = extended_x[this_indexes]
            if abs(del_phase) > 0.:
                d_phase = del_phase / del_x * dx
                this_range_piece = np.arange(phase0, phase1, d_phase)
                this_range_interp_piece = np.interp(this_interp_x_piece, this_x_piece, this_range_piece)
            else:
                this_range_interp_piece = np.ones(len(this_interp_x_piece)) * phase0
            phase_force[this_indexes] = this_range_interp_piece
    phase_force = wrap_around(phase_force, interp_x)
    return phase_force


def generate_rate_maps(simiter, extended_x, interp_x, extended_t, run_vel_gate=None, global_phase_offset=None):
    """

    :param simiter: int
    :param extended_x: array (just one item from the outer list)
    :param interp_x: array (just one item from the outer list)
    :param extended_t: array (just one item from the outer list)
    :param run_vel_gate: array (just one item from the outer list)
    :param global_phase_offset: float
    """
    local_random.seed(simiter)
    if run_vel_gate is None:
        run_vel_gate = np.ones_like(interp_x)
    if global_phase_offset is None:
        global_phase_offset = local_random.uniform(-np.pi, np.pi)
    extended_run_vel_gate = np.append(np.append(run_vel_gate, run_vel_gate), run_vel_gate)
    rate_maps = []
    for group in ['CA3']:  # stim_exc_syns:
        for i, syn in enumerate(stim_exc_syns[group]):
            peak_loc = peak_locs[group][i]
            gauss_force = excitatory_peak_rate[group] * np.exp(-((extended_x - peak_loc) / gauss_sigma)**2.)
            gauss_force = wrap_around(gauss_force, interp_x)
            gauss_force = np.multiply(gauss_force, extended_run_vel_gate)
            if group in excitatory_precession_range:
                phase_force = get_dynamic_theta_phase_force(excitatory_precession_range[group], peak_loc, extended_x,
                                                            interp_x, dx/100.)
                theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] * np.cos(phase_force +
                                        excitatory_theta_phase_offset[group] - 2. * np.pi * extended_t /
                                        global_theta_cycle_duration + global_phase_offset))
            else:
                theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] *
                                 np.cos(excitatory_theta_phase_offset[group] - 2. * np.pi * extended_t /
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


def get_expected_depolarization(rate_maps, weights, extended_t, interp_x):
    """
    Take pre-computed rate maps and weights for a set of inputs. Convolve with an EPSP kernel, sum, normalize, and
    downsample to 100 spatial bins.
    :param rate_maps: list of array
    :param weights: list of float
    :param extended_t: array (just one item from the outer list)
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
        this_weighted_rate_map = np.convolve(this_weighted_rate_map, epsp_filter)[:len(extended_t)]
        weighted_rate_maps.append(this_weighted_rate_map)
    expected_depolarization = np.sum(weighted_rate_maps, axis=0)
    expected_depolarization = low_pass_filter(expected_depolarization, 2., len(extended_t)*dt, dt, 1.)
    expected_depolarization = np.interp(x_bins, interp_x, expected_depolarization[len(interp_x): 2 * len(interp_x)])
    return expected_depolarization


NMDA_type = 'NMDA_KIN5'

dt = 1.  # ms
down_dt = 10.  # ms, to speed up optimization
equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # -0.5 * run_vel + 155.  # (ms)
input_field_width = 90.  # cm
track_length = 187.  # cm
dx = track_length / 100.  # cm
x = np.arange(0., track_length+dx/2., dx)
x_bins = x[:100] + dx/2.
default_run_vel = 30.  # cm/s
default_position_dt = dx / default_run_vel * 1000.
default_t = np.arange(0., len(x)*default_position_dt, default_position_dt)[:len(x)]
default_interp_t = np.arange(0., default_t[-1] + dt / 2., dt)
default_interp_x = np.interp(default_interp_t, default_t, x)

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
    ramp[1] = None
    position[1] = []
    induction_locs[1] = []
    induction_durs[1] = []
    t[1] = []
    position_dt = dx / default_run_vel * 1000.
    for i in range(5):
        position[1].append(np.array(x))
        t[1].append(np.arange(0., len(x)*position_dt, position_dt)[:len(x)])
        induction_locs[1].append(track_length * field1_loc + 5.)  # optimize for backward shift of peak of ramp
        induction_durs[1].append(300.)

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
        padded_x = np.insert(this_interp_x, 0, this_interp_x[-vel_window_bins-1:-1]-track_length)
        padded_x = np.insert(padded_x, -1, this_interp_x[1:vel_window_bins+1]+track_length)
        this_run_vel = []
        for j in range(vel_window_bins, vel_window_bins+len(this_interp_x)):
            this_run_vel.append(np.sum(np.diff(padded_x[j-vel_window_bins:j+vel_window_bins]))/vel_window_dur*1000.)
        this_run_vel = np.array(this_run_vel)
        this_run_vel_gate = np.zeros_like(this_run_vel)
        this_run_vel_gate[np.where(this_run_vel>5.)[0]] = 1.
        run_vel[induction].append(this_run_vel)
        run_vel_gate[induction].append(this_run_vel_gate)


track_equilibrate = 2. * global_theta_cycle_duration
track_duration = {induction:
                      [interp_t[induction][i][-1] for i in range(len(position[induction]))] for induction in position}
default_track_duration = default_interp_t[-1]
duration = {induction:
                [equilibrate + track_equilibrate + track_duration[induction][i]
                 for i in range(len(position[induction]))] for induction in position}
default_duration = equilibrate + track_equilibrate + default_track_duration

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
sim.parameters['duration'] = default_duration
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

extended_x = {}
stim_t = {}
extended_t = {}
for induction in position:
    extended_x[induction] = []
    stim_t[induction] = []
    extended_t[induction] = []
    for i in range(len(position[induction])):
        this_interp_x = interp_x[induction][i]
        this_extended_x = np.zeros(3*len(this_interp_x))
        this_extended_x[:len(this_interp_x)] = this_interp_x - track_length
        this_extended_x[len(this_interp_x):2*len(this_interp_x)] = this_interp_x
        this_extended_x[2*len(this_interp_x):] = this_interp_x + track_length
        extended_x[induction].append(this_extended_x)
        stim_t[induction].append(np.insert(interp_t[induction][i], 0, np.arange(-track_equilibrate, 0., dt)))
        this_extended_t = np.zeros(3*len(interp_t[induction][i]))
        this_extended_t[:len(interp_t[induction][i])] = interp_t[induction][i] - track_duration[induction][i]
        this_extended_t[len(interp_t[induction][i]):2*len(interp_t[induction][i])] = interp_t[induction][i]
        this_extended_t[2*len(interp_t[induction][i]):] = interp_t[induction][i] + track_duration[induction][i]
        extended_t[induction].append(this_extended_t)

default_extended_x = np.zeros(3*len(default_interp_x))
default_extended_x[:len(default_interp_x)] = default_interp_x - track_length
default_extended_x[len(default_interp_x):2*len(default_interp_x)] = default_interp_x
default_extended_x[2*len(default_interp_x):] = default_interp_x + track_length
default_stim_t = np.insert(default_interp_t, 0, np.arange(-track_equilibrate, 0., dt))
default_extended_t = np.zeros(3*len(default_interp_t))
default_extended_t[:len(default_interp_t)] = default_interp_t - default_track_duration
default_extended_t[len(default_interp_t):2*len(default_interp_t)] = default_interp_t
default_extended_t[2*len(default_interp_t):] = default_interp_t + default_track_duration

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

default_global_phase_offset = 0.

rate_maps = {}
j = 0
for induction in [1]:  # position:
    rate_maps[induction] = []
    for i in range(len(position[induction])):
        rate_maps[induction].append(generate_rate_maps(trial_seed+j, extended_x[induction][i],
                                                       interp_x[induction][i], extended_t[induction][i],
                                                       run_vel_gate[induction][i]))
        j += 1
default_rate_maps = generate_rate_maps(trial_seed, default_extended_x, default_interp_x, default_extended_t,
                                       global_phase_offset=default_global_phase_offset)


# modulate the weights of inputs with peak_locs along this stretch of the track
field_center1 = track_length * field1_loc

ramp_baseline_indexes = {}

baseline = None
for induction in position:
    if ramp[induction] is None:
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
        ramp[induction] = get_expected_depolarization(rate_maps[induction][i], cos_mod_weight['CA3'],
                                                    interp_x[induction][i], stim_t[induction][i],
                                                    track_duration[induction][i])
    if baseline is None:
        ramp_baseline_indexes = np.where(np.array(ramp[induction]) <= np.percentile(ramp[induction], 10.))[0]
        baseline = np.mean(ramp[induction][ramp_baseline_indexes])
    ignore = subtract_baseline(ramp[induction], baseline)

"""
for induction in ramp:
    plt.plot(x_bins, ramp[induction])
plt.show()
plt.close()
"""


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


def build_kernel1(x, plot=False):
    """
    Construct two kernels with exponential rise and decay:
    1) Local kernel that generates a plasticity signal at each spine
    2) Global kernal that generates a plasticity signal during dendritic calcium spiking
    :param x: array: [local_rise_tau, local_decay_tau, local_scale, global_rise_tau, global_decay_tau, global_scale]
    :param plot: bool
    :return: array, array
    """
    local_rise_tau = x[0]
    local_decay_tau = x[1]
    global_rise_tau = x[2]
    global_decay_tau = x[3]
    kernel_scale = x[4]

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


def calculate_plasticity_signal(x, local_kernel, global_kernel, induction, plot=False):
    """
    Given the local and global kernels, convolve each input rate_map with the local kernel, and convolve the
    current injection with the global kernel. The weight change for each input is proportional to the area under the
    product of the two signals. Incremental weight changes accrue across multiple induction trials.
    :param x: array: [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, kernel_scale]
    :param local_kernel: array
    :param global_kernel: array
    :param induction: int: key for dicts of arrays
    :param plot: bool
    :return: plasticity_signal: array
    """
    saturation_factor = 0.02
    kernel_scale = x[4]
    group = 'CA3'
    for attempt in range(2):
        plasticity_signal = np.zeros_like(peak_locs[group])
        max_local_signal = 0.
        max_global_signal = 0.
        for i in range(len(induction_locs[induction])):
            start_time = time.time()
            this_rate_maps = rate_maps[induction][i]
            this_induction_loc = induction_locs[induction][i]
            this_induction_dur = induction_durs[induction][i]
            this_extended_x = extended_x[induction][i]
            this_extended_t = extended_t[induction][i]
            this_interp_t = interp_t[induction][i]
            global_signal = np.zeros_like(this_extended_x)
            start_index = np.where(this_extended_x >= this_induction_loc)[0][0]
            end_index = start_index + int(this_induction_dur / dt)
            global_signal[start_index:end_index] = 1.
            global_signal = np.convolve(global_signal, global_kernel)[:len(this_extended_x)] * kernel_scale
            down_t = np.arange(this_extended_t[0], this_extended_t[-1] + down_dt / 2., down_dt)
            global_signal = np.interp(down_t, this_extended_t, global_signal)
            max_global_signal = max(max_global_signal, np.max(global_signal))
            filter_t = np.arange(0., len(local_kernel) * dt, dt)
            down_filter_t = np.arange(0., filter_t[-1] + down_dt / 2., down_dt)
            local_kernel_down = np.interp(down_filter_t, filter_t, local_kernel)

            for j, stim_force in enumerate(this_rate_maps):
                this_stim_force = np.interp(down_t, this_extended_t, stim_force)
                local_signal = np.convolve(0.001 * down_dt * this_stim_force, local_kernel_down)[:len(down_t)] / \
                               saturation_factor * kernel_scale / filter_ratio
                max_local_signal = max(max_local_signal, np.max(local_signal))
                this_signal = np.minimum(local_signal, global_signal)
                this_area = np.trapz(this_signal, dx=down_dt)
                plasticity_signal[j] += this_area
                if plot and j == int(len(this_rate_maps)/2) and i == 0 and attempt == 1:
                    ylim = max(np.max(local_signal), np.max(global_signal))
                    x_start = 0.25 + this_extended_t[start_index]/this_interp_t[-1]/2.
                    x_end = 0.25 + this_extended_t[end_index]/this_interp_t[-1]/2.
                    fig, axes = plt.subplots(1)
                    axes.plot(down_t, local_signal, label='Local signal', color='g')
                    axes.plot(down_t, global_signal, label='Global signal', color='k')
                    axes.fill_between(down_t, 0., this_signal, label='Overlap', facecolor='r', alpha=0.5)
                    axes.axhline(y=ylim*1.05, xmin=x_start, xmax=x_end, linewidth=3, c='k')
                    axes.legend(loc='best', frameon=False, framealpha=0.5)
                    axes.set_xlabel('Time (ms)')
                    axes.set_ylabel('Signal amplitude (a.u.)')
                    axes.set_xlim(-0.5*this_interp_t[-1], 1.5*this_interp_t[-1])
                    axes.set_ylim(-0.05*ylim, ylim*1.1)
                    axes.set_title('Induced plasticity signal')
                    clean_axes(axes)
                    plt.show()
                    plt.close()
        saturation_factor *= filter_ratio * max_local_signal / max_global_signal
            # print 'Computed induction trial', i, 'in %i s' % (time.time() - start_time)

    if plot:
        x_start = np.mean(induction_locs[induction]) / track_length
        ylim = np.max(plasticity_signal) + 1.
        fig, axes = plt.subplots(1)
        axes.scatter(peak_locs['CA3'], plasticity_signal+1., color='k')
        axes.axhline(y=ylim + 0.1, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='k')
        axes.set_ylabel('Relative synaptic weight')
        axes.set_xlabel('Location (cm)')
        axes.set_xlim(0., track_length)
        axes.set_title('Induced synaptic weight distribution')
        clean_axes(axes)
        plt.show()
        plt.close()

    return plasticity_signal


def ramp_error1(x, xmin, xmax, ramp, induction=None, orig_weights=None, baseline=None, plot=False, full_output=False):
    """
    Calculates a rule_waveform and set of weights to match the first place field induction.
    :param x: array [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, kernel_scale]
    :param xmin: array
    :param xmax: array
    :param ramp: array
    :param induction: int: key for dicts of arrays
    :param orig_weights: array
    :param baseline: float
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
    if not check_bounds(x, xmin, xmax):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    elif (x[1] < x[0]) or (x[3] < x[2]):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    start_time = time.time()
    formatted_x = '[' + ', '.join(['%.3f' % xi for xi in x]) + ']'
    print 'Trying x: %s' % formatted_x
    if induction is None:
        induction = 1
    local_kernel, global_kernel = build_kernel1(x, plot)
    this_weights = calculate_plasticity_signal(x, local_kernel, global_kernel, induction, plot)
    # model_ramp_baseline = np.mean(get_expected_depolarization(default_rate_maps, np.ones_like(this_weights),
    #                                                   default_extended_t, default_interp_x))
    model_ramp = get_expected_depolarization(default_rate_maps, this_weights + 1., default_extended_t, default_interp_x)
    model_ramp_baseline = np.mean(model_ramp[ramp_baseline_indexes])
    model_ramp -= model_ramp_baseline

    Err = 0.
    for j in range(len(ramp)):
        Err += ((ramp[j] - model_ramp[j]) / 0.05) ** 2.

    if plot:
        x_start = np.mean(induction_locs[induction])/track_length
        ylim = max(np.max(ramp), np.max(model_ramp))
        ymin = min(np.min(ramp), np.min(model_ramp))
        fig, axes = plt.subplots(1)
        axes.plot(x_bins, ramp, label='Experiment', color='r')
        axes.plot(x_bins, model_ramp, label='Model', color='k')
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

    formatted_x = '[' + ', '.join(['%.4f' % xi for xi in x]) + ']'
    print 'x: %s, Err: %.4E took %i s' % (formatted_x, Err, time.time()-start_time)
    if full_output:
        return local_kernel, global_kernel, this_weights, model_ramp
    else:
        hist.x.append(x)
        hist.Err.append(Err)
        return Err


def optimize_polish(x, xmin, xmax, error_function, ramp, induction=None, orig_weights=None, baseline=None, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param ramp: array
    :param induction: int: key for dicts of arrays
    :param orig_weights: array
    :param baseline: float
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 400

    result = optimize.minimize(error_function, x, method='Nelder-Mead', options={'ftol': 1e-3,
                                                    'xtol': 1e-3, 'disp': True, 'maxiter': maxfev},
                               args=(xmin, xmax, ramp, induction, orig_weights, baseline))
    formatted_x = '['+', '.join(['%.3E' % xi for xi in result.x])+']'
    print 'Process: %i completed optimize_polish after %i iterations with Error: %.4E and x: %s' % (os.getpid(),
                                                                            result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


def optimize_explore(x, xmin, xmax, error_function, ramp, induction=None, orig_weights=None, baseline=None,
                     maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param ramp: array
    :param induction: int: key for dicts of arrays
    :param orig_weights: array
    :param baseline: float
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 700

    take_step = Normalized_Step(x, xmin, xmax)
    minimizer_kwargs = dict(method=null_minimizer, args=(xmin, xmax, ramp, induction, orig_weights, baseline))
    result = optimize.basinhopping(error_function, x, niter=maxfev, niter_success=maxfev/2,
                                   disp=True, interval=min(20, int(maxfev/20)), minimizer_kwargs=minimizer_kwargs,
                                   take_step=take_step)
    formatted_x = '['+', '.join(['%.3E' % xi for xi in result.x])+']'
    print 'Process: %i completed optimize_explore after %i iterations with Error: %.4E and x: %s' % (os.getpid(),
                                                                            result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


local_signal = {}
global_signal = {}
delta_weights = {}
model_ramp = {}

x0 = {}
# x0['1'] = [300., 6.69E+02, 9.77E-02, 6.40E+01, 1.31E+02, 4.04E-03]
# x0['1'] = [223.606, 978.505, 0.117, 60.137, 106.736, 0.003]
# x0['1'] = [2.34E+02, 8.74E+02, 1.53E-01, 5.33E+01, 1.09E+02, 2.21E-03]
# x0['1'] = [3.00E+02, 7.18E+02, 1.91E-01, 2.97E+01, 1.00E+02, 2.10E-03]
# x0['1'] = [1.00E+01, 1.35E+03, 3.02E-01, 2.12E+01, 2.15E+02, 3.13E-03]

# x0['1'] = [3.61E+01, 1.24E+03, 1.29E+01, 4.91E+02, 3.21E-03]  # Err: 1.929E+04
x0['1'] = [1.000E+01, 1.669E+03, 1.000E+01, 3.366E+02, 4.820E-03]  # 4.6118E+04
# x0['2'] = [4.96E+02, 2.34E+03, 1.00E+01, 2.55E+01, 1.96E-02]  # Err: 9.360E+04
# x0['2'] = [5.000E+02, 2.343E+03, 2.042E+01, 2.500E+01, 8.808E-03]  # Err: 9.5359E+04
x0['2'] = [5.000E+02, 2.365E+03, 1.002E+01, 2.500E+01, 9.188E-03]  # Err: 9.4305E+04
# x0['3'] = [1.00E+01, 1.79E+03, 3.55E+01, 1.22E+02, 3.13E-03]  # Err: 9.785E+03
# x0['3'] = [1.001E+01, 1.650E+03, 2.943E+01, 1.349E+02, 1.536E-03]  # Err: 1.0000E+04
x0['3'] = [1.000E+01, 1.632E+03, 5.000E+01, 1.064E+02, 1.553E-03]  # Err: 9.8435E+03
# x0['4'] = [1.00E+01, 2.13E+03, 1.83E+01, 5.00E+02, 1.46E-03]  # Err: 2.792E+04
# x0['4'] = [1.126E+01, 2.124E+03, 2.196E+01, 4.926E+02, 1.001E-03]  # Err: 2.9893E+04
x0['4'] = [1.000E+01, 3.419E+03, 5.000E+01, 9.351E+02, 1.212E-03]  # Err: 2.6862E+04
# x0['5'] = [3.59E+02, 6.97E+02, 2.07E+01, 3.76E+02, 2.40E-03]  # Err: 3.785E+03
x0['5'] = [3.150E+02, 7.208E+02, 1.000E+01, 3.844E+02, 2.241E-03]  # Err: 3.7821E+03
# x0['6'] = [1.003E+02, 5.002E+02, 1.695E+01, 3.845E+02, 2.663E-03]  # Err: 1.5555E+04
x0['6'] = [1.540E+02, 3.449E+02, 1.000E+01, 3.081E+02, 2.731E-03]  # Err: 1.3552E+04
# x0['7'] = [2.395E+02, 7.743E+02, 2.788E+01, 5.000E+02, 1.871E-03]  # Err: 9.7416E+04
x0['7'] = [3.280E+02, 1.231E+03, 1.000E+01, 9.999E+02, 2.056E-03]  # Err: 8.2503E+04
# x0['8'] = [2.877E+02, 6.819E+02, 1.223E+01, 5.000E+02, 2.538E-03]  # Err: 1.8180E+04
x0['8'] = [4.823E+02, 4.823E+02, 1.489E+01, 5.154E+02, 2.526E-03]  # Err: 1.7002E+04
# x0['9'] = [1.000E+01, 6.480E+02, 1.619E+01, 5.000E+02, 2.296E-03]  # Err: 8.4610E+04
x0['9'] = [2.635E+01, 1.040E+03, 1.000E+01, 9.995E+02, 2.486E-03]  # Err: 6.1448E+04
# x0['10'] = [4.956E+02, 3.000E+03, 1.000E+01, 2.500E+01, 5.693E-03]  # Err: 7.9027E+04
x0['10'] = [5.000E+02, 4.000E+03, 1.401E+01, 2.500E+01, 6.154E-03]  # Err: 6.1108E+04
x0['11'] = [2.46E+02, 5.00E+02, 4.97E+01, 4.23E+02, 1.00E-03]  # Err: 2.3930E+05

# to avoid saturation and reduce variability of time courses across cells, impose the relative amplitude
# of global and local kernels:
filter_ratio = 1.5
# [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, kernel_scale]
# x1 = [300., 1000., 50., 500., 2.e-3]
if cell_id in x0:
    x1 = x0[cell_id]
else:
    x1 = x0['1']
xmin1 = [10., 300., 10., 25., 5.e-4]
xmax1 = [500., 4000., 50., 1000., 2.e-2]


induction = 1

# ramp_error1(x1, xmin1, xmax1, ramp[induction], induction, plot=True)

result = optimize_explore(x1, xmin1, xmax1, ramp_error1, ramp[induction], induction, maxfev=700)

polished_result = optimize_polish(result['x'], xmin1, xmax1, ramp_error1, ramp[induction], induction)

hist.report_best()
hist.export('112716_magee_data_optimization_long_cell'+cell_id)
"""
ramp_error1(polished_result['x'], xmin1, xmax1, ramp[induction], induction, plot=True)

local_signal[induction], global_signal[induction], delta_weights[induction], model_ramp[induction] = \
    ramp_error1(x1, xmin1, xmax1, ramp[induction], induction, plot=True, full_output=True)



colors = ['r', 'k', 'c', 'g']
x_start = [induction_loc/track_length for induction_loc in induction_locs]

ylim = max(np.max(ramp), np.max(model_ramp))
fig, axes = plt.subplots(1)
for i in range(len(model_ramp)):
    axes.plot(x_bins, ramp[i], color=colors[i+2], label='Experiment'+str(i+1))
    axes.plot(x_bins, model_ramp[i], color=colors[i], label='Model'+str(i+1))
    axes.axhline(y=ylim+0.3, xmin=x_start[i], xmax=x_start[i]+0.02, c=colors[i], linewidth=3., zorder=0)
axes.set_xlabel('Location (cm)')
axes.set_xlim(0., track_length)
axes.set_ylabel('Depolarization (mV)')
plt.legend(loc='best', frameon=False, framealpha=0.5)
plt.show()
plt.close()

fig, axes = plt.subplots(1)
t_waveform = np.arange(-rule_max_timescale, rule_max_timescale, dt)
for i in range(len(model_ramp)):
    axes.plot(t_waveform, rule_waveforms[i], label='Induction kernel '+str(i+1), color=colors[i])
axes.set_xlabel('Time (ms)')
axes.set_xlim(-rule_max_timescale, rule_max_timescale)
axes.set_ylabel('Change in synaptic weight per spike')
plt.legend(loc='best', frameon=False, framealpha=0.5)
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
plt.legend(loc='best', frameon=False, framealpha=0.5)
plt.show()
plt.close()

if len(delta_weights) > 1:
    fig, axes = plt.subplots(1)
    axes.scatter(delta_weights[0]+1., delta_weights[1]+1.)
    clean_axes(axes)
    plt.show()
    plt.close()
"""
