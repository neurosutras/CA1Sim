__author__ = 'milsteina'
from function_lib import  *
from ipyparallel import interactive

"""
These methods attempt to fit the experimental BTSP data with a simple 3-state markov-like process model with 
nonstationary rates dependent on the presence of two intracellular biochemical signals.

Assumptions:
1) Synaptic weights are all = 1 prior to induction 1.
"""

try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass


if len(sys.argv) > 1:
    cell_id = str(sys.argv[1])
else:
    cell_id = None

experimental_filenames = {'cell': '121216 magee lab first induction', 'spont_cell': '120216 magee lab spont'}

if len(sys.argv) > 2:
    label = str(sys.argv[2])
else:
    label = 'cell'

if label in ['cell', 'spont_cell']:
    experimental_filename = experimental_filenames[label]
else:
    raise Exception('Unknown label or category of data.')

experimental_file_dir = data_dir

# placeholder for parameters pushed from controller
# [local_decay_tau, global_decay_tau, local_kon, global_kon, global_koff, saturated_delta_weights]
x = [6.318E+02, 1.707E+02, 3.604E-01, 2.497E-01, 1.020E-01, 2.]
induction = 1


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


def generate_spatial_rate_maps():
    """
    Given a set of place field peak locations, return firing rate vs. location computed at a resolution of
    track_length/10000 bins.
    :return: array
    """
    spatial_rate_maps = []
    for i, peak_loc in enumerate(peak_locs):
        gauss_force = excitatory_peak_rate * np.exp(-((extended_x - peak_loc) / gauss_sigma) ** 2.)
        gauss_force = wrap_around_and_compress(gauss_force, generic_x)
        spatial_rate_maps.append(gauss_force)
    return spatial_rate_maps


def generate_complete_rate_maps(induction, spatial_rate_maps):
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
    for i, rate_map in enumerate(rate_maps):
        this_EPSP_map = np.interp(default_interp_x, this_interp_x, rate_map) * EPSP_scaling_factor
        this_EPSP_map = np.concatenate([this_EPSP_map for i in range(3)])
        this_EPSP_map = np.convolve(this_EPSP_map, epsp_filter)[:3 * len(default_interp_x)]
        EPSP_maps.append(this_EPSP_map[len(default_interp_x):2 * len(default_interp_x)])
    return np.array(EPSP_maps)


dt = 1.  # ms
down_dt = 10.  # ms; lower temporal resolution to speed up calculation during optimization
equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # ms
input_field_width = 90.  # cm
track_length = 187.  # cm

spatial_resolution = 1.  # CA3 place field inputs per cm of linear track
# generates a predicted 6 mV depolarization from gaussian weights with peak = 2.5
EPSP_scaling_factor = 3.7635E-03  # may need to be re-determined empirically when spatial resolution changes

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
        padded_t = np.insert(this_interp_t, 0, this_interp_t[-vel_window_bins-1:-1] - this_interp_t[-1])
        padded_t = np.append(padded_t, this_interp_t[1:vel_window_bins+1] + this_interp_t[-1])
        padded_x = np.insert(this_interp_x, 0, this_interp_x[-vel_window_bins-1:-1] - track_length)
        padded_x = np.append(padded_x, this_interp_x[1:vel_window_bins+1] + track_length)
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

excitatory_peak_rate = 40.
v_init = -67.

local_random = random.Random()

gauss_sigma = input_field_width / 3. / np.sqrt(2.)  # contains 99.7% gaussian area

num_inputs = int(track_length * spatial_resolution)
peak_locs = np.linspace(0., track_length, num_inputs+1)[:-1]

spatial_rate_maps = generate_spatial_rate_maps()  # x=generic_x

complete_t, complete_x, complete_rate_maps, complete_induction_gates = {}, {}, {}, {}
for induction in position:
    complete_t[induction], complete_x[induction], complete_rate_maps[induction] = \
        generate_complete_rate_maps(induction, spatial_rate_maps)
    complete_induction_gates[induction] = generate_complete_induction_gate(induction)

input_matrix = compute_EPSP_matrix(spatial_rate_maps, generic_x)  # x=default_interp_x

baseline = None
for induction in position:
    if baseline is None:
        ramp_baseline_indexes = np.where(np.array(ramp[induction]) <= np.percentile(ramp[induction], 10.))[0]
        baseline = np.mean(ramp[induction][ramp_baseline_indexes])
    ignore = subtract_baseline(ramp[induction], baseline)

print 'Process: %i ready to optimize %s %s' % (os.getpid(), label, cell_id)


@interactive
def build_kernels(x, plot=False):
    """
    Construct two kernels with exponential decay:
    1) Local kernel that generates a plasticity signal at each synapse
    2) Global kernal that generates a plasticity signal during dendritic calcium spiking
    :param x: array: [local_decay_tau, global_decay_tau]
    :param plot: bool
    :return: array, array
    """
    local_decay_tau = x[0]
    global_decay_tau = x[1]

    max_time_scale = max(local_decay_tau, global_decay_tau)
    filter_t = np.arange(0., 7.*max_time_scale, dt)
    local_filter = np.exp(-filter_t/local_decay_tau)
    decay_indexes = np.where(local_filter < 0.001)[0]
    if np.any(decay_indexes):
        local_filter = local_filter[:decay_indexes[0]]
    local_filter /= np.sum(local_filter)
    global_filter = np.exp(-filter_t / global_decay_tau)
    decay_indexes = np.where(global_filter < 0.001)[0]
    if np.any(decay_indexes):
        global_filter = global_filter[:decay_indexes[0]]
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


@interactive
def calculate_delta_weights(x, local_kernel, global_kernel, complete_rate_maps, induction, plot=False,
                            full_output=False):
    """
    Given a set of local and global signal kernels, convolve each input rate_map with the local kernel, and convolve the
    current injection with the global kernel. The weight change for each input results from a markov-like process with
    nonstationary transition rates.
    :param x: array: [local_decay_tau, global_decay_tau, local_kon, global_kon, global_koff, saturated_delta_weights]
    :param local_kernel: array
    :param global_kernel: array
    :param complete_rate_maps: dict of array
    :param induction: int: key for dicts of arrays
    :param plot: bool
    :param full_output: bool
    :return: delta_weights: array
    """
    down_sample = True
    local_kon = x[2]
    global_kon = x[3]
    global_koff = x[4]
    saturated_delta_weights = x[5]
    delta_weights = np.zeros_like(peak_locs)
    global_signal = np.convolve(complete_induction_gates[induction], global_kernel)[:len(complete_t[induction])]
    if down_sample:
        this_t = np.arange(complete_t[induction][0], complete_t[induction][-1] + down_dt / 2., down_dt)
        global_signal = np.interp(this_t, complete_t[induction], global_signal)
        this_dt = down_dt
    else:
        this_t = complete_t[induction]
        this_dt = dt
    states = {'U': 1., 'A': 0., 'P': 0.}
    rates = {'U': {'A': global_kon * global_signal},
            'A': {'U': global_koff}}
    state_machine = StateMachine(complete_t[induction][0], dt=this_dt, states=states, rates=rates)
    for j, stim_force in enumerate(complete_rate_maps[induction]):
        local_time = time.time()
        local_signal = np.convolve(0.001 * dt * stim_force, local_kernel)[:len(stim_force)]
        if down_sample:
            local_signal = np.interp(this_t, complete_t[induction], local_signal)
        # print 'Input: %i; convolution with local kernel took %.1f s' % (j, time.time() - local_time)
        local_time = time.time()
        state_machine.update_transition('A', 'P', local_kon * local_signal)
        state_machine.run()
        # print 'Input: %i; running state machine took %.1f s' % (j, time.time() - local_time)
        delta_weights[j] += saturated_delta_weights * state_machine.states['P']
        if plot and j == int(len(complete_rate_maps[induction])/2):
            buffer = 5000.
            # buffer = 0.
            # orig_font_size = mpl.rcParams['font.size']
            # orig_fig_size = mpl.rcParams['figure.figsize']
            # mpl.rcParams['font.size'] = 8.
            # mpl.rcParams['figure.figsize'] = 7.34, 3.25
            # fig1 = plt.figure()
            # gs1 = gridspec.GridSpec(2, 2)
            # axes = plt.subplot(gs1[0, 0])
            fig1, axes = plt.subplots(2, sharex=True)
            ylim = max(np.max(local_signal), np.max(global_signal))
            start_index = np.where(interp_x[induction][0] >= induction_locs[induction][0])[0][0]
            this_induction_start = interp_t[induction][0][start_index]
            this_induction_dur = induction_durs[induction][0]
            start_time = -buffer
            end_time = interp_t[induction][0][-1] + buffer
            this_duration = end_time - start_time
            x_start = (buffer + this_induction_start) / this_duration
            x_end = (buffer + this_induction_start + this_induction_dur) / this_duration
            axes[0].plot(this_t / 1000., global_signal, label='Global signal', color='b')
            axes[0].plot(this_t / 1000., local_signal, label='Local signal', color='k')
            axes[0].axhline(y=ylim*1.05, xmin=x_start, xmax=x_end, linewidth=1, c='k')
            axes[0].legend(loc='best', frameon=False, framealpha=0.5)
            axes[1].set_xlabel('Time (s)')
            axes[0].set_ylabel('Signal amplitude (a.u.)')
            axes[0].set_xlim(-buffer/1000., interp_t[induction][0][-1]/1000. + buffer/1000.)
            axes[0].set_ylim(-0.05*ylim, ylim*1.1)
            axes[0].set_title('Plasticity signals')
            axes[1].plot(state_machine.t_history/1000., state_machine.states_history['P'], color='r')
            axes[1].set_ylabel('Change in synaptic weight')
            clean_axes(axes)
            # gs1.tight_layout(fig1)
            fig1.tight_layout()
            plt.show()
            plt.close()
            # mpl.rcParams['font.size'] = orig_font_size
            # mpl.rcParams['figure.figsize'] = orig_fig_size
        # print 'Input: %i, numerical calculation of synaptic weight took %i s' % (j, time.time() - local_time)

    if full_output:
        state_machine_output = {}
        state_machine_output['time'] = np.array(state_machine.t_history)
        for state in state_machine.states:
            state_machine_output[state] = np.array(state_machine.states_history[state])
    else:
        state_machine_output = None
    return delta_weights, state_machine_output


@interactive
def ramp_error_parametric(local_x=None, local_induction=None, baseline=None, plot=False, full_output=False):
    """
    Given a set of model parameters, and run velocities during plasticity induction, this method calculates a set of 
    synaptic weights consistent with the provided experimental ramp data.
    :param local_x: array [local_decay_tau, global_decay_tau, local_kon, global_kon, global_koff, 
                            saturated_delta_weights]
    :param induction: int: key for dicts of arrays
    :param baseline: float
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
    if local_x is None:
        local_x = x
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in local_x]) + ']'
    if local_induction is None:
        local_induction = induction
    this_induction_loc = np.mean([induction_loc for induction_loc in induction_locs[local_induction] if
                                  induction_loc is not None])
    # print 'Process: %i trying x: %s for %s %s, induction_loc: %.1f' % (os.getpid(), formatted_x, label, cell_id,
                                                                       # this_induction_loc)
    sys.stdout.flush()
    exp_ramp = np.array(ramp[local_induction])
    start_time = time.time()
    local_kernel, global_kernel = build_kernels(x, plot)
    delta_weights, state_machine_output = calculate_delta_weights(local_x, local_kernel, global_kernel, complete_rate_maps,
                                                                  induction, plot, full_output)
    # print 'Process: %i calculated synaptic weights in %i s' % (os.getpid(), time.time() - start_time)
    sys.stdout.flush()
    amp, width, peak_shift, ratio, start_loc, end_loc = {}, {}, {}, {}, {}, {}
    amp['exp'], width['exp'], peak_shift['exp'], ratio['exp'], start_loc['exp'], end_loc['exp'] = \
        calculate_ramp_features(exp_ramp, this_induction_loc)

    weights = np.add(delta_weights, 1.)
    model_ramp = weights.dot(input_matrix)
    if baseline is None:
        model_baseline = subtract_baseline(model_ramp)
    else:
        model_baseline = baseline
        model_ramp -= model_baseline
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
        axes.plot(binned_x, model_ramp, label='Model Parametric)', color='k')
        axes.axhline(y=ylim + 0.2, xmin=x_start, xmax=x_start + 0.02, linewidth=1, c='k')
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
        axes1.plot(peak_locs, weights, c='r')
        axes1.axhline(y=ylim + 0.2, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='k')
        axes1.set_xlabel('Location (cm)')
        axes1.set_ylabel('Synaptic weights')
        axes1.set_title('Induction %i: Synaptic weights' % induction)
        axes1.set_xlim([0., track_length])
        axes1.set_ylim([math.floor(ymin), max(math.ceil(ylim), ylim + 0.4)])
        clean_axes(axes1)
        fig1.tight_layout()
        plt.show()
        plt.close()

    print 'Process: %i, %s %s, x: %s, Err: %.4E took %i s' % (os.getpid(), label, cell_id, formatted_x, Err,
                                                              time.time()-start_time)
    print 'exp: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, end_loc: %.1f' % \
          (amp['exp'], width['exp'], peak_shift['exp'], ratio['exp'], start_loc['exp'], end_loc['exp'])
    print 'model: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, end_loc: %.1f' % \
          (amp['model'], width['model'], peak_shift['model'], ratio['model'], start_loc['model'], end_loc['model'])
    sys.stdout.flush()
    hist.x.append(local_x)
    hist.Err.append(Err)
    if full_output:
        return local_kernel, global_kernel, weights, model_ramp, model_baseline, Err, state_machine_output
    else:
        return Err


def estimate_weights_nonparametric(ramp, input_matrix, induction=None, baseline=None, beta=10., plot=False,
                                   full_output=False):
    """
    Uses singular value decomposition to estimate a set of weights to match any arbitrary place field ramp, agnostic
    about underlying kernel, induction velocity, etc.
    :param ramp: dict of array
    :param input_matrix: array; x=default_interp_x
    :param induction: int: key for dicts of arrays
    :param baseline: float
    :param beta: regularization parameter
    :param plot: bool
    :param full_output: bool: whether to return all relevant objects (True), or just Err (False)
    :return: float
    """
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
    forced_weights = tuning_amp * np.cos(2. * np.pi / (input_field_width * 1.2) *
                                              (peak_locs - modulated_field_center)) + tuning_offset
    left = np.where(peak_locs >= modulated_field_center - input_field_width * 1.2 / 2.)[0][0]
    right = np.where(peak_locs > modulated_field_center + input_field_width * 1.2 / 2.)[0][0]
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

    D[np.where(np.eye(*D.shape))] = s / (s ** 2. + beta ** 2.)
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
    asymmetry = np.zeros_like(peak_locs)
    start_index = int(start_loc['exp']/track_length*len(peak_locs))
    induction_index = int(this_induction_loc/track_length*len(peak_locs))
    end_index = int(end_loc['exp']/track_length*len(peak_locs))
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
        axes1.plot(peak_locs, weights, c='r')
        axes1.axhline(y=ylim + 0.2, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='k')
        axes1.set_xlabel('Location (cm)')
        axes1.set_ylabel('Candidate synaptic weights')
        axes1.set_xlim([0., track_length])
        axes1.set_ylim([math.floor(ymin), max(math.ceil(ylim), ylim + 0.4)])
        clean_axes(axes1)
        fig1.tight_layout()
        plt.show()
        plt.close()

    print 'Process: %i, %s %s: SVD took %i s, Err: %.4E' % (os.getpid(), label, cell_id, time.time()-start_time, Err)
    print 'exp: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, end_loc: %.1f' % \
          (amp['exp'], width['exp'], peak_shift['exp'], ratio['exp'], start_loc['exp'], end_loc['exp'])
    print 'model: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, end_loc: %.1f' % \
          (amp['model'], width['model'], peak_shift['model'], ratio['model'], start_loc['model'], end_loc['model'])
    sys.stdout.flush()
    if full_output:
        return weights, model_ramp, model_baseline, asymmetry
    else:
        return Err


@interactive
def process_local(output_filename=None, plot=False, local_x=None):
    """
    
    :return: 
    """
    if local_x is None:
        x = x
    else:
        x = local_x

    local_kernel = {}
    global_kernel = {}
    weights_parametric = {}
    model_ramp_parametric = {}
    weights_SVD = {}
    model_ramp_SVD = {}
    model_baseline = {}
    asymmetry = {}
    state_machine_output = {}

    mean_induction_loc, mean_induction_dur = {}, {}

    for induction in position:
        if induction == 2 and 1 in position:
            this_model_baseline = model_baseline[1]
        else:
            this_model_baseline = None
        local_kernel[induction], global_kernel[induction], weights_parametric[induction], \
        model_ramp_parametric[induction], model_baseline[induction], state_machine_output[induction], Err = \
            ramp_error_parametric(induction=induction, baseline=this_model_baseline, plot=plot, full_output=True)
        mean_induction_loc[induction] = np.mean([induction_loc for induction_loc in induction_locs[induction] if
                                                 induction_loc is not None])
        mean_induction_dur[induction] = np.mean([induction_dur for induction_dur in induction_durs[induction] if
                                                 induction_dur is not None])

    for induction in ramp:
        if induction == 2:
            this_model_baseline = model_baseline[1]
        else:
            this_model_baseline = None
        weights_SVD[induction], model_ramp_SVD[induction], model_baseline[induction], asymmetry[induction] = \
            estimate_weights_nonparametric(ramp, input_matrix, induction, plot=False, full_output=True)

    if plot:
        fig1, axes1 = plt.subplots(1)
        fig2, axes2 = plt.subplots(1)
        fig3, axes3 = plt.subplots(1)

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
            axes1.plot(binned_x, model_ramp_parametric[induction], label='Model (Parametric)', c=colors[1])
            axes1.plot(binned_x, model_ramp_SVD[induction], label='Model (SVD)', c=colors[2])
            axes1.axhline(y=ylim1 + 0.25, xmin=x_start, xmax=x_end, linewidth=2, c='k')
            axes2.plot(peak_locs, weights_parametric[induction], label='Model (Parametric)', c=colors[1])
            axes2.plot(peak_locs, weights_SVD[induction], label='Model (SVD)', c=colors[2])
            axes2.axhline(y=ylim2 + 0.25, xmin=x_start, xmax=x_end, linewidth=2, c='k')
            for state in (state for state in state_machine_output[induction] if state != 'time'):
                axes3.plot(state_machine_output['time']/1000., state_machine_output[state])
        axes1.set_xlabel('Location (cm)')
        axes1.set_ylabel('Ramp depolarization (mV)')
        axes1.set_xlim([0., track_length])
        axes1.legend(loc='best', frameon=False, framealpha=0.5)
        axes1.set_title('Induction %s: Vm ramp' % induction)
        clean_axes(axes1)
        fig1.tight_layout()
        axes2.set_xlabel('Location (cm)')
        axes2.set_ylabel('Synaptic weights (a.u.)')
        axes2.set_title('Induction %s: Synaptic weights' % induction)
        axes2.legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes2)
        fig2.tight_layout()
        axes3.set_xlabel('Time (s)')
        axes3.set_ylabel('Occupancy')
        axes3.set_title('Induction %s: State occupancy' % induction)
        axes2.legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes3)
        fig3.tight_layout()

        plt.show()
        plt.close()

    if output_filename is not None:
        with h5py.File(data_dir + output_filename + '.hdf5', 'a') as f:
            for induction in position:
                if label not in f:
                    f.create_group(label)
                if 'position' not in f:
                    f.create_dataset('position', compression='gzip', compression_opts=9, data=binned_x)
                if 'peak_locs' not in f:
                    f.create_dataset('peak_locs', compression='gzip', compression_opts=9, data=peak_locs)
                if cell_id not in f[label]:
                    f[label].create_group(cell_id)
                if str(induction) not in f[label][cell_id]:
                    f[label][cell_id].create_group(induction)
                f[label][cell_id][induction].attrs['track_length'] = track_length
                f[label][cell_id][induction].attrs['induction_loc'] = mean_induction_loc[induction]
                f[label][cell_id][induction].attrs['induction_dur'] = mean_induction_dur[induction]
                f[label][cell_id][induction].attrs['parameters'] = x
                f[label][cell_id][induction].attrs['error'] = Err
                f[label][cell_id][induction].create_dataset('local_kernel', compression='gzip', compression_opts=9,
                                                  data=local_kernel[induction])
                f[label][cell_id][induction].create_dataset('global_kernel', compression='gzip', compression_opts=9,
                                                  data=global_kernel[induction])
                f[label][cell_id][induction].attrs['dt'] = dt
                f[label][cell_id][induction].create_dataset('exp_ramp', compression='gzip', compression_opts=9,
                                                            data=ramp[induction])
                f[label][cell_id][induction].create_dataset('model_ramp_parametric', compression='gzip',
                                                            compression_opts=9, data=model_ramp_parametric[induction])
                f[label][cell_id][induction].create_dataset('model_ramp_SVD', compression='gzip', compression_opts=9,
                                                  data=model_ramp_SVD[induction])
                f[label][cell_id][induction].create_dataset('weights_parametric', compression='gzip',
                                                            compression_opts=9, data=weights_parametric[induction])
                f[label][cell_id][induction].create_dataset('weights_SVD', compression='gzip', compression_opts=9,
                                                  data=weights_SVD[induction])
                if 'states' not in f[label][cell_id][induction]:
                    f[label][cell_id][induction].create_group('states')
                for key, value in state_machine_output[induction].iteritems():
                    f[label][cell_id][induction]['states'].create_dataset(key, compression='gzip', compression_opts=9,
                                                  data=value)

