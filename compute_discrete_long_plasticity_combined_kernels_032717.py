__author__ = 'milsteina'
from plot_results import *
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

if len(sys.argv) > 1:
    group = str(sys.argv[1])
else:
    group = 'long'
if len(sys.argv) > 2:
    cell_id = str(sys.argv[2])
else:
    cell_id = None

experimental_file_dir = data_dir
if group == 'long':
    experimental_filename = '121216 magee lab first induction'
elif group == 'long_spont':
    experimental_filename = '120216 magee lab spont'

mkl.set_num_threads(2)


def generate_complete_induction_gate(induction):
    """

    :param induction: int
    :return:
    """
    absolute_induction_start_times = []
    absolute_induction_end_times = []
    running_delta_t = 0.
    for i in range(len(position[induction])):
        this_interp_t = interp_t[induction][i]
        this_interp_x = interp_x[induction][i]
        if i == 0:
            complete_t = this_interp_t - len(this_interp_t) * dt
            complete_x = this_interp_x - track_length
            complete_induction_gate = np.zeros_like(this_interp_t)
        this_induction_loc = induction_locs[induction][i]
        this_induction_dur = induction_durs[induction][i]
        this_induction_gate = np.zeros_like(this_interp_t)
        start_index = np.where(this_interp_x >= this_induction_loc)[0][0]
        end_index = start_index + int(this_induction_dur / dt)
        this_induction_gate[start_index:end_index] = 1.
        complete_t = np.append(complete_t, this_interp_t + running_delta_t)
        complete_x = np.append(complete_x, this_interp_x + i * track_length)
        complete_induction_gate = np.append(complete_induction_gate, this_induction_gate)
        absolute_induction_start_times.append(this_interp_t[start_index] + running_delta_t)
        absolute_induction_end_times.append(this_interp_t[start_index] + running_delta_t)
        running_delta_t += len(this_interp_t) * dt
        if i == len(position[induction]) - 1:
            complete_t = np.append(complete_t, this_interp_t + running_delta_t)
            complete_x = np.append(complete_x, this_interp_x + (i + 1) * track_length)
            complete_induction_gate = np.append(complete_induction_gate, np.zeros_like(this_interp_t))
    return complete_induction_gate, complete_t, complete_x, absolute_induction_start_times, absolute_induction_end_times


dt = 1.  # ms

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

back, forward = 6000., 6000.  # minimum interval between plateaus across all cells
kernel_t = np.arange(-back, forward, dt)

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

complete_t, complete_x, complete_induction_gates, absolute_induction_start_times, \
    absolute_induction_end_times = {}, {}, {}, {}, {}
for induction in [1]:  # position:
    complete_induction_gates[induction], complete_t[induction], complete_x[induction], \
        absolute_induction_start_times[induction], absolute_induction_end_times[induction] = \
        generate_complete_induction_gate(induction)

output_filename = '032617 discrete plasticity summary2'
with h5py.File(data_dir+output_filename+'.hdf5', 'r') as f:
    peak_locs = f['peak_locs'][:]
    weights = f[group][cell_id]['weights_parametric'][:]


def calculate_plasticity_rule(induction=None, plot=False):
    """
    Given the local and global kernels, convolve a single spike with the local kernel, and convolve the
    current injection with the global kernel. Compute the STDP-like kernel for each relative spike time. Weight changes
    are proportional to the area under the overlap of the two signals.
    :param induction: int
    :param plot: bool
    :return: array
    """
    if induction is None:
        induction = 1
    t, x, start_times = complete_t[induction], complete_x[induction], absolute_induction_start_times[induction]

    peak_times = []

    delta_locs = []

    adjusted_weights = []

    if plot:
        fig, axes = plt.subplots(2, sharex=True)

    for trial in range(len(start_times)):
    # trial = 0
        peak_times.append([])
        delta_locs.append([])
        adjusted_weights.append([])
        start_index = np.where(t >= start_times[trial])[0][0]
        start_t = t[start_index]
        start_x = x[start_index]
        for i, input in enumerate(peak_locs):
            delta_loc = input + track_length * trial - start_x
            delta_locs[trial].append(delta_loc)
            delta_index = np.where(x >= delta_loc + start_x)[0][0]
            delta_t = t[delta_index] - start_t
            peak_times[trial].append(delta_t)
            adjusted_weights[trial].append(weights[i])
            if delta_loc <= 0.:
                delta_loc += track_length
            else:
                delta_loc -= track_length
            delta_locs[trial].append(delta_loc)
            delta_index = np.where(x >= delta_loc + start_x)[0][0]
            delta_t = t[delta_index] - start_t
            peak_times[trial].append(delta_t)
            adjusted_weights[trial].append(weights[i])
        indexes = range(len(delta_locs[trial]))
        indexes.sort(key=delta_locs[trial].__getitem__)
        delta_locs[trial] = map(delta_locs[trial].__getitem__, indexes)
        peak_times[trial] = map(peak_times[trial].__getitem__, indexes)
        adjusted_weights[trial] = map(adjusted_weights[trial].__getitem__, indexes)
        if plot:
            axes[0].plot(np.array(peak_times[trial])/1000., delta_locs[trial], color='grey')
    if plot:
        axes[0].plot(np.mean(peak_times, axis=0)/1000., np.mean(delta_locs, axis=0))
        axes[1].plot(np.mean(peak_times, axis=0)/1000., np.mean(adjusted_weights, axis=0))
        axes[1].set_xlabel('time from induction (s)')
        axes[0].set_ylabel('position (cm)')
        axes[1].set_ylabel('synaptic weight (a.u.)')
        clean_axes(axes)
        plt.show()
        plt.close()

    return np.mean(delta_locs, axis=0), np.mean(adjusted_weights, axis=0), np.mean(peak_times, axis=0)


induction = 1
shifted_peak_locs, shifted_weights, relative_peak_times = calculate_plasticity_rule(induction, plot=True)

"""
output_filename = '032617 discrete plasticity summary2'
with h5py.File(data_dir+output_filename+'.hdf5', 'a') as f:
    f[group][cell_id].create_dataset('shifted_peak_locs', compression='gzip', compression_opts=9,
                                      data=shifted_peak_locs)
    f[group][cell_id].create_dataset('shifted_weights', compression='gzip', compression_opts=9,
                                     data=shifted_weights)
    f[group][cell_id].create_dataset('relative_peak_times', compression='gzip', compression_opts=9,
                                     data=relative_peak_times)
print 'Exported data for %s cell %s' % (group, cell_id)
"""