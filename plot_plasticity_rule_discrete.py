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

# input_filename = 'output120520161957-pid18137-seed0-e1600-i600-induction_spine_voltage0'
input_filename = 'output120720161451-pid48517-seed0-e1600-i600-induction_spine_voltage0'

peak_locs = {}
spine_vm = {}
spike_trains = {}
successes = {}
dt = 1.

with h5py.File(data_dir+input_filename+'.hdf5', 'r') as f:
    trial = f.itervalues().next()
    induction_run_vel = trial.attrs['run_vel']
    track_length = trial.attrs['track_length']
    stim_dt = trial.attrs['stim_dt']
    equilibrate = trial.attrs['equilibrate']
    track_equilibrate = trial.attrs['track_equilibrate']
    duration = trial.attrs['duration']
    track_duration = duration - equilibrate - track_equilibrate
    interp_t = np.arange(0., duration, dt)
    track_t = np.arange(0., track_duration, dt)
    start = int((equilibrate+track_equilibrate)/dt)
    t = trial['time'][:]
    dx = dt * induction_run_vel / 1000.
    interp_x = np.arange(0., len(track_t)*dx, dx)[:len(track_t)]
    for rec in trial['rec'].itervalues():
        if 'spine' in rec.attrs['description']:
            index = rec.attrs['index']
            vm = np.interp(interp_t, t, rec[:])[start:]
            spine_vm[index] = vm
    for key in (key for key in trial['train'] if trial['train'][key].attrs['index'] in spine_vm.keys()):
        index = trial['train'][key].attrs['index']
        peak_locs[int(index)] = trial['train'][key].attrs['peak_loc']
        spike_trains[int(index)] = trial['train'][key][:]
        successes[int(index)] = trial['successes'][key][:]


field1_loc = 0.5
induction_dur = 300.


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


def get_saturation_factor(rule):
    """
    Given the local and global kernels, convolve each input spike train (successes) with the local kernel, and convolve
    the current injection with the global kernel. Plot the resulting plasticity signals for a single induction trial.
    :param x: array: [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio,
                    kernel_scale, saturation_factor]
    :param rule: str
    "param input_index: int
    """
    x = x0[rule]
    filter_ratio = x[4]
    kernel_scale = x[5]
    saturation_factor = 0.02
    this_local_kernel = local_kernel[rule]
    this_global_kernel = global_kernel[rule]

    induction_gate = np.zeros_like(track_t)
    start_index = np.where(interp_x >= induction_loc)[0][0]
    end_index = start_index + int(induction_dur / dt)
    induction_gate[start_index:end_index] = 1.
    global_signal = np.convolve(induction_gate, this_global_kernel)[:len(track_t)] * kernel_scale
    filter_t = np.arange(0., len(this_local_kernel) * dt, dt)
    for attempt in range(2):
        max_local_signal = []
        plasticity_signal = {}
        for input_index in spike_trains:
            success_train = successes[input_index]
            indexes = (np.array(success_train)/dt).astype(int)
            this_stim_force = np.zeros_like(track_t)
            this_stim_force[indexes] = 1.
            local_signal = np.convolve(this_stim_force, this_local_kernel)[:len(track_t)] / saturation_factor * kernel_scale / \
                           filter_ratio
            this_signal = np.minimum(local_signal, global_signal)
            this_area = np.trapz(this_signal, dx=dt)
            plasticity_signal[input_index] = this_area
            if attempt == 1:
                plt.plot(track_t, local_signal)
            max_local_signal.append(np.max(local_signal))
        if attempt == 0:
            saturation_factor *= filter_ratio * np.mean(max_local_signal) / np.max(global_signal)
    print 'Rule: %s, saturation_factor: %.3E' % (rule, saturation_factor)
    plt.plot(track_t, global_signal)
    plt.title('Plasticity signals (%s)' % rule)
    plt.show()
    plt.close()
    x_start = induction_loc / track_length
    ylim = np.max(plasticity_signal.values()) + 1.
    fig, axes = plt.subplots(1)
    axes.scatter(peak_locs.values(), np.add(plasticity_signal.values(), 1.), color='k')
    axes.axhline(y=ylim + 0.1, xmin=x_start, xmax=x_start + 0.02, linewidth=3, c='k')
    axes.set_ylabel('Relative synaptic weight')
    axes.set_xlabel('Location (cm)')
    axes.set_xlim(0., track_length)
    axes.set_title('Induced synaptic weight distribution (%s)' % rule)
    clean_axes(axes)
    plt.show()
    plt.close()


def plot_plasticity_signal_discrete(rule, input_index):
    """
    Given the local and global kernels, convolve each input spike train (successes) with the local kernel, and convolve
    the current injection with the global kernel. Plot the resulting plasticity signals for a single induction trial.
    :param rule: str
    "param input_index: int
    """
    x = x0[rule]
    filter_ratio = x[4]
    kernel_scale = x[5]
    saturation_factor = x[6]
    this_local_kernel = local_kernel[rule]
    this_global_kernel = global_kernel[rule]

    induction_gate = np.zeros_like(track_t)
    start_index = np.where(interp_x >= induction_loc)[0][0]
    end_index = start_index + int(induction_dur / dt)
    induction_gate[start_index:end_index] = 1.
    global_signal = np.convolve(induction_gate, this_global_kernel)[:len(track_t)] * kernel_scale
    filter_t = np.arange(0., len(this_local_kernel) * dt, dt)
    spike_train = spike_trains[input_index]
    success_train = successes[input_index]
    indexes = (np.array(success_train)/dt).astype(int)
    this_stim_force = np.zeros_like(track_t)
    this_stim_force[indexes] = 1.
    local_signal = np.convolve(this_stim_force, this_local_kernel)[:len(track_t)] / saturation_factor * kernel_scale / \
                   filter_ratio
    this_signal = np.minimum(local_signal, global_signal)
    x_start = interp_t[start_index] / track_duration
    x_end = interp_t[end_index] / track_duration
    fig, axes = plt.subplots(3, sharex=True)
    ymax0 = -30.  # np.max(spine_vm[input_index])
    ymin0 = np.min(spine_vm[input_index])
    ydepth0 = ymax0 - ymin0
    ymax1 = np.max(global_signal)
    axes[0].plot(track_t / 1000., spine_vm[input_index], label='Spine Vm', color='k')
    axes[0].scatter(spike_train / 1000., np.ones_like(spike_train) * (ymax0 - ydepth0 * 0.1), color='k')
    axes[0].scatter(success_train / 1000., np.ones_like(success_train) * (ymax0 - ydepth0 * 0.2), color='r')
    axes[0].set_ylabel('Spine voltage (mV)')
    axes[0].set_ylim((ymin0 - ydepth0 * 0.05), ymax0)
    axes[1].plot(track_t / 1000., local_signal, label='Local plasticity signal (%s)' % rule, color='r')
    axes[1].plot(track_t / 1000., global_signal, label='Global plasticity signal (%s)' % rule, color='k')
    axes[1].fill_between(track_t / 1000., 0., this_signal, label='Overlap', facecolor='r', alpha=0.5)
    axes[1].axhline(y=ymax1*0.95, xmin=x_start, xmax=x_end, linewidth=3, c='k')
    axes[1].legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel('Signal amplitude (a.u.)')
    axes[1].set_xlim(0., track_duration/1000.)
    axes[1].set_ylim(0., ymax1*1.05)
    # axes.set_title('Induced plasticity signal')
    clean_axes(axes)
    plt.show()
    plt.close()


# kernel_scale corrected to produce closer to 2.5 mean weight in center 10 spatial bins 6 mV with 30 cm/s during
# induction laps
# to avoid saturation and reduce variability of time courses across cells, constrain the relative amplitude
# of global and local kernels:
# [local_rise_tau, local_decay_tau, global_rise_tau, global_decay_tau, filter_ratio, kernel_scale, saturation_factor]
x0 = {}  # induced + spontaneous
x0['long'] = [2.614E+02, 1.103E+03, 2.686E+01, 4.574E+02, 1.388E+00, 4.533E-03*0.4619, 1.372E-02]
x0['short'] = [10., 100., 1.507E+01, 7.193E+01, 1.208E+00, 5.363E-03*1.0165, 2.716E-02]

"""
induction_loc = field1_loc*track_length
local_kernel, global_kernel = {}, {}
for rule in x0:
    local_kernel[rule], global_kernel[rule] = build_kernels(x0[rule], plot=False)


for rule in x0:
    get_saturation_factor(rule)

# for input_index in spine_vm:
for input_index in [18662]:
    # for rule in ['long']:
    for rule in x0:
        plot_plasticity_signal_discrete(rule, input_index)


short_weights_filename = '120716 short discrete weights'
long_weights_filename = '120716 long discrete weights'
peak_locs_filename = '120716 peak_locs discrete'
ramps_filename = '120716 plasticity discrete ramps'
ramp_x, long_ramp, short_ramp = read_from_pkl(data_dir+ramps_filename+'.pkl')
short_weights = read_from_pkl(data_dir+short_weights_filename+'.pkl')
long_weights = read_from_pkl(data_dir+long_weights_filename+'.pkl')
peak_locs_x = read_from_pkl(data_dir+peak_locs_filename+'.pkl')

x_start = induction_loc / track_length
start_index = np.where(interp_x >= induction_loc)[0][0]
end_index = start_index + int(induction_dur / dt)
x_end = interp_x[end_index] / track_length

ylim = max(np.max(short_weights), np.max(long_weights)) + 1.
fig, axes = plt.subplots(2, 2)
fig.set_size_inches(5.2, 3.9)
axes[0][0].scatter(peak_locs_x, short_weights+1., color='r', s=1.)
axes[0][0].scatter(peak_locs_x, long_weights+1., color='b', s=1.)
axes[0][0].axhline(y=ylim + 0.15, xmin=x_start, xmax=x_end, linewidth=2, c='k')
axes[0][0].set_ylabel('Relative synaptic weight')
axes[0][0].set_xlabel('Location (cm)')
axes[0][0].set_xlim(0., track_length)
axes[0][0].set_ylim(0.5, ylim + 0.3)
axes[0][0].set_title('Induced synaptic weight distribution', fontsize=12.)
clean_axes(axes[0])
plt.tight_layout()
plt.show()
plt.close()

for ramp in long_ramp, short_ramp:
    baseline = np.mean(ramp[np.where(ramp <= np.percentile(ramp, 10.))[0]])
    ramp -= baseline

ylim = max(np.max(long_ramp), np.max(short_ramp))
fig, axes = plt.subplots(2, 2)
fig.set_size_inches(5.2, 3.9)
axes[0][0].plot(ramp_x, short_ramp, color='r')
axes[0][0].plot(ramp_x, long_ramp, color='b')
axes[0][0].axhline(y=ylim + 0.25, xmin=x_start, xmax=x_end, linewidth=2, c='k')
axes[0][0].set_ylabel('Depolarization (mV)')
axes[0][0].set_xlabel('Location (cm)')
axes[0][0].set_xlim(0., track_length)
axes[0][0].set_ylim(-0.5, ylim + 0.5)
axes[0][0].set_title('Induced Vm ramp', fontsize=12.)
clean_axes(axes[0])
plt.tight_layout()
plt.show()
plt.close()
"""