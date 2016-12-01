__author__ = 'milsteina'
from function_lib import *
from ipyparallel import interactive
import os

"""
In this version of the simulation, phase precession of CA3 inputs is implemented using the method from Chadwick et al.,
Elife, 2015, which uses a circular gaussian with a phase sensitivity factor that effectively compresses the range of
phases within each theta cycle that each input is active, which will reduce jitter across within-cycle input sequences.

"""

# [local_rise_tau, local_decay_tau, local_scale, global_rise_tau, global_decay_tau, global_scale]
x = None
local_kernel = None
global_kernel = None
induction = 1
peak_locs = {}
rate_maps = {}
induction_locs = {}
induction_durs = {}
extended_x = {}
extended_t = {}
dt = None
down_dt = None

@interactive
def calculate_plasticity_signal_piece(trial):
    """
    Given the local and global kernels, convolve each input rate_map with the local kernel, and convolve the
    current injection with the global kernel. The weight change for each input is proportional to the area under the
    product of the two signals. Incremental weight changes accrue across multiple induction trials.
    :param trial: int
    :return: plasticity_signal: array
    """
    local_scale = x[2]
    global_scale = x[5]
    group = 'CA3'
    plasticity_signal = np.zeros_like(peak_locs[group])
    i = trial
    start_time = time.time()
    this_rate_maps = rate_maps[induction][i]
    this_induction_loc = induction_locs[induction][i]
    this_induction_dur = induction_durs[induction][i]
    this_extended_x = extended_x[induction][i]
    this_extended_t = extended_t[induction][i]
    global_signal = np.zeros_like(this_extended_x)
    start_index = np.where(this_extended_x >= this_induction_loc)[0][0]
    end_index = start_index + int(this_induction_dur/dt)
    global_signal[start_index:end_index] = 1.
    global_signal = np.convolve(global_signal, global_kernel)[:len(this_extended_x)] * global_scale
    down_t = np.arange(this_extended_t[0], this_extended_t[-1] + down_dt / 2., down_dt)
    global_signal = np.interp(down_t, this_extended_t, global_signal)
    filter_t = np.arange(0., len(local_kernel)*dt, dt)
    down_filter_t = np.arange(0., filter_t[-1]+down_dt /2., down_dt)
    this_local_kernel = np.interp(down_filter_t, filter_t, local_kernel)

    for j, stim_force in enumerate(this_rate_maps):
        # time0 = time.time()
        this_stim_force = np.interp(down_t, this_extended_t, stim_force)
        local_signal = np.convolve(0.001 * down_dt * this_stim_force, this_local_kernel)[:len(down_t)] * local_scale
        # time1 = time.time()
        # print 'local signal: %.4f s' % (time1 - time0)
        this_signal = np.minimum(local_signal, global_signal)
        # time2 = time.time()
        # print 'this signal: %.4f s' % (time2 - time1)
        this_area = np.trapz(this_signal, dx=down_dt)
        # time3 = time.time()
        # print 'local signal: %.4f s' % (time3 - time2)
        plasticity_signal[j] += this_area

    # print 'Process:', os.getpid(), 'computed induction trial', trial, 'in %i s' % (time.time() - start_time)
    return plasticity_signal