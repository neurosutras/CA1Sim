__author__ = 'Grace Ng'
import parallel_optimize_spike_adaptation_engine
import os
import sys
from numpy import *
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import mkl

"""
Optimizes gCA factor and gCadepK factor, and gKm for soma and gKm factor for axon initial segment and hill

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

mkl.set_num_threads(1)

def spike_adaptation_error(x, target_spike_num=6):
    """

    :param spike_num: int
    :param axon_th: float
    :param start: float
    :param stop: float
    :return: float
    """
    if not check_bounds.within_bounds(x, 'spike_adaptation'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        hist.x_values.append(x)
        Err = 1e9
        hist.error_values.append(Err)
        return Err
    start_time = time.time()
    dv['x'] = x
    hist.x_values.append(x)
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    print 'Process %i using current x: %s: %s' % (os.getpid(), str(xlabels['spike_adaptation']), formatted_x)
    #rheobase: the current to cross threshold for a single spike; uses a 100 ms injection
    result = c[0].apply(parallel_optimize_spike_adaptation_engine.adjust_spike_number, 1)
    last = ''
    while not result.ready():
        time.sleep(1.)
        clear_output()
        stdout = result.stdout
        if stdout:
            line = stdout.splitlines()[-1]
            if line != last:
                print line
                last = line
        sys.stdout.flush()
    result = result.get()
    if result is None:
        print 'Cell is spontaneously firing, or parameters are out of bounds, or no reasonable rheobase.'
        Err = 1e9
        hist.error_values.append(Err)
        return Err
    final_result = result
    rheo_spike_times = result['spike_times']
    rheobase = result['amp']
    print 'rheobase: %.2f' % (rheobase)
    #Calculate more frequencies using stim duration of 500 ms
    result = v.map_async(parallel_optimize_spike_adaptation_engine.sim_spike_times,
                         [rheobase + amp_incr for amp_incr in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]])
    last = []
    while not result.ready():
        time.sleep(1.)
        clear_output()
        for i, stdout in enumerate([stdout for stdout in result.stdout if stdout][-3:]):
            line = stdout.splitlines()[-1]
            if line not in last:
                print line
                last.append(line)
        if len(last) > len(x):
            last = last[-len(x):]
        sys.stdout.flush()
    result = result.get()
    freq = []
    target_freq = []
    adi = -1.
    Err = 0.
    stim_dur = parallel_optimize_spike_adaptation_engine.stim_dur
    for trial in result:
        if trial is not None:
            spike_times = trial['spike_times']
            freq.append(len(spike_times)/stim_dur*1000.)
            if adi is -1.:
                if len(spike_times) >= target_spike_num:
                    adi = adapt_index(spike_times)
                    adi_error = ((target_val['adi'][str(target_spike_num)] - adi)/target_range['adi'][str(target_spike_num)])**2.
                    Err += adi_error
            target_freq.append(target_val['slope']*(trial['amp']-rheobase) + 1./stim_dur*1000.)
    FI_err = 0.
    for i in range(len(freq)):
        FI_err += ((target_freq[i] - freq[i])/(0.01*target_freq[i]))**2.
    FI_err *= FI_err_factor
    Err += FI_err
    for target in final_result:
        if target not in hist.features:
            hist.features[target] = []
        hist.features[target].append(final_result[target])
    print 'Simulation took %i s' % (time.time() - start_time)
    print 'Process : [soma.gCa factor', 'soma.gCadepK factor', 'soma.gKm', 'axon(ais, hill).gKm factor]: ' \
          '[%.2f, %.2f, %.2f, %.2f]' % (x[0], x[1], x[2], x[3])
    print 'Process %i: adi error %.4E, FI error %.4E, Total error: %.4E' % (os.getpid(), adi_error, FI_err, Err)
    hist.error_values.append(Err)
    sys.stdout.flush()
    return Err

def plot_best(x=None, discard=True):
    """
    Run simulations on the engines with the last best set of parameters, have the engines export their results to .hdf5,
    and then read in and plot the results.
    :param x: array
    """
    if x is None:
        if hist.x_values:
            dv['x'] = hist.report_best()
            c[0].apply(parallel_optimize_spike_adaptation_engine.adjust_spike_number, 6)
        else:
            raise Exception('Please specify input parameters (history might be empty).')
    else:
        dv['x'] = x
        c[0].apply(parallel_optimize_spike_adaptation_engine.adjust_spike_number, 6)
    dv.execute('export_sim_results()')
    # Does this export the most recent sim object? Does the sim object get re-written every time we run a new current injection?
    rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir + filename + '.hdf5')]
    for i, rec_filename in enumerate(rec_file_list):
        with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
            plt.figure(i)
            for trial in f.itervalues():
                amplitude = trial.attrs['amp']
                fig, axes = plt.subplots(1)
                for rec in trial['rec'].itervalues():
                    axes.plot(trial['time'], rec, label=rec.attrs['description'])
                axes.legend(loc='best', frameon=False, framealpha=0.5)
                axes.set_xlabel('Time (ms)')
                axes.set_ylabel('Vm (mV)')
                axes.set_title('Optimize Vm: I_inj amplitude %.2f' % amplitude)
                fig.tight_layout()
                clean_axes(axes)
    plt.show()
    plt.close()
    if discard:
        for rec_filename in rec_file_list:
            os.remove(data_dir + rec_filename + '.hdf5')

def adapt_index(spike_times):
    """

    :param spike_times: list of the times at which there are spike peaks
    :return: adi is a large value for high spike adaptation (large differences in lengths of interspike intervals)
    """
    if len(spike_times) < 4:
        return None
    adi = 0
    count = 0
    isi = []
    for i in range(len(spike_times) - 1):
        isi.append(spike_times[i + 1] - spike_times[i])
    for i in range(len(isi) - 1):
        adi += 1.0 * (isi[i + 1] - isi[i]) / (isi[i + 1] + isi[i])
        count += 1
    adi /= count
    return adi

#The target values and acceptable ranges
target_val = {}
target_range = {}

#This is the target adaptation index for a series of 6 spikes
target_val['adi'] = {'6': 0.118989}
target_val['slope'] = 50.
#50 spikes/s/nA
target_val['corr'] = 1.
target_range['adi'] = {'6': 0.001}
target_range['corr'] = 0.25
#Do we need to modify the target range??

FI_err_factor = 0.5

x0 = {}
xlabels = {}
xmin = {}
xmax = {}


if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    mech_filename = '041317 GC optimizing excitability'
if len(sys.argv) > 3:
    cluster_id = sys.argv[3]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

# These values will now be saved in the mech dictionary that is updated by previous round of optimization
x0['spike_adaptation'] = [1., 1., 0.0015, 5.]
xlabels['spike_adaptation'] = ['soma.gCa factor', 'soma.gCadepK factor', 'soma.gKm', 'axon(ais, hill).gKm factor']
xmin['spike_adaptation'] = [0.5, 0.5, 0.0005, 1.]
xmax['spike_adaptation'] = [2., 2., 0.003, 5.]

max_niter = 2100  # max number of iterations to run
niter_success = 400  # max number of interations without significant progress before aborting optimization
take_step = Normalized_Step(x0['spike_adaptation'], xmin['spike_adaptation'], xmax['spike_adaptation'])
minimizer_kwargs = dict(method=null_minimizer)

#Check bounds to make sure soma gKm and axon gKm are within a reasonable range
check_bounds = CheckBounds(xmin, xmax)
hist = optimize_history()
hist.xlabels = xlabels
history_filename = '040717 spike adaptation'

dv = c[:]
dv.block = True
global_start_time = time.time()


dv.execute('run parallel_optimize_spike_adaptation_engine %i \"%s\"' % (int(spines), mech_filename))
# time.sleep(60)
v = c.load_balanced_view()


result = optimize.basinhopping(spike_adaptation_error, x0['spike_adaptation'], niter=max_niter,
                               niter_success=niter_success, disp=True, interval=40,
                               minimizer_kwargs=minimizer_kwargs, take_step=take_step)
print result

"""
best_x = hist.report_best()
hist.export_to_pkl(history_filename)
#dv['x'] = best_x
dv['x'] = x0['spike_adaptation']
c[0].apply(parallel_optimize_spike_adaptation_engine.update_mech_dict)
plot_best(x0['spike_adaptation'])
"""