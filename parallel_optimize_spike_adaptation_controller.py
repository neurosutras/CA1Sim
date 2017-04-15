__author__ = 'Grace Ng'
import parallel_optimize_spike_adaptation_engine
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
# import mkl

"""
Optimizes gCA factor and gCadepK factor, and gkmbar for soma and axon

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

# mkl.set_num_threads(1)


def get_adaptation_index(spike_times):
    """
    A large value indicates large degree of spike adaptation (large increases in interspike intervals during a train)
    :param spike_times: list of float
    :return: float
    """
    if len(spike_times) < 3:
        return None
    isi = []
    adi = []
    for i in range(len(spike_times) - 1):
        isi.append(spike_times[i + 1] - spike_times[i])
    for i in range(len(isi) - 1):
        adi.append((isi[i + 1] - isi[i]) / (isi[i + 1] + isi[i]))
    return np.mean(adi)


def get_adaptation_index_error(spike_times):
    """

    :param spike_times: array of float 
    :return: float 
    """
    if len(spike_times) < 3:
        return None, None
    elif len(spike_times) > len(experimental_spike_times):
        adi = get_adaptation_index(spike_times[:len(experimental_spike_times)])
        exp_adi = experimental_adaptation_indexes[len(experimental_spike_times) - 3]
    else:
        adi = get_adaptation_index(spike_times)
        exp_adi = experimental_adaptation_indexes[len(spike_times) - 3]
    adi_error = ((adi - exp_adi) / (0.01 * exp_adi)) ** 2.
    return adi, adi_error


def spike_adaptation_error(x, full_output=False):
    """
    
    :param x: array
    :param full_output: bool
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
    formatted_x = '[' + ', '.join(['%.4E' % xi for xi in x]) + ']'
    print 'Controller: using current x: %s: %s' % (str(xlabels['spike_adaptation']), formatted_x)
    # rheobase: the current to cross threshold for a single spike; uses a 100 ms injection
    result = c[0].apply(parallel_optimize_spike_adaptation_engine.get_rheobase)
    last_buffer_len = 0
    while not result.ready():
        time.sleep(1.)
        clear_output()
        stdout = result.stdout
        if stdout:
            lines = stdout.splitlines()
            if len(lines) > last_buffer_len:
                for line in lines[last_buffer_len:]:
                    print line
                last_buffer_len = len(lines)
        sys.stdout.flush()
    result = result.get()
    if result is None:
        print 'Cell is spontaneously firing, or rheobase is outside target range.'
        Err = 1e9
        hist.error_values.append(Err)
        return Err
    rheobase = result
    final_result = {'rheobase': rheobase}
    # Calculate firing rates for a range of I_inj amplitudes using a stim duration of 500 ms
    result = v.map_async(parallel_optimize_spike_adaptation_engine.sim_f_I,
                         [rheobase + i_inj_increment * (i + 1) for i in range(num_increments)])
    last_buffer_len = []
    while not result.ready():
        time.sleep(1.)
        clear_output()
        for i, stdout in enumerate(result.stdout):
            if i + 1 > len(last_buffer_len):
                last_buffer_len.append(0)
            if stdout:
                lines = stdout.splitlines()
                if len(lines) > last_buffer_len[i]:
                    for line in lines[last_buffer_len[i]:]:
                        print line
                    last_buffer_len[i] = len(lines)
        sys.stdout.flush()
    result = result.get()
    final_result['amp'] = []
    final_result['adi'] = []
    final_result['f_I'] = []
    stim_dur = parallel_optimize_spike_adaptation_engine.stim_dur
    adi_Err = 0.
    for i, trial in enumerate(result):
        final_result['amp'].append(trial['amp'])
        spike_times = trial['spike_times']
        this_adi, this_adi_Err = get_adaptation_index_error(spike_times)
        if this_adi_Err is not None:
            adi_Err += this_adi_Err
        final_result['adi'].append(this_adi)
        this_rate = len(spike_times) / stim_dur * 1000.
        final_result['f_I'].append(this_rate)
    indexes = range(len(final_result['f_I']))
    indexes.sort(key=final_result['amp'].__getitem__)
    final_result['amp'] = map(final_result['amp'].__getitem__, indexes)
    final_result['adi'] = map(final_result['adi'].__getitem__, indexes)
    final_result['f_I'] = map(final_result['f_I'].__getitem__, indexes)
    f_I_Err = 0.
    for i, this_rate in enumerate(final_result['f_I']):
        f_I_Err += ((this_rate - target_f_I[i]) / (0.01 * target_f_I[i])) ** 2.
    Err = adi_Err + f_I_Err
    for feature in final_result:
        if feature not in hist.features:
            hist.features[feature] = []
        hist.features[feature].append(final_result[feature])
    print 'Simulation took %i s' % (time.time() - start_time)
    print '[soma.gCa factor, soma.gCadepK factor, soma.gkmbar]: %s' % formatted_x
    print 'Process %i: adi error %.4E, f_I error %.4E, Total error: %.4E' % (os.getpid(), adi_Err, f_I_Err, Err)
    hist.error_values.append(Err)
    sys.stdout.flush()
    if full_output:
        return Err, final_result
    else:
        return Err


def plot_best(x=None, discard=True):
    """
    Run simulations on the engines with the last best set of parameters, have the engines export their results to .hdf5,
    and then read in and plot the results.
    :param x: array
    """
    if x is None:
        if hist.x_values:
            spike_adaptation_error(hist.report_best())
        else:
            raise Exception('Please specify input parameters (history might be empty).')
    else:
        spike_adaptation_error(x)
    dv.execute('export_sim_results()')
    rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir + filename + '.hdf5')]
    for i, rec_filename in enumerate(rec_file_list):
        with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
            for trial in f.itervalues():
                amplitude = trial.attrs['amp']
                fig, axes = plt.subplots(1)
                for rec in trial['rec'].itervalues():
                    axes.plot(trial['time'], rec, label=rec.attrs['description'])
                axes.legend(loc='best', frameon=False, framealpha=0.5)
                axes.set_xlabel('Time (ms)')
                axes.set_ylabel('Vm (mV)')
                axes.set_title('Optimize f_I and spike adaptation: I_inj amplitude %.2f' % amplitude)
                fig.tight_layout()
                clean_axes(axes)
    plt.show()
    plt.close()
    if discard:
        for rec_filename in rec_file_list:
            os.remove(data_dir + rec_filename + '.hdf5')


# GC experimental spike adaptation data from Brenner...Aldrich, Nat. Neurosci., 2005
experimental_spike_times = [0., 8.57331572, 21.79656539, 39.24702774, 60.92470277, 83.34214003, 109.5640687,
                            137.1598415, 165.7067371, 199.8546896, 236.2219287, 274.3857332, 314.2404227, 355.2575958,
                            395.8520476, 436.7635403]
experimental_adaptation_indexes = []
for i in range(3, len(experimental_spike_times)+1):
    experimental_adaptation_indexes.append(get_adaptation_index(experimental_spike_times[:i]))

# 50 spikes/s/nA
experimental_f_I_slope = 50.
i_inj_increment = 0.1
num_increments = 8
target_f_I = [experimental_f_I_slope * i_inj_increment * (i + 1) for i in range(num_increments)]

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
    mech_filename = '041417 GC optimizing excitability'
if len(sys.argv) > 3:
    cluster_id = sys.argv[3]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

# These values will now be saved in the mech dictionary that is updated by previous round of optimization
x0['spike_adaptation'] = [1., 1., 0.0015]
xlabels['spike_adaptation'] = ['soma.gCa factor', 'soma.gCadepK factor', 'soma.gkmbar']
xmin['spike_adaptation'] = [0.5, 0.5, 0.0005]
xmax['spike_adaptation'] = [2., 2., 0.003]

max_niter = 2100  # max number of iterations to run
niter_success = 400  # max number of interations without significant progress before aborting optimization
take_step = Normalized_Step(x0['spike_adaptation'], xmin['spike_adaptation'], xmax['spike_adaptation'])
minimizer_kwargs = dict(method=null_minimizer)

# Check bounds to make sure soma gKm and axon gKm are within a reasonable range
check_bounds = CheckBounds(xmin, xmax)
hist = optimize_history()
hist.xlabels = xlabels
history_filename = '041417 spike adaptation'

dv = c[:]
dv.block = True
global_start_time = time.time()


dv.execute('run parallel_optimize_spike_adaptation_engine %i \"%s\"' % (int(spines), mech_filename))
# time.sleep(60)
v = c.load_balanced_view()

# Err, final_result = spike_adaptation_error(x0['spike_adaptation'], True)


result = optimize.basinhopping(spike_adaptation_error, x0['spike_adaptation'], niter=max_niter,
                               niter_success=niter_success, disp=True, interval=40,
                               minimizer_kwargs=minimizer_kwargs, take_step=take_step)
print result
"""

best_x = hist.report_best()
hist.export_to_pkl(history_filename)
dv['x'] = best_x
# dv['x'] = x0['spike_adaptation']
c[0].apply(parallel_optimize_spike_adaptation_engine.update_mech_dict)
"""
# plot_best(x0['spike_adaptation'])
