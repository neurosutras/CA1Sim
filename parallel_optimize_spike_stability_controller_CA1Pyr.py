__author__ = 'Grace Ng'
import parallel_optimize_spike_stability_engine_CA1Pyr
import sys
import os
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import scipy.optimize as optimize
# import mkl

"""
Aims for spike initiation at initial segment by increasing nax density and decreasing activation V1/2 relative to soma,
axon_hill, and axon compartments. Extend linear kap gradient into basals and obliques, aim for 60% spike attenuation
at bifurcation of trunk and tuft.

Optimizes gbar_nax/nas/sh/sha, gkabar_kap/d, gkdrbar for target na spike threshold, AHP amp, and vm stability

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

# mkl.set_num_threads(1)


def na_ka_stability_error(x, plot=0):
    """

    :param x: array [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor, soma.sh_nas/x,
                    axon.gkdrbar factor, dend.gkabar factor]
    :param plot: int
    :return: float
    """
    hist.x_values.append(x)
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    print 'Trying x: %s: %s' % (str(xlabels['na_ka_stability']), formatted_x)
    if not check_bounds.within_bounds(x, 'na_ka_stability'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        Err = 1e9
        hist.error_values.append(Err)
        return Err
    start_time = time.time()
    dv['x'] = x
    result = c[0].apply(parallel_optimize_spike_stability_engine_CA1Pyr.compute_spike_shape_features)
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
        print 'Cell is spontaneously firing, or parameters are out of bounds.'
        Err = 1e9
        hist.error_values.append(Err)
        return Err
    final_result = result
    rheobase = result['amp']

    result = v.map_async(parallel_optimize_spike_stability_engine_CA1Pyr.compute_spike_stability_features,
                         [[rheobase+0.05, 300.], [rheobase+0.5, 100.]])
    last_buffer_len = []
    while not result.ready():
        time.sleep(1.)
        clear_output()
        for i, stdout in enumerate(result.stdout):
            if (i + 1) > len(last_buffer_len):
                last_buffer_len.append(0)
            if stdout:
                lines = stdout.splitlines()
                if len(lines) > last_buffer_len[i]:
                    for line in lines[last_buffer_len[i]:]:
                        print line
                    last_buffer_len[i] = len(lines)
        sys.stdout.flush()
    result = result.get()
    temp_dict = {}
    temp_dict['amp'] = []
    temp_dict['rate'] = []
    for i, this_dict in enumerate(result):
        temp_dict['amp'].append(this_dict['amp'])
        temp_dict['rate'].append(this_dict['rate'])
        if 'stability' not in final_result:
            final_result['stability'] = this_dict['stability']
        else:
            final_result['stability'] += this_dict['stability']
        if 'slow_depo' not in final_result:
            final_result['slow_depo'] = this_dict['v_min_late'] - final_result['v_th']
        else:
            final_result['slow_depo'] += this_dict['v_min_late'] - final_result['v_th']
    indexes = range(len(temp_dict['rate']))
    indexes.sort(key=temp_dict['amp'].__getitem__)
    temp_dict['amp'] = map(temp_dict['amp'].__getitem__, indexes)
    temp_dict['rate'] = map(temp_dict['rate'].__getitem__, indexes)
    target_f_I = experimental_f_I_slope * np.log(temp_dict['amp'][0] / rheobase)
    final_result['rate'] = temp_dict['rate'][0]
    f_I_Err = ((temp_dict['rate'][0] - target_f_I) / (0.01 * target_f_I))**2.
    Err = f_I_Err
    for target in final_result:
        if target not in hist.features:
            hist.features[target] = []
        hist.features[target].append(final_result[target])
    for target in ['v_th', 'ADP', 'AHP', 'stability', 'slow_depo', 'dend_amp']:
        # don't penalize AHP or slow_depo less than target
        if not ((target == 'AHP' and final_result[target] < target_val['na_ka'][target]) or
                (target == 'slow_depo' and final_result[target] < target_val['na_ka'][target])):
            Err += ((target_val['na_ka'][target] - final_result[target])/target_range['na_ka'][target])**2.
            if target not in hist.features:
                hist.features[target] = []
            hist.features[target].append(final_result[target])
    if 'rate' not in hist.features:
        hist.features['rate'] = []
    hist.features['rate'].append(final_result['rate'])

    print 'Simulation took %i s' % (time.time()-start_time)
    print 'Process %i: [soma.gkabar, soma.gkdrbar, soma.sh_nas/x, axon.gkdrbar factor, dend.gkabar factor, ' \
          'soma.gCa factor, soma.gCadepK factor, soma.gkmbar]: %s, amp: %.3f, v_rest: %.1f, threshold: %.1f, ' \
          'ADP: %.1f, AHP: %.1f, stability: %.2f, slow_depo: %.2f, dend_amp: %.2f, rate %.3f' % \
          (os.getpid(), formatted_x, final_result['amp'], final_result['v_rest'], final_result['v_th'],
           final_result['ADP'], final_result['AHP'], final_result['stability'], final_result['slow_depo'],
           final_result['dend_amp'], final_result['rate'])
    print 'Process %i: f_I Error: %.4E, Error: %.4E' % (os.getpid(), f_I_Err, Err)
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
            na_ka_stability_error(hist.report_best())
        else:
            raise Exception('Please specify input parameters (history might be empty).')
    else:
        na_ka_stability_error(x)
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
                axes.set_title('Optimize spike stability: I_inj amplitude %.2f' % amplitude)
                fig.tight_layout()
                clean_axes(axes)
    plt.show()
    plt.close()
    if discard:
        for rec_filename in rec_file_list:
            os.remove(data_dir + rec_filename + '.hdf5')


v_init = -77.
soma_ek = -77.

# the target values and acceptable ranges
target_val = {}
target_range = {}
target_val['na_ka'] = {'v_rest': v_init, 'v_th': -51., 'soma_peak': 40., 'ADP': 0., 'AHP': 3.,
                       'stability': 0., 'ais_delay': 0., 'slow_depo': 20., 'dend_amp': 0.6}
target_range['na_ka'] = {'v_rest': 0.25, 'v_th': .05, 'soma_peak': 2., 'ADP': 0.01, 'AHP': .01,
                         'stability': 1., 'ais_delay': 0.001, 'slow_depo': 0.5, 'dend_amp': 0.005}

experimental_f_I_slope = 12.  # Hz/ln(pA); rate = slope * ln(current - rheobase)
# CA1Pyr experimental f-I data from Kowalski J...Pernia-Andrade AJ, Hippocampus, 2016

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
    mech_filename = '041817 CA1Pyr optimizing spike stability'
if len(sys.argv) > 3:
    cluster_id = sys.argv[3]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

check_bounds = CheckBounds(xmin, xmax)
xlabels['na_ka_stability'] = ['soma.gkabar', 'soma.gkdrbar', 'soma.sh_nas/x', 'axon.gkdrbar factor',
                              'dend.gkabar factor', 'soma.gCa factor', 'soma.gCadepK factor', 'soma.gkmbar']
hist = optimize_history()
hist.xlabels = xlabels['na_ka_stability']


# [soma.gkabar, soma.gkdrbar, soma.sh_nas/x, axon.gkdrbar factor, dend.gkabar factor,
#            'soma.gCa factor', 'soma.gCadepK factor', 'soma.gkmbar']

# x0['na_ka_stability'] = [0.0305, 0.0478, 1.7, 1., 2.8, 2.006E+00, 2.995E+00, 0.0015]
x0['na_ka_stability'] = [1.136E-02, 3.995E-02, 3.968E+00, 1.323E+00, 3.007E+00, 1.357E+00, 4.207E+00, 3.880E-03]
# Err: 1.360E+03
xmin['na_ka_stability'] = [0.01, 0.01, 0.1, 1., 1., 1., 1., 0.0005]
xmax['na_ka_stability'] = [0.05, 0.05, 6., 2., 5., 5., 5., 0.005]

max_niter = 1500  # max number of iterations to run
ninterval = max_niter / 50

take_step = Normalized_Step(x0['na_ka_stability'], xmin['na_ka_stability'], xmax['na_ka_stability'])
minimizer_kwargs = dict(method=null_minimizer)

dv = c[:]
dv.block = True
global_start_time = time.time()


dv.execute('run parallel_optimize_spike_stability_engine_CA1Pyr %i \"%s\"' % (int(spines), mech_filename))
# time.sleep(60)
v = c.load_balanced_view()
"""
result = optimize.basinhopping(na_ka_stability_error, x0['na_ka_stability'], niter=max_niter,
                               niter_success=niter_success, disp=True, interval=ninterval,
                               minimizer_kwargs=minimizer_kwargs, take_step=take_step)
print result

history_filename = '041417 spike stability optimization history'
best_x = hist.report_best()
# hist.export_to_pkl(history_filename)
dv['x'] = best_x
# dv['x'] = x0['na_ka_stability']
c[0].apply(parallel_optimize_spike_stability_engine_CA1Pyr.update_mech_dict)
"""
# plot_best(best_x)
plot_best(x0['na_ka_stability'])