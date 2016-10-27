__author__ = 'Grace Ng'
import parallel_optimize_dendritic_excitability_engine
import sys
import os
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import scipy.optimize as optimize

"""
Aims for spike initiation at initial segment by increasing nax density and decreasing activation V1/2 relative to soma,
axon_hill, and axon compartments. Extend linear kap gradient into basals and obliques, aim for 60% spike attenuation
at bifurcation of trunk and tuft.

Hierarchical optimization:
I) optimize g_pas for target rinp at soma, trunk bifurcation, and tuft bifurcation [without h].
II) optimize ghbar_h for target rinp at soma, trunk bifurcation, and tuft bifurcation, while also optimizing for v_rest
offset between soma and tuft, and EPSP shape changes between proximal and distal synapses measured at the soma.
III) optimize gbar_nax/nas/sh/sha, gkabar_kap/d, gkdrbar for target na spike threshold, AHP amp, and vm stability

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""


def check_bounds(x, param_name):
    """
    Check that the current set of parameters for optimization are within the desired bounds.
    :param x: array
    :param param_name: str
    :return: bool
    """
    for i in range(len(x)):
        if ((xmin[param_name][i] is not None and x[i] < xmin[param_name][i]) or
                (xmax[param_name][i] is not None and x[i] > xmax[param_name][i])):
            return False
    return True


class History(object):
    def __init__(self):
        """

        """
        self.x_values = []
        self.error_values = []
        self.Rinp_values = {section: [] for section in ['soma', 'dend']}

    def report_best(self):
        """
        Report the input parameters and output values with the lowest error.
        """
        lowest_Err = min(self.error_values)
        index = self.error_values.index(lowest_Err)
        best_x = self.x_values[index]
        best_Rinp_values = {section: self.Rinp_values[section][index] for section in self.Rinp_values}
        formatted_x = '[' + ', '.join(['%.2E' % xi for xi in best_x]) + ']'
        print 'best x: %s' % formatted_x
        print 'lowest Err: %.3E' % lowest_Err
        print 'Rinp:', ['%s: %.1f' % (section, Rinp) for (section, Rinp) in best_Rinp_values.iteritems()]
        return best_x

    def export_to_pkl(self, hist_filename):
        """
        Save the history to .pkl
        :param hist_filename: str
        """
        saved_history = {'x_values': self.x_values, 'error_values': self.error_values, 'Rinp_values': self.Rinp_values}
        write_to_pkl(data_dir+hist_filename+'.pkl', saved_history)

    def import_from_pkl(self, hist_filename):
        """
        Update a history object with data from a .pkl file
        :param hist_filename: str
        """
        previous_history = read_from_pkl(data_dir+hist_filename +'.pkl')
        self.x_values = previous_history['x_values']
        self.error_values = previous_history['error_values']
        self.Rinp_values = previous_history['Rinp_values']

    def plot(self):
        """
        Remember to : also plot each value in x against error, and against input resistance
        """
        for section in self.Rinp_values:
            plt.plot(self.Rinp_values[section], self.error_values, label=section)
            plt.xlabel("Rinp values")
            plt.ylabel("Error values")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()


hist = History()


def pas_error(x):
    """
    Distribute simulations across available engines for optimization of leak conductance density gradient.
    :param x: array (soma.g_pas, dend.g_pas slope, dend.g_pas tau, dend.g_pas xhalf)
    :return: float
    """
    if not check_bounds(x, 'pas'):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    start_time = time.time()
    dv['x'] = x
    hist.x_values.append(x)

    sec_list = ['soma', 'dend']
    print 'Process %i using current x: [soma.g_pas, dend.g_pas slope, dend.g_pas tau, dend.g_pas xhalf]: ' \
          '[%.2E, %.2E, %.1f, %.1f]' % (os.getpid(), x[0], x[1], x[2], x[3])
    result = v.map_async(parallel_optimize_dendritic_excitability_engine.get_Rinp_for_section, sec_list)
    last_printed = ''
    while not result.ready():
        clear_output()
        for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2] and lines[-2] != last_printed:
                print lines[-2]
        sys.stdout.flush()
        time.sleep(1.)
    result = result.get()

    Err = 0.

    final_result = {}
    for dict in result:
        final_result.update(dict)
    for section in final_result:
        Err += ((target_val['pas'][section] - final_result[section]) / target_range['pas'][section]) ** 2.
        hist.Rinp_values[section].append(final_result[section])
    hist.error_values.append(Err)

    print('Simulation took %.3f s' % (time.time() - start_time))
    print 'Process %i: [soma.g_pas, dend.g_pas slope, dend.g_pas tau, dend.g_pas xhalf]: [%.2E, %.2E, %.1f, %.1f], ' \
          'soma R_inp: %.1f, dend R_inp: %.1f' % (os.getpid(), x[0], x[1], x[2], x[3], final_result['soma'],
                                                  final_result['dend'])
    print 'Error: %.3E' % Err
    return Err


def plot_best(x=None):
    """
    Run simulations on the engines with the last best set of parameters, have the engines export their results to .hdf5,
    and then read in and plot the results.
    :param x: array
    """
    if x is None:
        if hist.x_values:
            pas_error(hist.report_best())
        else:
            raise Exception('Please specify input parameters (history might be empty).')
    else:
        pas_error(x)
    dv.execute('export_sim_results()')
    rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir + filename + '.hdf5')]
    for rec_filename in rec_file_list:
        with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
            for trial in f.itervalues():
                target = trial.attrs['target']
                section = trial.attrs['section']
                optimization = trial.attrs['optimization']
                fig, axes = plt.subplots(1)
                for rec in trial['rec'].itervalues():
                    axes.plot(trial['time'], rec, label=rec.attrs['description'])
                axes.legend(loc='best', frameon=False, framealpha=0.5)
                axes.set_xlabel('Time (ms)')
                axes.set_ylabel('Vm (mV)')
                axes.set_title('Optimize %s: %s (%s)' % (optimization, target, section))
                clean_axes(axes)
                plt.show()
                plt.close()
    for rec_filename in rec_file_list:
        os.remove(data_dir + rec_filename + '.hdf5')


v_init = -67.
soma_ek = -77.

# the target values and acceptable ranges
target_val = {}
target_range = {}
target_val['pas'] = {'soma': 295., 'dend': 375.}
target_range['pas'] = {'soma': 10., 'dend': 10.}
target_val['v_rest'] = {'soma': v_init, 'tuft_offset': 0.}
target_range['v_rest'] = {'soma': 0.25, 'tuft_offset': 0.1}
target_val['na_ka'] = {'v_rest': v_init, 'th_v': -51., 'soma_peak': 40., 'trunk_amp': 0.6, 'ADP': 0., 'AHP': 4.,
                       'stability': 0., 'ais_delay': 0., 'slow_depo': 25.}
target_range['na_ka'] = {'v_rest': 0.25, 'th_v': .2, 'soma_peak': 2., 'trunk_amp': 0.01, 'ADP': 0.01, 'AHP': .2,
                         'stability': 1., 'ais_delay': 0.001, 'slow_depo': 1.}


x0 = {}
xmin = {}
xmax = {}

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    cluster_id = sys.argv[2]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

if spines:
    # x0['pas'] = [1.52E-06, 5.63E-04, 67., 220.]
    # x0['pas'] = [1.165e-06, 5.317e-04, 2.07e+02, 2.89e+01]
    # x0['pas'] = [1.85E-07, 5.15E-03, 2.13E+01, 3.85E+02]  # Err: 7.588
    x0['pas'] = [1.39E-07, 7.85E-03, 12.2, 370.8]  # Err: 7.510E-10
    xmin['pas'] = [1.0E-8, 1.0E-07, 5., 25.]
    xmax['pas'] = [2.0E-4, 2.0E-02, 400., 500.]
else:
    # x0['pas'] = [1.50E-06, 3.33E-04, 1.50E+02, 1.04E+02]
    # x0['pas'] = [1.66E-06, 3.83E-04, 1.71E+02, 1.69E+02]
    x0['pas'] = [2.80E-08, 5.25E-03, 9.46E+00, 3.58E+02]  # Err: 1.004E-07
    xmin['pas'] = [1.0E-9, 1.0E-07, 5., 25.]
    xmax['pas'] = [2.0E-4, 2.0E-02, 400., 500.]

# [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor]
x0['na_ka_stability'] = [0.0305, 0.0478, 2.97, 2.34]  # Error: 4.2416E+02
xmin['na_ka_stability'] = [0.01, 0.01, 1., 2.]
xmax['na_ka_stability'] = [0.05, 0.05, 5., 5.]

# [soma.sh_nas, trunk.ka factor]
# x0['na_ka_dend'] = [0., 2.9]
x0['na_ka_dend'] = [1.7, 2.8]
xmin['na_ka_dend'] = [0., 1.5]
xmax['na_ka_dend'] = [4., 5.]

# [ais.sha_nas, ais.gbar_nax factor]
x0['ais_delay'] = [-3.6, 4.4]

# [soma.e_pas, tuft.e_pas]
x0['v_rest'] = [-63., -77.]
xmin['v_rest'] = [v_init, soma_ek]
xmax['v_rest'] = [-63., -63.]

explore_niter = 1000  # max number of iterations to run
polish_niter = 200
take_step = Normalized_Step(x0['pas'], xmin['pas'], xmax['pas'])
minimizer_kwargs = dict(method=null_minimizer)

dv = c[:]
dv.clear()
dv.block = True
global_start_time = time.time()
dv.execute('run parallel_optimize_dendritic_excitability_engine %i' % int(spines))

v = c.load_balanced_view()

"""
result = optimize.basinhopping(pas_error, x0['pas'], niter=explore_niter, niter_success=400,
                               disp=True, interval=40, minimizer_kwargs=minimizer_kwargs, take_step=take_step)
print result


polished_result = optimize.minimize(pas_error, result.x, method='Nelder-Mead', options={'ftol': 1e-3,
                                                    'xtol': 1e-3, 'disp': True, 'maxiter': polish_niter})
print polished_result


polished_result = optimize.minimize(pas_error, x0['pas'], method='Nelder-Mead', options={'ftol': 1e-3,
                                                    'xtol': 1e-3, 'disp': True, 'maxiter': polish_niter})
print polished_result
"""