__author__ = 'Grace Ng'
import parallel_optimize_leak_engine
import sys
import os
import math
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
from specify_cells2 import *
import scipy.optimize as optimize
import mkl

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

mkl.set_num_threads(1)

neurotree_filename = '121516_DGC_trees.pkl'
neurotree_dict = read_from_pkl(morph_dir+neurotree_filename)


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
        self.xlabels = []
        self.x_values = []
        self.error_values = []
        self.Rinp_values = {}

    def report_best(self):
        """
        Report the input parameters and output values with the lowest error.
        """
        lowest_Err = min(self.error_values)
        index = self.error_values.index(lowest_Err)
        best_x = self.x_values[index]
        best_Rinp_values = {section: self.Rinp_values[section][index] for section in self.Rinp_values}
        formatted_x = '[' + ', '.join(['%.3E' % xi for xi in best_x]) + ']'
        print 'best x: %s' % formatted_x
        print 'lowest Err: %.3E' % lowest_Err
        print 'Rinp:', ['%s: %.1f' % (section, Rinp) for (section, Rinp) in best_Rinp_values.iteritems()]
        return best_x

    def export_to_pkl(self, hist_filename):
        """
        Save the history to .pkl
        :param hist_filename: str
        """
        saved_history = {'xlabels': self.xlabels, 'x_values': self.x_values, 'error_values': self.error_values,
                         'Rinp_values': self.Rinp_values}
        write_to_pkl(data_dir+hist_filename+'.pkl', saved_history)

    def import_from_pkl(self, hist_filename):
        """
        Update a history object with data from a .pkl file
        :param hist_filename: str
        """
        previous_history = read_from_pkl(data_dir+hist_filename +'.pkl')
        self.xlabels = previous_history['xlabels']
        self.x_values = previous_history['x_values']
        self.error_values = previous_history['error_values']
        self.Rinp_values = previous_history['Rinp_values']

    def plot(self):
        """
        Remember to : also plot each value in x against error, and against input resistance
        """
        num_x_param = len(self.xlabels)
        num_plot_rows = math.floor(math.sqrt(num_x_param))
        print(num_plot_rows)
        num_plot_cols = math.ceil(num_x_param/num_plot_rows)
        print(num_plot_cols)

        #plot x-values against error
        plt.figure(1)
        for i, x_param in enumerate(self.xlabels):
            plt.subplot(num_plot_rows, num_plot_cols, i+1)
            x_param_vals = [x_val[i] for x_val in self.x_values]
            range_param_vals = max(x_param_vals) - min(x_param_vals)
            plt.scatter(x_param_vals, self.error_values)
            plt.xlim((min(x_param_vals)-0.1*range_param_vals, max(x_param_vals)+0.1*range_param_vals))
            plt.xlabel(x_param)
            plt.ylabel("Error values")
        #plt.show()
        #plt.close()

        plt.figure(2)
        colors = ['r', 'g', 'b', 'gray', 'darkviolet', 'goldenrod']
        #plot x-values against input resistance
        for i, x_param in enumerate(self.xlabels):
            plt.subplot(num_plot_rows, num_plot_cols, i+1)
            x_param_vals = [x_val[i] for x_val in self.x_values]
            range_param_vals = max(x_param_vals) - min(x_param_vals)
            for j, section in enumerate(self.Rinp_values):
                plt.scatter([x_param_vals], self.Rinp_values[section], label=section, color = colors[j])
            plt.xlim((min(x_param_vals) - 0.1 * range_param_vals, max(x_param_vals) + 0.1 * range_param_vals))
            plt.xlabel(x_param)
            plt.ylabel("Rinp values")
            plt.legend(loc='upper right', scatterpoints = 1, frameon=False, framealpha=0.5)
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

    sec_list = ['soma', 'dend', 'distal_dend']
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    print 'Process %i using current x: %s: %s' % (os.getpid(), str(xlabels['pas']), formatted_x)
    result = v.map_async(parallel_optimize_leak_engine.get_Rinp_for_section, sec_list)
    last = []
    while not result.ready():
        time.sleep(1.)
        clear_output()
        for i, stdout in enumerate([stdout for stdout in result.stdout if stdout][-len(sec_list):]):
            line = stdout.splitlines()[-1]
            if line not in last:
                print line
                last.append(line)
        if len(last) > len(sec_list):
            last = last[-len(sec_list):]
        sys.stdout.flush()
    result = result.get()

    Err = 0.

    final_result = {}
    for dict in result:
        final_result.update(dict)
    for section in final_result:
        if section not in hist.Rinp_values:
            hist.Rinp_values[section] = []
    for section in target_val['pas']:
        Err += ((target_val['pas'][section] - final_result[section]) / target_range['pas'][section]) ** 2.
        hist.Rinp_values[section].append(final_result[section])
    section = 'distal_dend'
    hist.Rinp_values[section].append(final_result[section])
    # add catch for decreasing terminal end input resistance too much
    if final_result['distal_dend'] < final_result['dend']:
        Err += ((final_result['dend'] - final_result['distal_dend']) / target_range['pas'][section]) ** 2.
    hist.error_values.append(Err)

    print('Simulation took %.3f s' % (time.time() - start_time))
    print 'Process %i: %s: %s; soma R_inp: %.1f, dend R_inp: %.1f, distal_dend R_inp: %.1f; Err: %.3E' % (os.getpid(),
                                                                                str(xlabels['pas']), formatted_x,
                                                                                final_result['soma'],
                                                                                final_result['dend'],
                                                                                final_result['distal_dend'], Err)
    return Err


def plot_best(x=None, discard=True):
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


def update_pas_exp(x):
    """

    x0 = [2.28e-05, 1.58e-06, 58.4]
    :param x: array [soma.g_pas, dend.g_pas slope, dend.g_pas tau]
    """
    if spines is False:
        cell.reinit_mechanisms(reset_cable=True)
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.modify_mech_param('apical', 'pas', 'g', origin='soma', slope=x[1], tau=x[2])
    for sec_type in ['basal', 'axon_hill', 'axon', 'ais', 'trunk', 'apical', 'tuft', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')
    if spines is False:
        cell.correct_for_spines()


v_init = -67.
soma_ek = -77.

# the target values and acceptable ranges
target_val = {}
target_range = {}
target_val['pas'] = {'soma': 295., 'dend': 375.}
target_range['pas'] = {'soma': 0.5, 'dend': 1.}
target_val['v_rest'] = {'soma': v_init, 'tuft_offset': 0.}
target_range['v_rest'] = {'soma': 0.25, 'tuft_offset': 0.1}

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
    # mech_filename = '121516 DG_GC pas spines'
    mech_filename = '012416 GC optimizing excitability'
if len(sys.argv) > 3:
    cluster_id = sys.argv[3]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

# xlabels['pas'] = ['soma.g_pas', 'dend.g_pas slope', 'dend.g_pas tau']
xlabels['pas'] = ['soma.g_pas', 'dend.g_pas slope', 'dend.g_pas tau']

if spines:
    x0['pas'] = [3.80E-08, 8.08E-07, 6.78E+01]  # Err: 2.527E-09
    xmin['pas'] = [1.0E-18, 1.0E-12, 25.]
    xmax['pas'] = [1.0E-7, 1.0E-4, 400.]
else:
    # x0['pas'] = [4.94E-08, 3.74E-06, 9.67E+01]  # Err: 3.995E-11
    x0['pas'] = [1.88486629e-07, 1.15951655e-06, 8.73328512e+01]  # Error: 2.231E+04
    xmin['pas'] = [1.0E-18, 1.0E-12, 25.]
    xmax['pas'] = [1.0E-6, 1.0E-4, 400.]

hist.xlabels = xlabels['pas']

explore_niter = 700  # max number of iterations to run
polish_niter = 400
take_step = Normalized_Step(x0['pas'], xmin['pas'], xmax['pas'])
minimizer_kwargs = dict(method=null_minimizer)

dv = c[:]
dv.clear()
dv.block = True
global_start_time = time.time()
dv.execute('run parallel_optimize_leak_engine %i \"%s\"' % (int(spines), mech_filename))
# time.sleep(60)
v = c.load_balanced_view()


result = optimize.basinhopping(pas_error, x0['pas'], niter=explore_niter, niter_success=explore_niter,
                               disp=True, interval=20, minimizer_kwargs=minimizer_kwargs, take_step=take_step)

polished_result = optimize.minimize(pas_error, result.x, method='Nelder-Mead', options={'ftol': 1e-5,
                                                    'disp': True, 'maxiter': polish_niter})
"""

polished_result = optimize.minimize(pas_error, x0['pas'], method='Nelder-Mead', options={'ftol': 1e-5, 'disp': True,
                                                                                         'maxiter': polish_niter})
"""
print polished_result
best_x = hist.report_best()
cell = DG_GC(neurotree_dict=neurotree_dict[0], mech_filename=mech_filename, full_spines=spines)
update_pas_exp(best_x)
# update_pas_exp(x0['pas'])
cell.export_mech_dict(cell.mech_filename)
