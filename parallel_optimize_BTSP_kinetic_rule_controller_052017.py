__author__ = 'milsteina'
from function_lib import *
import parallel_optimize_BTSP_kinetic_rule_engine_052017
from ipyparallel import Client
from IPython.display import clear_output

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
    experimental_file_dir = sys.argv[1]
else:
    experimental_file_dir = data_dir

if len(sys.argv) > 2:
    cluster_id = sys.argv[2]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

experimental_filenames = {'cell': '121216 magee lab first induction', 'spont_cell': '120216 magee lab spont'}

cell_ids = []
labels = []

for label, experimental_filename in experimental_filenames.iteritems():
    with h5py.File(experimental_file_dir+experimental_filename+'.hdf5') as f:
        cell_ids.extend(f.keys())
        labels.extend([label for i in range(len(f))])

dv = c[:]
dv.block = True
global_start_time = time.time()

for i, engine in enumerate(c):
    engine.execute('run parallel_optimize_BTSP_kinetic_rule_engine_052017 \"%s\" \"%s\"' % (cell_ids[i], labels[i]))
time.sleep(60)

dv['experimental_file_dir'] = experimental_file_dir


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


def ramp_population_error(x, xmin, xmax, induction):
    """
    Push the specified parameters to the engines. Each engine processes a different cell from an experimental dataset.
    Calculate the mean error across engines.
    :param x: array [local_decay_tau, global_decay_tau, local_kon, global_kon, global_koff, saturated_delta_weights]
    :param xmin: array
    :param xmax: array
    :return: float
    """
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    if induction is None:
        induction = 1
    print 'Controller: trying x: %s' % formatted_x
    if not check_bounds(x, xmin, xmax):
        print 'Aborting: Invalid parameter values.'
        hist.x.append(x)
        hist.Err.append(1e9)
        return 1e9
    start_time = time.time()
    dv['x'] = x
    dv['induction'] = induction
    result = c[:].apply(parallel_optimize_BTSP_kinetic_rule_engine_052017.ramp_error_parametric)
    last_buffer_len = []
    ready = False
    while not ready:
        try:
            ready = result.ready()
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
        except:
            ready = True
    result = result.get()
    Err = np.mean(result)
    print result
    print 'Controller: Mean error: %.4E; calculation of error across all cells took %i s' % (Err,
                                                                                             time.time() - start_time)
    sys.stdout.flush()
    hist.x.append(x)
    hist.Err.append(Err)
    return Err


def optimize_polish(x, xmin, xmax, error_function, induction=None, maxfev=None, method='Nelder-Mead'):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param induction: int: key for dicts of arrays
    :param maxfev: int
    :param method: str
    :return: dict
    """
    if maxfev is None:
        maxfev = 600

    if method == 'Nelder-Mead':
        result = optimize.minimize(error_function, x, method='Nelder-Mead',
                                   options={'fatol': 1e-3, 'xatol': 1e-3, 'disp': True, 'maxiter': maxfev,
                                            'maxfev': maxfev},
                                   args=(xmin, xmax, induction))
    elif method == 'Powell':
        result = optimize.minimize(error_function, x, method='Powell',
                          options={'ftol': 1e-3, 'xtol': 1e-3, 'disp': True, 'maxiter': maxfev,
                                   'maxfev': maxfev},
                          args=(xmin, xmax, induction))
    formatted_x = '['+', '.join(['%.3E' % xi for xi in result.x])+']'
    print 'Contoller: completed optimize_polish for all cells after %i iterations with Error: %.4E and x: %s' % \
          (result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


def optimize_explore(x, xmin, xmax, error_function, induction=None, maxfev=None):
    """

    :param x: array
    :param xmin: array
    :param xmax: array
    :param error_function: callable
    :param induction: int: key for dicts of arrays
    :param maxfev: int
    :return: dict
    """
    if maxfev is None:
        maxfev = 700

    take_step = Normalized_Step(x, xmin, xmax)
    minimizer_kwargs = dict(method=null_minimizer, args=(xmin, xmax, induction))
    result = optimize.basinhopping(error_function, x, niter=maxfev, niter_success=maxfev/2,
                                   disp=True, interval=min(20, int(maxfev/20)), minimizer_kwargs=minimizer_kwargs,
                                   take_step=take_step)
    formatted_x = '['+', '.join(['%.3E' % xi for xi in result.x])+']'
    print 'Contoller: completed optimize_explore for all cells after %i iterations with Error: %.4E and x: %s' % \
          (result.nit, result.fun, formatted_x)
    return {'x': result.x, 'Err': result.fun}


x0 = [6.318E+02, 1.707E+02, 3.604E-01, 2.497E-01, 1.020E-01, 2.]  # Error: 4.274E+04 (max weight = 3)


xmin1 = [100., 100., 0., 0., 0., 1.5]
xmax1 = [6000., 6000., 1., 1., 1., 2.5]

induction = 1

result = optimize_explore(x0, xmin1, xmax1, ramp_population_error, induction, maxfev=1000)
x1 = result['x']
"""
polished_result = optimize_polish(x1, xmin1, xmax1, ramp_population_error, induction, maxfev=600, method='Nelder-Mead')  
    # method='Powell')
x1 = polished_result['x']
"""