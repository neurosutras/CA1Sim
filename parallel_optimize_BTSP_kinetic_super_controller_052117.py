__author__ = 'milsteina'
from function_lib import *
import parallel_optimize_BTSP_kinetic_sub_controller_052117
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
    cluster_id = sys.argv[1]
    c = Client(cluster_id=cluster_id)
else:
    cluster_id = None
    c = Client()

if len(sys.argv) > 2:
    sub_controller_id = int(sys.argv[2])
else:
    sub_controller_id = 0

if len(sys.argv) > 3:
    group_size = int(sys.argv[3])
else:
    group_size = 1


global_start_time = time.time()
dv = c[::group_size]
dv.block = True

for id in dv.targets:
    c[id].execute('run parallel_optimize_BTSP_kinetic_sub_controller_052117 %s %i %i ' %
                  (cluster_id, id, group_size), block=True)


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


result = dv.apply_sync(parallel_optimize_BTSP_kinetic_sub_controller_052117.report, x0, induction)
print result

"""
result = optimize_explore(x0, xmin1, xmax1, ramp_population_error, induction, maxfev=1000)
x1 = result['x']

polished_result = optimize_polish(x1, xmin1, xmax1, ramp_population_error, induction, maxfev=600, method='Nelder-Mead')  
    # method='Powell')
x1 = polished_result['x']
"""
