import numpy as np
from copy import deepcopy
from plot_results import *


class RelativeBoundedStep(object):
    """
    Step-taking method for use with BGen. Steps each parameter within specified bounds. Explores the range in log10
    space when the range is >= 2 orders of magnitude. Uses the log-modulus transformation (John & Draper, 1980) as an
    approximation that tolerates ranges that span zero. If bounds are not provided for some parameters, the default is
    (0.1 * x0, 10. * x0).
    """
    def __init__(self, x0, param_names, bounds=None, rel_bounds=None, stepsize=0.5, wrap=False, random=None, disp=False,
                 **kwargs):
        """

        :param x0: array
        :param param_names: list
        :param bounds: list of tuple
        :param rel_bounds: list of lists
        :param stepsize: float in [0., 1.]
        :param wrap: bool  # whether or not to wrap around bounds
        :param random: int or :class:'np.random.RandomState'
        :param disp: bool
        """
        self.disp = disp
        self.wrap = wrap
        self.stepsize = stepsize
        if x0 is None and bounds is None:
            raise ValueError('RelativeBoundedStep: Either starting parameters or bounds are missing.')
        if random is None:
            self.random = np.random
        else:
            self.random = random
        if bounds is None:
            xmin = [None for xi in x0]
            xmax = [None for xi in x0]
        else:
            xmin = [bound[0] for bound in bounds]
            xmax = [bound[1] for bound in bounds]
        if x0 is None:
            x0 = [None for i in xrange(len(bounds))]
        for i in xrange(len(x0)):
            if x0[i] is None:
                if xmin[i] is None or xmax[i] is None:
                    raise ValueError('RelativeBoundedStep: Either starting parameters or bounds are missing.')
                else:
                    x0[i] = 0.5 * (xmin[i] + xmax[i])
            if xmin[i] is None:
                if x0[i] > 0.:
                    xmin[i] = 0.1 * x0[i]
                elif x0[i] == 0.:
                    xmin[i] = -1.
                else:
                    xmin[i] = 10. * x0[i]
            if xmax[i] is None:
                if x0[i] > 0.:
                    xmax[i] = 10. * x0[i]
                elif x0[i] == 0.:
                    xmax[i] = 1.
                else:
                    xmax[i] = 0.1 * x0[i]
        self.x0 = np.array(x0)
        self.param_names = param_names
        self.param_indexes = {param: i for i, param in enumerate(param_names)}
        self.xmin = np.array(xmin)
        self.xmax = np.array(xmax)
        self.x_range = np.subtract(self.xmax, self.xmin)
        self.logmod = lambda x, offset, factor: np.log10(x * factor + offset)
        self.logmod_inv = lambda logmod_x, offset, factor: ((10. ** logmod_x) - offset) / factor
        self.abs_order_mag = []
        for i in xrange(len(xmin)):
            xi_logmin, xi_logmax, offset, factor = self.logmod_bounds(xmin[i], xmax[i])
            self.abs_order_mag.append(xi_logmax - xi_logmin)
        self.rel_bounds = rel_bounds

    def __call__(self, current_x, stepsize=None, wrap=None):
        """
        Take a step within bounds. If stepsize or wrap is specified for an individual call, it overrides the default.
        :param current_x: array
        :param stepsize: float in [0., 1.]
        :param wrap: bool
        :return: array
        """
        if stepsize is None:
            stepsize = self.stepsize
        if wrap is None:
            wrap = self.wrap
        x = np.array(current_x)
        for i in xrange(len(x)):
            if not self.xmax[i] >= self.xmin[i]:
                raise Exception('Bounds for paramter %d: max is not >= to min.') %i
            new_xi = self.generate_param(x[i], i, self.xmin[i], self.xmax[i], stepsize, wrap, self.disp)
            x[i] = new_xi
        if self.rel_bounds is not None:
            x = self.apply_rel_bounds(x, stepsize, wrap, self.rel_bounds)
        return x

    def logmod_bounds(self, xi_min, xi_max):
        """

        :param xi_min: float
        :param xi_max: float
        :return: xi_logmin, xi_logmax, offset, factor
        """
        if xi_min < 0.:
            if xi_max < 0.:
                offset = 0.
                factor = -1.
            elif xi_max == 0.:
                offset = 0.1
                factor = -1.
            else:
                #If xi_min and xi_max are opposite signs, do not sample in log space; do linear sampling
                return 0., 0., None, None
            xi_logmin = self.logmod(xi_max, offset, factor)  # When the sign is flipped, the max and min will reverse
            xi_logmax = self.logmod(xi_min, offset, factor)
        elif xi_min == 0.:
            if xi_max == 0.:
                return 0., 0., None, None
            else:
                offset = 0.1
                factor = 1.
                xi_logmin = self.logmod(xi_min, offset, factor)
                xi_logmax = self.logmod(xi_max, offset, factor)
        else:
            offset = 0.
            factor = 1.
            xi_logmin = self.logmod(xi_min, offset, factor)
            xi_logmax = self.logmod(xi_max, offset, factor)
        return xi_logmin, xi_logmax, offset, factor

    def logmod_inv_bounds(self, xi_logmin, xi_logmax, offset, factor):
        """

        :param xi_logmin: float
        :param xi_logmax: float
        :param offset: float
        :param factor: float
        :return: xi_min, xi_max
        """
        if factor < 0.:
            xi_min = self.logmod_inv(xi_logmax, offset, factor)
            xi_max = self.logmod_inv(xi_logmin, offset, factor)
        else:
            xi_min = self.logmod_inv(xi_logmin, offset, factor)
            xi_max = self.logmod_inv(xi_logmax, offset, factor)
        return xi_min, xi_max

    def generate_param(self, xi, i, xi_min, xi_max, stepsize, wrap, disp=False):
        """

        :param xi: float
        :param i: int
        :param min: float
        :param max: float
        :param stepsize: float
        :param wrap: bool
        :return:
        """
        if xi_min == xi_max:
            return xi_min
        xi_logmin, xi_logmax, offset, factor = self.logmod_bounds(xi_min, xi_max)
        if offset is None and factor is None:
            new_xi = self.linear_step(xi, i, xi_min, xi_max, stepsize, wrap, disp)
        else:
            order_mag = min(xi_logmax - xi_logmin, self.abs_order_mag[i] * stepsize)
            if order_mag <= 1.:
                new_xi = self.linear_step(xi, i, xi_min, xi_max, stepsize, wrap, disp)
            else:
                new_xi = self.log10_step(xi, i, xi_logmin, xi_logmax, offset, factor, stepsize, wrap, disp)
        return new_xi

    def linear_step(self, xi, i, xi_min, xi_max, stepsize=None, wrap=None, disp=False):
        """
        Steps the specified parameter within the bounds according to the current stepsize.
        :param xi: float
        :param i: int
        :param stepsize: float in [0., 1.]
        :param wrap: bool
        :return: float
        """
        if stepsize is None:
            stepsize = self.stepsize
        if wrap is None:
            wrap = self.wrap
        step = stepsize * self.x_range[i] / 2.
        if disp:
            print 'Before: xi: %.4f, step: %.4f, xi_min: %.4f, xi_max: %.4f' % (xi, step, xi_min, xi_max)
        if wrap:
            step = min(step, xi_max - xi_min)
            delta = self.random.uniform(-step, step)
            new_xi = xi + delta
            if xi_min > new_xi:
                new_xi = max(xi_max - (xi_min - new_xi), xi_min)
            elif xi_max < new_xi:
                new_xi = min(xi_min + (new_xi - xi_max), xi_max)
        else:
            xi_min = max(xi_min, xi - step)
            xi_max = min(xi_max, xi + step)
            new_xi = self.random.uniform(xi_min, xi_max)
        if disp:
            print 'After: xi: %.4f, step: %.4f, xi_min: %.4f, xi_max: %.4f' % (new_xi, step, xi_min, xi_max)
        linear_steps['steps'][i] += 1
        return new_xi

    def log10_step(self, xi, i, xi_logmin, xi_logmax, offset, factor, stepsize=None, wrap=None, disp=False):
        """
        Steps the specified parameter within the bounds according to the current stepsize.
        :param xi: float
        :param i: int
        :param xi_logmin: float
        :param xi_logmax: float
        :param offset: float.
        :param factor: float
        :param stepsize: float in [0., 1.]
        :param wrap: bool
        :return: float
        """
        if stepsize is None:
            stepsize = self.stepsize
        if wrap is None:
            wrap = self.wrap
        xi_log = self.logmod(xi, offset, factor)
        step = stepsize * self.abs_order_mag[i] / 2.
        if disp:
            print 'Before: log_xi: %.4f, step: %.4f, xi_logmin: %.4f, xi_logmax: %.4f' % (xi_log, step, xi_logmin,
                                                                                          xi_logmax)
        if wrap:
            step = min(step, xi_logmax - xi_logmin)
            delta = np.random.uniform(-step, step)
            step_xi_log = xi_log + delta
            if xi_logmin > step_xi_log:
                step_xi_log = max(xi_logmax - (xi_logmin - step_xi_log), xi_logmin)
            elif xi_logmax < step_xi_log:
                step_xi_log = min(xi_logmin + (step_xi_log - xi_logmax), xi_logmax)
            new_xi = self.logmod_inv(step_xi_log, offset, factor)
            log_steps['steps'][i] += 1
        else:
            step_xi_logmin = max(xi_logmin, xi_log - step)
            step_xi_logmax = min(xi_logmax, xi_log + step)
            new_xi_log = self.random.uniform(step_xi_logmin, step_xi_logmax)
            new_xi = self.logmod_inv(new_xi_log, offset, factor)
            log_steps['steps'][i] += 1
        if disp:
            print 'After: xi: %.4f, step: %.4f, xi_logmin: %.4f, xi_logmax: %.4f' % (new_xi, step, xi_logmin,
                                                                                      xi_logmax)
        return new_xi

    def apply_rel_bounds(self, x, stepsize, wrap, rel_bounds=None, disp=False):
        """

        :param x: array
        :param rel_bounds: list of lists
        :return:
        """
        if disp:
            print 'orig x: %s' % str(x)
        new_x = np.array(x)
        new_min = deepcopy(self.xmin)
        new_max = deepcopy(self.xmax)
        if rel_bounds is not None:
            for i, rel_bound_rule in enumerate(rel_bounds):
                dep_param = rel_bound_rule[0]  #Dependent param: name of the parameter that may be modified
                dep_param_ind = self.param_indexes[dep_param]
                if dep_param_ind >= len(x):
                    raise Exception('Dependent parameter index is out of bounds for rule %d.' %i)
                factor = rel_bound_rule[2]
                ind_param = rel_bound_rule[3]  #Independent param: name of the parameter that sets the bounds
                ind_param_ind = self.param_indexes[ind_param]
                if ind_param_ind >= len(x):
                    raise Exception('Independent parameter index is out of bounds for rule %d.' %i)
                if rel_bound_rule[1] == "=":
                    new_xi = factor * new_x[ind_param_ind]
                    if (new_xi >= self.xmin[dep_param_ind]) and (new_xi < self.xmax[dep_param_ind]):
                        new_x[dep_param_ind] = new_xi
                    else:
                        raise Exception('Relative bounds rule %d contradicts fixed parameter bounds.' %i)
                    continue
                if disp:
                    print 'Before rel bound rule %i. xi: %.4f, min: %.4f, max: %.4f' % (i, new_x[dep_param_ind],
                                                                                        new_min[dep_param_ind],
                                                                                        new_max[dep_param_ind])

                if rel_bound_rule[1] == "<":
                    rel_max = factor * new_x[ind_param_ind]
                    new_max[dep_param_ind] = max(min(new_max[dep_param_ind], rel_max), new_min[dep_param_ind])
                elif rel_bound_rule[1] == "<=":
                    rel_max = factor * new_x[ind_param_ind]
                    new_max[dep_param_ind] = max(min(new_max[dep_param_ind], np.nextafter(rel_max, rel_max + 1)),
                                                 new_min[dep_param_ind])
                elif rel_bound_rule[1] == ">=":
                    rel_min = factor * new_x[ind_param_ind]
                    new_min[dep_param_ind] = min(max(new_min[dep_param_ind], rel_min), new_max[dep_param_ind])
                elif rel_bound_rule[1] == ">":
                    rel_min = factor * new_x[ind_param_ind]
                    new_min[dep_param_ind] = min(max(new_min[dep_param_ind], np.nextafter(rel_min, rel_min + 1)),
                                                 new_max[dep_param_ind])
                if not (new_x[dep_param_ind] >= new_min[dep_param_ind] and new_x[dep_param_ind] < new_max[dep_param_ind]):
                    new_xi = max(new_x[dep_param_ind], new_min[dep_param_ind])
                    new_xi = min(new_xi, new_max[dep_param_ind])
                    if disp:
                        print 'After rel bound rule %i. xi: %.4f, min: %.4f, max: %.4f' % (i, new_xi,
                                                                                           new_min[dep_param_ind],
                                                                                           new_max[dep_param_ind])
                    new_x[dep_param_ind] = self.generate_param(new_xi, dep_param_ind, new_min[dep_param_ind],
                                                               new_max[dep_param_ind], stepsize, wrap, disp=disp)
        return new_x




param_names = ['dend.gbar_nas min', 'dend.gbar_nas', 'soma.gbar_nas', 'axon.gbar_nax', 'ais.gbar_nax',
               'axon.gkabar', 'soma.gkabar', 'dend.gkabar', 'logrange8', 'logrange9', 'log10', 'log11']
"""
Rules:
DEP                 IND
soma.gbar_nas     > dend.gbar_nas min
dend.gbar_nas     < soma.gbar_nas
dend.gbar_nas min < dend.gbar_nas
axon.gbar_nax     > soma.gbar_nas
ais.gbar_nax      > 2. * axon.gbar_nax
dend.gkabar       > soma.gkabar
axon.gkabar       < 3. * soma.gkabar
"""

rel_bounds = [['soma.gbar_nas', ">", 1., 'dend.gbar_nas min'], ['dend.gbar_nas', "<", 1., 'soma.gbar_nas'],
              ['dend.gbar_nas', ">", 1., 'dend.gbar_nas min'], ['axon.gbar_nax', ">", 1., 'soma.gbar_nas'],
              ['ais.gbar_nax', ">", 2., 'axon.gbar_nax'], ['dend.gkabar', ">", 1., 'soma.gkabar'],
              ['axon.gkabar', "<", 3., 'soma.gkabar'], ['logrange8', "<=", 1., 'logrange9'], ['log10', ">=", 2., 'log11']]

#rel_bounds = []
bounds = [(0., 0.), (0.01, 0.05), (0.01, 0.05), (0.02, 0.1), (0.02, 0.5), (0.01, 0.18), (0.01, 0.05), (0.01, 0.25),
          (0.01, 100.), (0.01, 100.), (0.01, 500.), (0.001, 100.)]
x0 = np.array([0., 0.03, 0.03, 0.06, 0.1681, 0.05266, 0.02108, 0.04, 0.05, 0.05, 0.2, 0.12])

"""
param_names = ['log0', 'log1']
rel_bounds = [['log0', "<=", -2., 'log1']]
#rel_bounds = []
bounds = [(-20000., 100.), (-10., 100.)]
#x0 = np.array([0.2, 0.12])
x0 = np.array([-100., 100.])
"""

linear_steps = {'order_mag': [], 'steps': [0 for key in param_names]}
log_steps = {'order_mag': [], 'steps': [0 for key in param_names]}


def check_abs_bounds(x, bounds):
    for i, xi in enumerate(x):
        if (xi < bounds[i][0]):
            print 'Paramater %d: value %.3f is less than minimum bound' % (i, xi)
            return False
        if (xi == bounds[i][0] and xi == bounds[i][1]):
            return True
        if (xi >= bounds[i][1]):
            print 'Paramater %d: value %.3f is greater than maximum bound' % (i, xi)
            return False
    return True


def check_rel_bounds(x, rel_bounds, param_indexes):
    for r, rule in enumerate(rel_bounds):
        dep_param_ind = param_indexes[rule[0]]  # Dependent param. index: index of the parameter that may be modified
        if dep_param_ind >= len(x):
            raise Exception('Dependent parameter index is out of bounds for rule %d.' % r)
        factor = rule[2]
        ind_param_ind = param_indexes[rule[3]]  # Independent param. index: index of the parameter that sets the bounds
        if ind_param_ind >= len(x):
            raise Exception('Independent parameter index is out of bounds for rule %d.' % r)
        if rule[1] == "=":
            operator = lambda x, y: x == y
        elif rule[1] == "<":
            operator = lambda x, y: x < y
        elif rule[1] == "<=":
            operator = lambda x, y: x <= y
        elif rule[1] == ">=":
            operator = lambda x, y: x >= y
        elif rule[1] == ">":
            operator = lambda x, y: x > y
        if not operator(x[dep_param_ind], factor * x[ind_param_ind]):
            print 'Parameter %d: value %.3f did not meet relative bound in rule %d.' % \
                  (dep_param_ind, x[dep_param_ind], r)
            return False
    return True

num = 20000
this_step_size = 0.1
param_history = lambda i: np.array([this_p[i] for this_p in x_history])


step = RelativeBoundedStep(x0, param_names, bounds, rel_bounds, stepsize=this_step_size, wrap=False)
prev_x = x0
x_history = []
abs_bounds_failed = 0
rel_bounds_failed = 0
for i in range(num):
    new_x = step(prev_x)
    if not check_abs_bounds(new_x, bounds):
        print 'Absolute bounds failed for x: %s' % str(new_x)
        abs_bounds_failed += 1
    if not check_rel_bounds(new_x, rel_bounds, step.param_indexes):
        print 'Relative bounds failed for x: %s' % str(new_x)
        rel_bounds_failed += 1
    prev_x = new_x
    x_history.append(prev_x)
print 'Relative bound test (no wrap) ended after %i iterations with %i abs failures and %i rel failures.' % \
      ((i + 1), abs_bounds_failed, rel_bounds_failed)


step2 = RelativeBoundedStep(x0, param_names, bounds, rel_bounds, stepsize=this_step_size, wrap=True)
prev_x = x0
x_history_wrap = []
failed = 0
for i in range(num):
    new_x = step2(prev_x)
    if not check_abs_bounds(new_x, bounds):
        print 'Absolute bounds failed for x: %s' % str(new_x)
        abs_bounds_failed += 1
    if not check_rel_bounds(new_x, rel_bounds, step.param_indexes):
        print 'Relative bounds failed for x: %s' % str(new_x)
        rel_bounds_failed += 1
    prev_x = new_x
    x_history_wrap.append(prev_x)
print 'Relative bound test (wrap) ended after %i iterations with %i abs failures and %i rel failures.' % \
      ((i + 1), abs_bounds_failed, rel_bounds_failed)


p0 = param_history(0)
p1 = param_history(1)

get_num_samples = lambda p, low, high: len(np.where((p >= low) & (p < high))[0])
