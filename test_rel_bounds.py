from moopgen import *
"""
param_names = ['dend.gbar_nas min', 'dend.gbar_nas', 'soma.gbar_nas', 'axon.gbar_nax', 'ais.gbar_nax',
               'axon.gkabar', 'soma.gkabar', 'dend.gkabar', 'logrange8', 'logrange9', 'log10', 'log11']

Rules:
DEP                 IND
soma.gbar_nas     > dend.gbar_nas min
dend.gbar_nas     < soma.gbar_nas
dend.gbar_nas min < dend.gbar_nas
axon.gbar_nax     > soma.gbar_nas
ais.gbar_nax      > 2. * axon.gbar_nax
dend.gkabar       > soma.gkabar
axon.gkabar       < 3. * soma.gkabar

rel_bounds = [['soma.gbar_nas', ">", 1., 'dend.gbar_nas min'], ['dend.gbar_nas', "<", 1., 'soma.gbar_nas'],
              ['dend.gbar_nas', ">", 1., 'dend.gbar_nas min'], ['axon.gbar_nax', ">", 1., 'soma.gbar_nas'],
              ['ais.gbar_nax', ">", 2., 'axon.gbar_nax'], ['dend.gkabar', ">", 1., 'soma.gkabar'],
              ['axon.gkabar', "<", 3., 'soma.gkabar'], ['logrange8', "<=", 1., 'logrange9'], ['log10', ">=", 2., 'log11']]
bounds = [(0., 0.), (0.01, 0.05), (0.01, 0.05), (0.02, 0.1), (0.02, 0.5), (0.01, 0.18), (0.01, 0.05), (0.01, 0.25),
          (0.01, 100.), (0.01, 100.), (0.01, 1000.), (0.001, 500.)]
x0 = np.array([0., 0.03, 0.03, 0.06, 0.1681, 0.05266, 0.02108, 0.04, 0.05, 0.05, 0.2, 0.12])
"""

param_names = ['log10', 'log11']
rel_bounds = [['log10', ">=", 2., 'log11']]
bounds = [(0.01, 1.), (0.001, 0.5)]
x0 = np.array([0.2, 0.12])


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



step = RelativeBoundedStep(x0, param_names, bounds, rel_bounds, stepsize=0.001, wrap=False)
prev_x = x0
x_history = []
abs_bounds_failed = 0
rel_bounds_failed = 0
for i in range(20000):
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


step2 = RelativeBoundedStep(x0, param_names, bounds, rel_bounds, stepsize=0.001, wrap=True)
prev_x = x0
x_history_wrap = []
failed = 0
for i in range(20000):
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
