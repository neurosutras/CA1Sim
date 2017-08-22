from moopgen import *

param_names = {0: 'dend.gbar_nas min', 1: 'dend.gbar_nas', 2: 'soma.gbar_nas', 3: 'axon.gbar_nax', 4: 'ais.gbar_nax',
               5: 'axon.gkabar', 6: 'soma.gkabar', 7: 'dend.gkabar', 8: 'soma.gkdrbar'}
"""
Rules:
DEP                 IND
dend.gbar_nas     < soma.gbar_nas
dend.gbar_nas min < dend.gbar_nas
axon.gbar_nax     > soma.gbar_nas
ais.gbar_nax      > 2. * axon.gbar_nax
dend.gkabar       > soma.gkabar
axon.gkabar       < 3. * soma.gkabar
"""
rel_bounds = [[1, "<", 1., 2], [0, "<", 1., 1], [3, ">", 1., 2], [4, ">", 2., 3], [7, ">", 1., 6], [5, "<", 3., 6]]
bounds = [(0., 0.015), (0.01, 0.05), (0.01, 0.05), (0.02, 0.1), (0.02, 0.5), (0.01, 0.18), (0.01, 0.05), (0.01, 0.25),
          (0.01, 0.06)]
x0 = np.array([0., 0.03, 0.03, 0.06, 0.1681, 0.05266, 0.02108, 0.04, 0.04299])


def check_bounds(x, bounds, rel_bonds):
    for i, xi in enumerate(x):
        if (xi < bounds[i][0]):
            print 'Paramater %d: value %.3f is less than minimum bound' % (i, xi)
        if (xi >= bounds[i][1]):
            print 'Paramater %d: value %.3f is greater than maximum bound' % (i, xi)
    for r, rule in enumerate(rel_bounds):
        dep_param_ind = rule[0]  # Dependent param. index: index of the parameter that may be modified
        if dep_param_ind >= len(x):
            raise Exception('Dependent parameter index is out of bounds for rule %d.' % r)
        factor = rule[2]
        ind_param_ind = rule[3]  # Independent param. index: index of the parameter that sets the bounds
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
        if operator(x[dep_param_ind], factor * x[ind_param_ind]) is False:
            print 'Parameter %d: value %.3f did not meet relative bound in rule %d.' %(dep_param_ind, x[dep_param_ind], r)


step = RelativeBoundedStep(x0, bounds, rel_bounds, wrap=False)
prev_x = x0
for i in range(50):
    new_x = step(prev_x)
    if check_bounds(new_x, bounds, rel_bounds) is False:
        break
    prev_x = new_x
print 'Completed no wrap test.'

step2 = RelativeBoundedStep(x0, bounds, rel_bounds, wrap=True)
prev_x = x0
for i in range(50):
    new_x = step2(prev_x)
    if check_bounds(new_x, bounds, rel_bounds) is False:
        break
    prev_x = new_x
print 'Completed wrap test.'

