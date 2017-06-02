import numpy as np
from moopgen import *
import matplotlib as mpl

problem = 'complex'


def simple_problem(parameters):
    x, y = parameters

    f1 = 0.125 * (20 + x + y) ** 2
    f2 = 100 + x * y

    return {'f1': f1, 'f2': f2}


def complex_problem(parameters):
    num_params = len(parameters)
    f1 = parameters[0]
    result = {'f1': f1}
    g = 1. + 9. / (num_params - 1.) * np.sum(parameters[1:])
    h = 1. - np.sqrt(f1 / g)
    f2 = g * h
    result['f2'] = f2
    return result


if problem == 'simple':

    # Define the variables and their lower and upper search bounds
    param_names = ['x', 'y']
    bounds = [(-10., 10.), (-10., 10.)]
    x0 = [0., 0.]

    # Define the list of objectives
    objective_names = ['f1', 'f2']

    get_objectives = simple_problem

elif problem == 'complex':

    num_params = 5
    param_names = ['x%i' % i for i in range(num_params)]
    bounds = [(0., 1.) for i in range(num_params)]
    x0 = [0.5 * (xmin + xmax) for (xmin, xmax) in bounds]
    objective_names = ['f1', 'f2']
    get_objectives = complex_problem


bgen = BGen(x0, param_names, objective_names, 50, bounds=bounds, evaluate=evaluate_basinhopping,
                 seed=0, max_gens=100, adaptive_step_interval=3, adaptive_step_factor=0.9, ngen_success=None,
                 survival_rate=0.25, disp=True)

for param_list in bgen():
    result = map(get_objectives, param_list)
    bgen.set_objectives(result)

bgen.storage.plot_history()