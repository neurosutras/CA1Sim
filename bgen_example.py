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
    """

    :param parameters: array
    :return: dict
    """
    features = {}
    num_params = len(parameters)
    f1 = parameters[0]
    features['f1'] = f1
    g = 1. + 9. / (num_params - 1.) * np.sum(parameters[1:])
    features['g'] = g
    h = 1. - np.sqrt(f1 / g)
    features['h'] = h
    return features


def get_objectives(features):
    """

    :param features: dict
    :return: dict
    """
    objectives = {}
    for feature_name in ['f1', 'g', 'h']:
        objective_name = feature_name
        objectives[objective_name] = features[feature_name]
    f2 = features['g'] * features['h']
    objectives['f2'] = f2
    return objectives


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
    feature_names = ['g', 'h']
    objective_names = ['f1', 'f2']
    get_features = complex_problem
    get_objectives = get_objectives


bgen = BGen(x0, param_names, feature_names, objective_names, 200, bounds=bounds, evaluate=evaluate_basinhopping,
                 seed=0, max_gens=30, adaptive_step_interval=1, adaptive_step_factor=0.9, ngen_success=None,
                 survival_rate=0.10, disp=True)

for param_list in bgen():
    features = map(get_features, param_list)
    objectives = map(get_objectives, features)
    bgen.update_population(features, objectives)

bgen.storage.plot_history()