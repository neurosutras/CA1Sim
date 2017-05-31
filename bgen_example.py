import numpy as np
from moopgen import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib as mpl

# Define the variables and their lower and upper search bounds
param_names = ['x', 'y']
bounds = [(-10., 10.), (-10., 10.)]

# Define the list of objectives
objective_names = ['f1', 'f2']


# This is the function which is going to be minimized
def func_to_optimize(parameters):
    x, y = parameters

    f1 = 0.125 * (20 + x + y) ** 2
    f2 = 100 + x * y

    return {'f1': f1, 'f2': f2}


def plot_history(storage):
    for this_attr in ['fitness', 'energy', 'distance', 'survivor']:
        plt.figure()
        color = iter(cm.rainbow(np.linspace(0, 1, len(storage.history))))
        for population in storage.history:
            plt.scatter([indiv.rank for indiv in population], [getattr(indiv, this_attr) for indiv in population],
                        c=next(color))
        plt.title(this_attr)
    for i in range(len(storage.param_names)):
        this_attr = 'x'
        plt.figure()
        color = iter(cm.rainbow(np.linspace(0, 1, len(storage.history))))
        for population in storage.history:
            plt.scatter([indiv.rank for indiv in population], [getattr(indiv, this_attr)[i] for indiv in population],
                        c=next(color))
        plt.title(this_attr+str(i))
    for i in range(len(storage.param_names)):
        this_attr = 'objectives'
        plt.figure()
        color = iter(cm.rainbow(np.linspace(0, 1, len(storage.history))))
        for population in storage.history:
            plt.scatter([indiv.rank for indiv in population], [getattr(indiv, this_attr)[i] for indiv in population],
                        c=next(color))
        plt.title(this_attr+str(i))
    plt.show()
    plt.close()

x0 = [0., 0.]

bgen = BGen(x0, param_names, objective_names, 50, bounds=bounds, evaluate=evaluate_basinhopping,
                 seed=0, max_gens=100, adaptive_step_interval=10, adaptive_step_factor=0.9, ngen_success=None,
                 survival_rate=0.5, disp=True)

for param_list in bgen():
    result = map(func_to_optimize, param_list)
    bgen.set_objectives(result)

plot_history(bgen.storage)