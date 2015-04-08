__author__ = 'milsteina'
from specify_cells import *
from function_lib import *
import numpy as np
import scipy.optimize as optimize
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
"""

"""
#rec_filename = '032415 kap_kad_ih_scale kd pas no_na - EB1 - AMPAR_scaling'
rec_filename = '032615 kap_kad_ih_scale kd pas no_na - EB2 - AMPAR_scaling'


def fit_exp_linear(t, y, A0=0):
    """
    Removes the specified offset and fits the exponential as a line in log space.
    :param t:
    :param y:
    :param A0:
    :return:
    """
    y = np.log(y[1:]-A0)
    tau, logA = np.polyfit(t[1:], y, 1)
    return np.exp(logA), 1/tau


def fit_synaptic_parameter_distribution(rec_filename, sec_type, syn_type, param_name):
    """
    Expects file to contain dendritic locations and values of synaptic point_processs parameters. Fits the specified
    param_name along the specified sec_type to an exponential distribution.
    :param rec_filename: str
    :param sec_type: str
    :param param_name: str
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
        if f.attrs['syn_type'] == syn_type and sec_type in f and param_name in f[sec_type]:
            dataset = f[sec_type][param_name]
            distances = f[sec_type]['distances']
            indexes = range(len(distances))
            indexes.sort(key=distances.__getitem__)
            sorted_distances = np.array(map(distances.__getitem__, indexes))
            sorted_dataset = np.array(map(dataset.__getitem__, indexes))
            interp_distances = np.arange(0, sorted_distances[-1], 1.)
            A0 = sorted_dataset[0]
            A, tau = fit_exp_linear(sorted_distances, sorted_dataset, A0)
            A1 = min(A0, sorted_dataset[0]-A)
            A, tau = fit_exp_linear(sorted_distances, sorted_dataset, A1)
            fit = model_scaled_exp(interp_distances, A, tau, A1)
            plt.plot(sorted_distances, sorted_dataset, label=syn_type+':'+param_name)
            plt.plot(interp_distances, fit, label='fit')
            gfit = model_scaled_exp(interp_distances, x_av[0], x_av[1], x_av[2])
            plt.plot(interp_distances, gfit, label='global fit')
            plt.legend(loc="best")
            plt.show()
        else:
            raise Exception('rec_filename is not formatted correctly or does not contain the specified data.')
    return [A, tau, A1]
    #return sorted_distances, sorted_dataset, interp_distances

x_EB1 = [0.00039168828539039226, 46.765092844753902, 0.00089820146712824803]
x_EB2 = [0.00031993063270730745, 43.44372490545944, 0.00085101040031032116]
x_av = [  3.55809459e-04,   4.51044089e+01,   8.74605934e-04]
#x, y, fit = fit_synaptic_parameter_distribution(rec_filename, 'trunk', 'AMPA_KIN', 'gmax')
x = fit_synaptic_parameter_distribution(rec_filename, 'trunk', 'AMPA_KIN', 'gmax')
print x