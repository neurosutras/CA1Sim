__author__ = 'milsteina'
from specify_cells import *
from function_lib import *
import numpy as np
import scipy.optimize as optimize
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import scipy.stats as stats

# rec_filename = '032415 kap_kad_ih_scale kd pas no_na - EB1 - AMPAR_scaling'
# rec_filename = '032615 kap_kad_ih_scale kd pas no_na - EB2 - AMPAR_scaling'
# rec_filename = '043015 pas_exp_scale kdr ka_scale ih_sig_scale - EB2 - AMPAR_scaling'
# rec_filename = '072515 optimized basal ka_scale dend_sh_ar_nas - EB2 - AMPAR_scaling'
# rec_filename = '102915 interim dendritic excitability - AMPAR_scaling'
# rec_filename = '012816 altered intrinsic properties - AMPAR_scaling'
rec_filename = '020516 altered km2 rinp - AMPAR_scaling'

def fit_exp_linear(t, y, y0=0):
    """
    Removes the specified offset and fits the exponential as a line in log space.
    :param t:
    :param y:
    :param A0:
    :return:
    """
    y = np.log(y[1:] - y0)
    tau, logA = np.polyfit(t[1:], y, 1)
    return np.exp(logA), 1/tau


def exp_offset(x, y0, slope, tau):
    """

    :param x:
    :param y0:
    :param slope:
    :return:
    """
    return y0 + slope * (np.exp(x/tau) - 1.)

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

            y0 = sorted_dataset[0]
            A, tau = fit_exp_linear(sorted_distances, sorted_dataset, y0)
            y1 = y0-A
            A, tau = fit_exp_linear(sorted_distances, sorted_dataset, y1)
            fit = (y0 - A) + A * np.exp(interp_distances/tau)

            plt.scatter(sorted_distances, sorted_dataset, label=syn_type+': '+param_name)
            plt.plot(interp_distances, fit, label='fit')
            plt.xlabel('Distance from Soma (um)')
            plt.ylabel('Peak Conductance (uS)')
            plt.title('Mechanism Parameter Distribution')
            plt.legend(loc="best", scatterpoints=1, frameon=False, framealpha=0.5)
            plt.show()
        else:
            raise Exception('rec_filename is not formatted correctly or does not contain the specified data.')
    return [y0, A, tau]


def fit_synaptic_parameter_distribution2(rec_filename, sec_type, syn_type, param_name):
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

            popt, pcov = optimize.curve_fit(exp_offset, sorted_distances, sorted_dataset, p0=[7.e-4, 1.e-4, 89.])

            y0, A, tau = popt
            fit = (y0 - A) + A * np.exp(interp_distances/tau)

            plt.scatter(sorted_distances, sorted_dataset, label=syn_type+': '+param_name)
            plt.plot(interp_distances, fit, label='fit')
            plt.xlabel('Distance from Soma (um)')
            plt.ylabel('Peak Conductance (uS)')
            plt.title('Mechanism Parameter Distribution')
            plt.legend(loc="best", scatterpoints=1, frameon=False, framealpha=0.5)
            plt.show()
        else:
            raise Exception('rec_filename is not formatted correctly or does not contain the specified data.')
    return [y0, A, tau]


"""
# trunk parameters from both morphologies:
trunk_EB2 = [3.37e-5, 75.1, 2.30e-4]  # from '043015 pas_exp_scale kdr ka_scale ih_sig_scale - EB2' for 'AMPA_coop'
trunk_EB2 = [0.000101, 75.2, 0.000693]  # from '043015 pas_exp_scale kdr ka_scale ih_sig_scale - EB2' for 'AMPA_KIN'
trunk_av = np.mean([trunk_EB1, trunk_EB2], axis=0)
# trunk_av = [  3.55809459e-04,   4.51044089e+01,   8.74605934e-04]

basal_EB1 = [0.00040857888906237249, 52.796040616290718, 0.0007019697520230507]
basal_EB2 = [0.0007769851617510558, 70.413706295645554, 0.00053225698367113418]
basal_av = np.mean([basal_EB1, basal_EB2], axis=0)
# basal_av = [  5.92782025e-04,   6.16048735e+01,   6.17113368e-04]
x_av = basal_av
"""

result = [y0, A, tau] = fit_synaptic_parameter_distribution2(rec_filename, 'trunk', 'AMPA_KIN', 'gmax')
#x = fit_synaptic_parameter_distribution(rec_filename, 'basal', 'AMPA_KIN', 'gmax')
print result
