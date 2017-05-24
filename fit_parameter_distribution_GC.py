__author__ = 'milsteina' and 'Grace Ng'
from specify_cells2 import *
from function_lib import *
import numpy as np
import scipy.optimize as optimize
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import scipy.stats as stats

rec_filename = '042617 GC optimizing spike stability - AMPAR_scaling'


"""
def exp_offset(x, y0, x0, slope, tau):


    :param x:
    :param y0:
    :param x0:
    :param slope:
    :return:

    return y0 + slope * (np.exp((x-x0)/tau) - 1.)
"""

def exp_offset(x, y0, slope, tau):
    """

    :param x:
    :param y0:
    :param x0:
    :param slope:
    :return:
    """
    return y0 + slope * (np.exp((x)/tau) - 1.)

def select_synpases_for_fitting(rec_filename, syn_type, param_name, branch_cutoff):
    """
    Adds synapses that are non-terminal and have a branch order less than the cut-off to a list; the synapses in this
    list will be fitted to a distribution.
    :param rec_filename:
    :param syn_type:
    :param param_name:
    :return:
    """
    node_list = []
    with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
        if 'syn_type' in f.itervalues().next().attrs.keys():
            if f.itervalues().next().attrs['syn_type'] != syn_type:
                raise Exception('Synapses in this hdf5 file are not of the specified synapse type.')
        if param_name not in f.itervalues().next().attrs.keys():
            raise Exception('This parameter is not saved in this hdf5 file.')
        dataset = []
        distances = []
        for sim in f.itervalues():
            input_loc = sim.attrs['input_loc']
            is_terminal = sim['rec']['2'].attrs['is_terminal']
            branch_order = sim['rec']['2'].attrs['branch_order']
            if branch_order <= branch_cutoff and not is_terminal:
            #found that all branches with order >=5 are terminal
            #if not is_terminal:
            #if branch_order <= branch_cutoff or not is_terminal:
            # more error
                node_list.append(sim['rec']['2'].attrs['index'])
                distances.append(sim['rec']['2'].attrs['soma_distance'])
                dataset.append(sim.attrs[param_name])
        indexes = range(len(distances))
        indexes.sort(key=distances.__getitem__)
        sorted_distances = np.array(map(distances.__getitem__, indexes))
        sorted_dataset = np.array(map(dataset.__getitem__, indexes))
    print set(node_list)
    return sorted_dataset, sorted_distances, syn_type, param_name


def fit_synaptic_parameter_distribution(sorted_dataset, sorted_distances, syn_type, param_name):
    """
    Expects file to contain dendritic locations and values of synaptic point_processs parameters. Fits the specified
    param_name along the specified sec_type to an exponential distribution.
    :param rec_filename: str
    :param sec_type: str
    :param param_name: str
    """
    interp_distances = np.arange(0, sorted_distances[-1], 1.)
    popt, pcov = optimize.curve_fit(exp_offset, sorted_distances, sorted_dataset, p0=[7.e-4, 0., 89.])

    y0, A, tau = popt
    fit = (y0 - A) + A * np.exp((interp_distances)/tau)

    plt.scatter(sorted_distances, sorted_dataset, label=syn_type+': '+param_name)
    plt.plot(interp_distances, fit, label='fit')
    plt.xlabel('Distance from Soma (um)')
    plt.ylabel('Peak Conductance (uS)')
    plt.title('Mechanism Parameter Distribution')
    plt.legend(loc="best", scatterpoints=1, frameon=False, framealpha=0.5)
    plt.show()
    print popt, pcov
    return [y0, A, tau]


# DG Granule Cells
sorted_dataset, sorted_distances, syn_type, param_name = select_synpases_for_fitting(rec_filename, 'AMPA_KIN', 'gmax', 5)
result = [y0, A, tau] = fit_synaptic_parameter_distribution(sorted_dataset, sorted_distances, syn_type, param_name)
print result

#[y0, x0, A, tau] = [0.0017237199429116546, 14.361692003965528, 7.977391921662206e-05, 79.081981369526119]
#[y0, A, tau] = [  1.71052379e-03   6.65158925e-05   7.90793466e+01]  # 051917
# min_distance = 140.625, y(min_distance) = 0.0020378