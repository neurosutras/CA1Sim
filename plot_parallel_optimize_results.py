__author__ = 'milsteina'
from function_lib import *
from matplotlib import cm
import matplotlib.lines as mlines
import matplotlib as mpl
import numpy as np
import scipy.signal as signal
import scipy.stats as stats

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.size'] = 16.
# mpl.rcParams['font.size'] = 14.
#mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.sans-serif'] = 'Calibri'
# mpl.rcParams['font.sans-serif'] = 'Myriad Pro'
mpl.rcParams['text.usetex'] = False
#mpl.rcParams['figure.figsize'] = 6, 4.3
"""
mpl.rcParams['axes.labelsize'] = 'larger'
mpl.rcParams['axes.titlesize'] = 'xx-large'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'
mpl.rcParams['legend.fontsize'] = 'x-large'
"""


def plot_exported_BTSP_model_ramp_features3(processed_export_file_path):
    """

    :param processed_export_file_path: str (path)
    """
    orig_fontsize = mpl.rcParams['font.size']
    # mpl.rcParams['font.size'] = 20.
    with h5py.File(processed_export_file_path, 'r') as f:
        description = 'model_ramp_context'
        group = f[description]
        peak_locs = group['peak_locs'][:]
        binned_x = group['binned_x'][:]
        local_signal_range = group['local_signal_range'][:]
        norm_depot_rate = group['norm_depot_rate'][:]
        peak_weight = group.attrs['peak_weight']
        fig, axes = plt.subplots(1)
        axes.plot(local_signal_range, norm_depot_rate)
        axes.set_xlabel('Local signal amplitude (a.u.)')
        axes.set_ylabel('Normalized rate')
        axes.set_title('Depotentiation rate')
        clean_axes(axes)
        fig.tight_layout()
        plt.show()
        plt.close()
        description = 'model_ramp_features'
        for cell_key in f[description]:
            fig, axes = plt.subplots(2, 2)
            ymin = -1.
            ymax = 10.
            for induction_key in f[description][cell_key]:
                i = int(induction_key) - 1
                group = f[description][cell_key][induction_key]
                weights = group['weights'][:]
                initial_weights = group['initial_weights'][:]
                delta_weights = np.subtract(weights, initial_weights)
                target_ramp = group['target_ramp'][:]
                model_ramp = group['model_ramp'][:]
                axes[i][1].plot(peak_locs, delta_weights)
                axes[i][0].plot(binned_x, target_ramp, label='Experiment')
                axes[i][0].plot(binned_x, model_ramp, label='Model')
                axes[i][1].set_ylabel('Change in synaptic weight (a.u.)')
                axes[i][1].set_xlabel('Location (cm)')
                axes[i][0].set_ylabel('Ramp amplitude (mV)')
                axes[i][0].set_xlabel('Location (cm)')
                axes[i][1].set_title('Induction: %s' % induction_key, fontsize=mpl.rcParams['font.size'])
                axes[i][0].legend(loc='best', frameon=False, framealpha=0.5)
                ymin = min(ymin, np.min(model_ramp) - 1., np.min(target_ramp) - 1.)
                ymax = max(ymax, np.max(model_ramp) + 1., np.max(target_ramp) + 1.)
            for induction_key in f[description][cell_key]:
                i = int(induction_key) - 1
                axes[i][0].set_ylim([ymin, ymax])
                axes[i][1].set_ylim([-peak_weight, peak_weight])
            clean_axes(axes)
            fig.suptitle('Cell_id: %s' % cell_key, fontsize=mpl.rcParams['font.size'])
            fig.tight_layout()
            plt.show()
            plt.close()
    mpl.rcParams['font.size'] = orig_fontsize


def plot_exported_BTSP_model_ramp_features4(processed_export_file_path):
    """

    :param processed_export_file_path: str (path)
    """
    orig_fontsize = mpl.rcParams['font.size']
    # mpl.rcParams['font.size'] = 20.
    with h5py.File(processed_export_file_path, 'r') as f:
        description = 'model_ramp_context'
        group = f[description]
        peak_locs = group['peak_locs'][:]
        binned_x = group['binned_x'][:]
        peak_weight = group.attrs['peak_weight']
        description = 'model_ramp_features'
        for cell_key in f[description]:
            fig, axes = plt.subplots(2, 2)
            ymin = -1.
            ymax = 10.
            for induction_key in f[description][cell_key]:
                i = int(induction_key) - 1
                group = f[description][cell_key][induction_key]
                weights = group['weights'][:]
                initial_weights = group['initial_weights'][:]
                delta_weights = np.subtract(weights, initial_weights)
                target_ramp = group['target_ramp'][:]
                model_ramp = group['model_ramp'][:]
                axes[i][1].plot(peak_locs, delta_weights)
                axes[i][0].plot(binned_x, target_ramp, label='Experiment')
                axes[i][0].plot(binned_x, model_ramp, label='Model')
                axes[i][1].set_ylabel('Change in synaptic weight (a.u.)')
                axes[i][1].set_xlabel('Location (cm)')
                axes[i][0].set_ylabel('Ramp amplitude (mV)')
                axes[i][0].set_xlabel('Location (cm)')
                axes[i][1].set_title('Induction: %s' % induction_key, fontsize=mpl.rcParams['font.size'])
                axes[i][0].legend(loc='best', frameon=False, framealpha=0.5)
                ymin = min(ymin, np.min(model_ramp) - 1., np.min(target_ramp) - 1.)
                ymax = max(ymax, np.max(model_ramp) + 1., np.max(target_ramp) + 1.)
            for induction_key in f[description][cell_key]:
                i = int(induction_key) - 1
                axes[i][0].set_ylim([ymin, ymax])
                axes[i][1].set_ylim([-peak_weight, peak_weight])
            clean_axes(axes)
            fig.suptitle('Cell_id: %s' % cell_key, fontsize=mpl.rcParams['font.size'])
            fig.tight_layout()
            plt.show()
            plt.close()
    mpl.rcParams['font.size'] = orig_fontsize


def plot_exported_BTSP_model_ramp_features1(processed_export_file_path):
    """

    :param processed_export_file_path: str (path)
    """
    orig_fontsize = mpl.rcParams['font.size']
    # mpl.rcParams['font.size'] = 20.
    with h5py.File(processed_export_file_path, 'r') as f:
        description = 'model_ramp_context'
        group = f[description]
        peak_locs = group['peak_locs'][:]
        binned_x = group['binned_x'][:]
        local_signal_range = group['local_signal_range'][:]
        norm_depot_rate = group['norm_depot_rate'][:]
        peak_weight = group.attrs['peak_weight']
        fig, axes = plt.subplots(1)
        axes.plot(local_signal_range, norm_depot_rate)
        axes.set_xlabel('Local signal amplitude (a.u.)')
        axes.set_ylabel('Normalized rate')
        axes.set_title('Depotentiation rate')
        clean_axes(axes)
        fig.tight_layout()
        plt.show()
        plt.close()
        description = 'model_ramp_features'
        for cell_key in f[description]:
            fig, axes = plt.subplots(2, 2)
            ymin = -1.
            ymax = 10.
            for induction_key in f[description][cell_key]:
                i = int(induction_key) - 1
                group = f[description][cell_key][induction_key]
                weights = group['weights'][:]
                initial_weights = group['initial_weights'][:]
                delta_weights = np.subtract(weights, initial_weights)
                target_ramp = group['target_ramp'][:]
                model_ramp = group['model_ramp'][:]
                axes[i][1].plot(peak_locs, delta_weights)
                axes[i][0].plot(binned_x, target_ramp, label='Experiment')
                axes[i][0].plot(binned_x, model_ramp, label='Model')
                axes[i][1].set_ylabel('Change in synaptic weight (a.u.)')
                axes[i][1].set_xlabel('Location (cm)')
                axes[i][0].set_ylabel('Ramp amplitude (mV)')
                axes[i][0].set_xlabel('Location (cm)')
                axes[i][1].set_title('Induction: %s' % induction_key, fontsize=mpl.rcParams['font.size'])
                axes[i][0].legend(loc='best', frameon=False, framealpha=0.5)
                ymin = min(ymin, np.min(model_ramp) - 1., np.min(target_ramp) - 1.)
                ymax = max(ymax, np.max(model_ramp) + 1., np.max(target_ramp) + 1.)
            for induction_key in f[description][cell_key]:
                i = int(induction_key) - 1
                axes[i][0].set_ylim([ymin, ymax])
                axes[i][1].set_ylim([-peak_weight, peak_weight])
            clean_axes(axes)
            fig.suptitle('Cell_id: %s' % cell_key, fontsize=mpl.rcParams['font.size'])
            fig.tight_layout()
            plt.show()
            plt.close()
    mpl.rcParams['font.size'] = orig_fontsize