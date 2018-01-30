__author__ = 'milsteina'
from function_lib import *
from matplotlib import cm
import matplotlib.lines as mlines
import matplotlib as mpl
import numpy as np
import scipy.signal as signal
import scipy.stats as stats

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.size'] = 14.
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


def plot_exported_BTSP_model_ramp_features5(processed_export_file_path):
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
        signal_xrange = group['signal_xrange'][:]
        depot_rate = group['depot_rate'][:]
        pot_rate = group['pot_rate'][:]
        peak_weight = group.attrs['peak_weight']
        fig, axes = plt.subplots(1)
        axes.plot(signal_xrange, pot_rate, label='Potentiation rate')
        axes.plot(signal_xrange, depot_rate, label='Depotentiation rate')
        axes.set_xlabel('Normalized plasticity signal amplitude (a.u.)')
        axes.set_ylabel('Normalized rate')
        axes.set_title('Plasticity signal transformations')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes)
        fig.tight_layout()
        # plt.show()
        # plt.close()
        description = 'model_ramp_features'
        for cell_key in f[description]:
            fig, axes = plt.subplots(2, 2)
            ymin = -1.
            ymax = 10.
            for induction_key in f[description][cell_key]:
                i = int(float(induction_key)) - 1
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
                axes[i][1].set_title('Induction: %i' % int(float(induction_key)), fontsize=mpl.rcParams['font.size'])
                axes[i][0].legend(loc='best', frameon=False, framealpha=0.5)
                ymin = min(ymin, np.min(model_ramp) - 1., np.min(target_ramp) - 1.)
                ymax = max(ymax, np.max(model_ramp) + 1., np.max(target_ramp) + 1.)
            for induction_key in f[description][cell_key]:
                i = int(float(induction_key)) - 1
                axes[i][0].set_ylim([ymin, ymax])
                axes[i][1].set_ylim([-peak_weight, peak_weight])
            clean_axes(axes)
            fig.suptitle('Cell_id: %i' % int(float(cell_key)), fontsize=mpl.rcParams['font.size'])
            fig.tight_layout()
            plt.show()
            plt.close()
    mpl.rcParams['font.size'] = orig_fontsize


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
        dual_signal_product_range = group['dual_signal_product_range'][:]
        norm_depot_rate = group['norm_depot_rate'][:]
        peak_weight = group.attrs['peak_weight']
        fig, axes = plt.subplots(1)
        axes.plot(dual_signal_product_range, norm_depot_rate)
        axes.set_xlabel('Plasticity signal amplitude (a.u.)')
        axes.set_ylabel('Normalized rate')
        axes.set_title('Depotentiation rate')
        clean_axes(axes)
        fig.tight_layout()
        # plt.show()
        # plt.close()
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
                i = int(float(induction_key)) - 1
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
                i = int(float(induction_key)) - 1
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
                i = int(float(induction_key)) - 1
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
    mpl.rcParams['font.size'] = 12.
    with h5py.File(processed_export_file_path, 'r') as f:
        description = 'model_ramp_context'
        group = f[description]
        peak_locs = group['peak_locs'][:]
        binned_x = group['binned_x'][:]
        # local_signal_range = group['skewnorm_range'][:]
        local_signal_range = group['local_signal_range'][:]
        norm_depot_rate = group['norm_depot_rate'][:]
        peak_weight = group.attrs['peak_weight']
        """
        fig, axes = plt.subplots(1)
        axes.plot(local_signal_range, norm_depot_rate)
        axes.set_xlabel('Local signal amplitude (a.u.)')
        axes.set_ylabel('Normalized rate')
        axes.set_title('Depotentiation rate')
        clean_axes(axes)
        fig.tight_layout()
        plt.show()
        plt.close()
        """
        description = 'model_ramp_features'
        for cell_key in f[description]:
            fig, axes = plt.subplots(4, 2)
            ybounds = [(0., 0.), (0., peak_weight * 1.2), (-peak_weight, peak_weight * 1.2), (0., 0.)]
            ybar = [0., 1.1 * peak_weight, 1.1 * peak_weight, 0.]
            this_signal_max = max(
                max([np.max(f[description][cell_key][induction_key]['example_local_signal'][:]) for
                 induction_key in f[description][cell_key]]),
                max([np.max(f[description][cell_key][induction_key]['global_signal'][:]) for
                 induction_key in f[description][cell_key]]))
            ybounds[0] = (0., 1.2 * this_signal_max)
            ybar[0] = 1.1 * this_signal_max
            this_ramp_max = max(
                max([np.max(f[description][cell_key][induction_key]['target_ramp'][:]) for
                 induction_key in f[description][cell_key]]),
                max([np.max(f[description][cell_key][induction_key]['model_ramp'][:]) for
                 induction_key in f[description][cell_key]]))
            this_ramp_min = min(
                min([np.min(f[description][cell_key][induction_key]['target_ramp'][:]) for
                 induction_key in f[description][cell_key]]),
                min([np.min(f[description][cell_key][induction_key]['model_ramp'][:]) for
                 induction_key in f[description][cell_key]]))
            ybounds[3] = (this_ramp_min - 0.1 * this_ramp_max, 1.2 * this_ramp_max)
            ybar[3] = 1.1 * this_ramp_max
            for induction_key in f[description][cell_key]:
                i = int(induction_key) - 1
                group = f[description][cell_key][induction_key]
                weights = group['weights'][:]
                initial_weights = group['initial_weights'][:]
                delta_weights = np.subtract(weights, initial_weights)
                target_ramp = group['target_ramp'][:]
                model_ramp = group['model_ramp'][:]
                example_local_signal = group['example_local_signal'][:]
                global_signal = group['global_signal'][:]
                down_t = group['down_t'][:] / 1000.  # ms to s
                example_weight_dynamics = group['example_weight_dynamics'][:]
                mean_induction_start_loc = group.attrs['mean_induction_start_loc']
                mean_induction_stop_loc = group.attrs['mean_induction_stop_loc']
                induction_start_times = group.attrs['induction_start_times']
                induction_stop_times = group.attrs['induction_stop_times']
                track_start_times = group.attrs['track_start_times']
                track_stop_times = group.attrs['track_stop_times']
                axes[0][i].plot(down_t, example_local_signal, c='k', label='Local synaptic signal')
                axes[0][i].plot(down_t, global_signal, c='c', label='Global dendritic signal')
                axes[0][i].set_xlim([track_start_times[0] / 1000., track_stop_times[0] / 1000.])
                axes[0][i].set_ylim(ybounds[0])
                axes[0][i].hlines([ybar[0]] * len(induction_start_times),
                               xmin=induction_start_times / 1000.,
                               xmax=induction_stop_times / 1000., linewidth=3, colors='r')
                axes[0][i].set_xlabel('Time (s)')
                axes[0][i].set_ylabel('Plasticity signal\namplitude')
                axes[1][i].plot(down_t, example_weight_dynamics, c='k')
                track_stop_index = np.where(down_t >= track_stop_times[0] / 1000.)[0][0]
                ymax1 = np.max(example_weight_dynamics[:track_stop_index])
                ymin1 = np.min(example_weight_dynamics[:track_stop_index])
                ydeltamin1 = ymax1 - ymin1
                ybar1 = ymax1 + 0.1 * ydeltamin1
                axes[1][i].set_ylim([ymin1 - 0.1 * ydeltamin1, ymax1 + 0.2 * ydeltamin1])
                axes[1][i].set_xlim([track_start_times[0] / 1000., track_stop_times[0] / 1000.])
                axes[1][i].hlines([ybar1] * len(induction_start_times),
                               xmin=induction_start_times / 1000.,
                               xmax=induction_stop_times / 1000., linewidth=3, colors='r')
                axes[1][i].set_ylabel('Synaptic weight\n(single input)')
                axes[1][i].set_xlabel('Time (s)')
                axes[2][i].plot(peak_locs, delta_weights, c='k')
                axes[2][i].hlines(ybar[2], xmin=mean_induction_start_loc, xmax=mean_induction_stop_loc, linewidth=3,
                                  colors='r')
                axes[3][i].plot(binned_x, target_ramp, label='Experiment', c='grey')
                axes[3][i].plot(binned_x, model_ramp, label='Model', c='k')
                axes[3][i].hlines(ybar[3], xmin=mean_induction_start_loc, xmax=mean_induction_stop_loc, linewidth=3,
                                  colors='r')
                axes[2][i].set_ylabel('Change in\nsynaptic weight\n(input population)')
                axes[2][i].set_xlabel('Peak location of presynaptic input (cm)')
                axes[3][i].set_ylabel('Place field\namplitude (mV)')
                axes[3][i].set_xlabel('Location (cm)')
                # axes[2][i].set_title('Induction: %s' % induction_key, fontsize=mpl.rcParams['font.size'])
                axes[3][i].set_ylim(ybounds[3])
                axes[2][i].set_ylim(ybounds[2])
            axes[0][0].set_title('First place field induction')
            axes[0][1].set_title('Second place field induction')
            handles, labels = axes[0][0].get_legend_handles_labels()
            handles.append(plt.Line2D((0,1),(0,0), color='r'))
            labels.append('Dendritic plateau potential')
            axes[0][0].legend(handles, labels, loc='best', frameon=False, framealpha=0.5, handlelength=1)
            axes[3][0].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
            clean_axes(axes)
            # fig.suptitle('Cell_id: %s' % cell_key, fontsize=mpl.rcParams['font.size'])
            fig.tight_layout()
            plt.subplots_adjust(hspace=0.7, wspace=0.5)
            plt.show()
            plt.close()
    mpl.rcParams['font.size'] = orig_fontsize

import nested
storage = nested.PopulationStorage()
storage.get_best()

"""
sig = lambda baseline, slope, xhalf, tau, signal: baseline - slope / (1. + np.exp(xhalf / tau)) + \
                                                      slope / (1. + np.exp((xhalf - signal) / tau))
"""