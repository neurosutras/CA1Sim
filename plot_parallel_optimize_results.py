__author__ = 'milsteina'
from BTSP_utils import *
from nested.optimize_utils import *
from matplotlib import cm
import matplotlib.lines as mlines
import matplotlib as mpl

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


def plot_BTSP_model_ramp_traces_9(file_path, self_consistent=False, savefig_dir=None, show=True,
                                  full_output=False):
    """

    :param file_path: str (path)
    :param self_consistent: bool : whether or not to also plot model_ramp_self_consistent_features when available
    :param savefig_dir: str : path to output figures
    :param show: bool : actually display plots?
    :param full_output : return summary dicts
    return dict (optional)
    """
    date_stamp = datetime.datetime.today().strftime('%Y%m%d%H%M')
    ramp_amp, ramp_width, peak_shift, ratio, min_val = {}, {}, {}, {}, {}
    orig_fontsize = mpl.rcParams['font.size']
    # mpl.rcParams['font.size'] = 24.
    with h5py.File(file_path, 'r') as f:
        shared_context_key = 'shared_context'
        group = f[shared_context_key]
        peak_locs = group['peak_locs'][:]
        binned_x = group['binned_x'][:]
        signal_xrange = group['signal_xrange'][:]
        param_names = group['param_names'][:]
        exported_data_key = 'exported_data'
        for cell_key in f[exported_data_key]:
            for parameter in ramp_amp, ramp_width, peak_shift, ratio, min_val:
                if cell_key not in parameter:
                    parameter[cell_key] = {}
            description = 'model_ramp_features'
            group = f[exported_data_key][cell_key].itervalues().next()[description]
            param_array = group['param_array'][:]
            param_dict = param_array_to_dict(param_array, param_names)
            peak_weight = param_dict['peak_delta_weight'] + 1.
            depot_rate = group['depot_rate'][:]
            pot_rate = group['pot_rate'][:]
            local_signal_filter_t = group['local_signal_filter_t'][:]
            local_signal_filter = group['local_signal_filter'][:]
            global_filter_t = group['global_filter_t'][:]
            global_filter = group['global_filter'][:]

            fig, axes = plt.subplots(1, figsize=[8., 6.])
            axes.plot(signal_xrange, pot_rate, label='Potentiation rate', c='k')
            axes.plot(signal_xrange, depot_rate, label='Depotentiation rate', c='r')
            axes.set_xlabel('Normalized signal amplitude (a.u.)')
            axes.set_ylabel('Normalized rate')
            axes.set_title('Plasticity signal transformations', fontsize=mpl.rcParams['font.size'])
            axes.legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
            clean_axes(axes)
            fig.tight_layout()

            fig2, axes2 = plt.subplots(1, figsize=[8., 6.])
            axes2.plot(local_signal_filter_t / 1000., local_signal_filter / np.max(local_signal_filter), color='r',
                       label='Local plasticity signal filter')
            axes2.plot(global_filter_t / 1000., global_filter / np.max(global_filter), color='k',
                       label='Global plasticity signal filter')
            axes2.set_xlabel('Time (s)')
            axes2.set_ylabel('Normalized\nfilter amplitude')
            axes2.set_title('Plasticity signal filters', fontsize=mpl.rcParams['font.size'])
            axes2.legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
            axes2.set_xlim(-0.5, max(5000., local_signal_filter_t[-1], global_filter_t[-1]) / 1000.)
            clean_axes(axes2)
            fig2.tight_layout()

            for induction_key in f[exported_data_key][cell_key]:
                for parameter in ramp_amp, ramp_width, peak_shift, ratio, min_val:
                    if induction_key not in parameter[cell_key]:
                        parameter[cell_key][induction_key] = {}
                i = int(float(induction_key)) - 1
                description = 'model_ramp_features'
                if i + 1 == 2 and self_consistent and \
                        'model_ramp_self_consistent_features' in f[exported_data_key][cell_key][induction_key]:
                    description = 'model_ramp_self_consistent_features'
                group = f[exported_data_key][cell_key][induction_key][description]

                ramp_amp[cell_key][induction_key]['target'] = group.attrs['target_ramp_amp']
                ramp_width[cell_key][induction_key]['target'] = group.attrs['target_ramp_width']
                peak_shift[cell_key][induction_key]['target'] = group.attrs['target_peak_shift']
                ratio[cell_key][induction_key]['target'] = group.attrs['target_ratio']
                min_val[cell_key][induction_key]['target'] = group.attrs['target_min_val']
                ramp_amp[cell_key][induction_key]['model'] = group.attrs['model_ramp_amp']
                ramp_width[cell_key][induction_key]['model'] = group.attrs['model_ramp_width']
                peak_shift[cell_key][induction_key]['model'] = group.attrs['model_peak_shift']
                ratio[cell_key][induction_key]['model'] = group.attrs['model_ratio']
                min_val[cell_key][induction_key]['model'] = group.attrs['model_min_val']

                induction_start_times = group.attrs['induction_start_times']
                induction_stop_times = group.attrs['induction_stop_times']
                down_t = group['down_t'][:]
                example_local_signal = group['example_local_signal'][:]
                global_signal = group['global_signal'][:]
                example_weight_dynamics = group['example_weight_dynamics'][:]

                fig4, axes4 = plt.subplots(2, sharex=True, figsize=[8., 8.])
                ymax0 = max(np.max(example_local_signal), np.max(global_signal))
                bar_loc0 = ymax0 * 1.05
                axes4[0].plot(down_t / 1000., example_local_signal, c='r', label='Local plasticity signal')
                axes4[0].plot(down_t / 1000., global_signal, c='k', label='Global signal')
                axes4[0].set_ylim([-0.1 * ymax0, 1.1 * ymax0])
                axes4[0].hlines([bar_loc0] * len(induction_start_times),
                                xmin=induction_start_times / 1000.,
                                xmax=induction_stop_times / 1000., linewidth=2)
                axes4[0].set_xlabel('Time (s)')
                axes4[0].set_ylabel('Plasticity\nsignal\namplitudes')
                axes4[0].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
                axes4[1].plot(down_t / 1000., example_weight_dynamics)
                axes4[1].set_ylim([0., peak_weight * 1.1])
                axes4[1].hlines([peak_weight * 1.05] * len(induction_start_times),
                                xmin=induction_start_times / 1000.,
                                xmax=induction_stop_times / 1000., linewidth=2)
                axes4[1].set_ylabel('Synaptic weight\n(example\nsingle input)')
                axes4[1].set_xlabel('Time (s)')
                fig4.suptitle('Cell: %i, Induction: %i' % (int(float(cell_key)), int(float(induction_key))),
                              fontsize=mpl.rcParams['font.size'])
                clean_axes(axes4)
                fig4.tight_layout()
                fig4.subplots_adjust(top=0.9, hspace=0.5)

            fig3, axes3 = plt.subplots(2, 3, figsize=[16, 6])
            ymin = -1.
            ymax = 10.
            for induction_key in f[exported_data_key][cell_key]:
                i = int(float(induction_key)) - 1
                description = 'model_ramp_features'
                if i + 1 == 2 and self_consistent and \
                        'model_ramp_self_consistent_features' in f[exported_data_key][cell_key][induction_key]:
                    description = 'model_ramp_self_consistent_features'
                group = f[exported_data_key][cell_key][induction_key][description]
                model_weights = group['model_weights'][:]
                initial_weights = group['initial_weights'][:]
                delta_weights = np.subtract(model_weights, initial_weights)
                initial_model_ramp = group['initial_model_ramp'][:]
                target_ramp = group['target_ramp'][:]
                model_ramp = group['model_ramp'][:]
                ymin = min(ymin, np.min(model_ramp) - 1., np.min(target_ramp) - 1.)
                ymax = max(ymax, np.max(model_ramp) + 1., np.max(target_ramp) + 1.)
                if 'initial_exp_ramp' in group:
                    initial_exp_ramp = group['initial_exp_ramp'][:]
                    axes3[i][0].plot(binned_x, initial_exp_ramp, label='Before', c='darkgrey')
                    ymin = min(ymin, np.min(initial_exp_ramp) - 1.)
                    ymax = max(ymax, np.max(initial_exp_ramp) + 1.)
                axes3[i][0].plot(binned_x, target_ramp, label='After', c='r')
                axes3[i][0].set_title('Induction %i\nExperiment Vm:' % (i + 1), fontsize=mpl.rcParams['font.size'])
                axes3[i][1].plot(binned_x, initial_model_ramp, label='Before', c='darkgrey')
                axes3[i][1].plot(binned_x, model_ramp, label='After', c='c')
                axes3[i][1].set_title('\nModel Vm', fontsize=mpl.rcParams['font.size'])
                axes3[i][2].plot(peak_locs, delta_weights, c='k')
                axes3[i][2].set_title('Change in\nSynaptic Weights', fontsize=mpl.rcParams['font.size'])

                axes3[i][0].set_xlabel('Location (cm)')
                axes3[i][0].set_ylabel('Depolarization\namplitude (mV)')
                axes3[i][1].set_xlabel('Location (cm)')
                axes3[i][1].set_ylabel('Depolarization\namplitude (mV)')
                axes3[i][2].set_xlabel('Location (cm)')
                axes3[i][2].set_ylabel('Change in synaptic\nweight (a.u.)')
                xmin, xmax = axes3[i][2].get_xlim()
                axes3[i][2].plot([xmin, xmax], [0., 0.], c='darkgrey', alpha=0.5, ls='--')

            for induction_key in f[exported_data_key][cell_key]:
                i = int(float(induction_key)) - 1
                description = 'model_ramp_features'
                if i + 1 == 2 and self_consistent and \
                        'model_ramp_self_consistent_features' in f[exported_data_key][cell_key][induction_key]:
                    description = 'model_ramp_self_consistent_features'
                group = f[exported_data_key][cell_key][induction_key][description]
                mean_induction_start_loc = group.attrs['mean_induction_start_loc']
                mean_induction_stop_loc = group.attrs['mean_induction_stop_loc']
                axes3[i][0].set_ylim([ymin, ymax * 1.05])
                axes3[i][1].set_ylim([ymin, ymax * 1.05])
                axes3[i][2].set_ylim([-peak_weight, peak_weight * 1.1])
                axes3[i][0].hlines(ymax, xmin=mean_induction_start_loc, xmax=mean_induction_stop_loc)
                axes3[i][1].hlines(ymax, xmin=mean_induction_start_loc, xmax=mean_induction_stop_loc)
                axes3[i][2].hlines(peak_weight * 1.05, xmin=mean_induction_start_loc, xmax=mean_induction_stop_loc)
                axes3[i][0].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
                axes3[i][1].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
            clean_axes(axes3)
            fig3.suptitle('Cell: %i' % int(float(cell_key)), fontsize=mpl.rcParams['font.size'])
            fig3.tight_layout()
            fig3.subplots_adjust(top=0.9, hspace=0.75)
            if savefig_dir is not None:
                if not os.path.isdir(savefig_dir):
                    raise IOError('plot_BTSP_model_ramp_traces: cannot save figure to invalid directory: %s' %
                                  savefig_dir)
                else:
                    fig3.savefig('%s/%s_BTSP2_cell_%s_ramp_shape_traces.svg' %
                                 (savefig_dir, date_stamp, cell_key), format='svg')
            if show:
                plt.show()
    plt.close('all')
    mpl.rcParams['font.size'] = orig_fontsize
    plot_BTSP_model_ramp_summary(ramp_amp, ramp_width, peak_shift, ratio, min_val, date_stamp=date_stamp,
                                 savefig_dir=savefig_dir, show=show)
    if full_output:
        return ramp_amp, ramp_width, peak_shift, ratio, min_val


def plot_BTSP_model_ramp_traces_8(file_path, self_consistent=False, savefig_dir=None, show=True,
                                             full_output=False):
    """

    :param file_path: str (path)
    :param self_consistent: bool : whether or not to also plot model_ramp_self_consistent_features when available
    :param savefig_dir: str : path to output figures
    :param show: bool : actually display plots?
    :param full_output : return summary dicts
    return dict (optional)
    """
    date_stamp = datetime.datetime.today().strftime('%Y%m%d%H%M')
    ramp_amp, ramp_width, peak_shift, ratio, min_val = {}, {}, {}, {}, {}
    orig_fontsize = mpl.rcParams['font.size']
    mpl.rcParams['font.size'] = 24.
    with h5py.File(file_path, 'r') as f:
        shared_context_key = 'shared_context'
        group = f[shared_context_key]
        peak_locs = group['peak_locs'][:]
        binned_x = group['binned_x'][:]
        signal_xrange = group['signal_xrange'][:]
        param_names = group['param_names'][:]
        exported_data_key = 'exported_data'
        for cell_key in f[exported_data_key]:
            for parameter in ramp_amp, ramp_width, peak_shift, ratio, min_val:
                if cell_key not in parameter:
                    parameter[cell_key] = {}
            description = 'model_ramp_features'
            group = f[exported_data_key][cell_key].itervalues().next()[description]
            param_array = group['param_array'][:]
            param_dict = param_array_to_dict(param_array, param_names)
            peak_weight = param_dict['peak_delta_weight'] + 1.
            depot_rate = group['depot_rate'][:]
            pot_rate = group['pot_rate'][:]
            local_pot_filter_t = group['local_pot_filter_t'][:]
            local_pot_filter = group['local_pot_filter'][:]
            local_depot_filter_t = group['local_depot_filter_t'][:]
            local_depot_filter = group['local_depot_filter'][:]
            global_filter_t = group['global_filter_t'][:]
            global_filter = group['global_filter'][:]
            
            fig, axes = plt.subplots(1, figsize=[8., 6.])
            axes.plot(signal_xrange, pot_rate, label='Potentiation rate', c='k')
            axes.plot(signal_xrange, depot_rate, label='Depotentiation rate', c='r')
            axes.set_xlabel('Normalized signal amplitude (a.u.)')
            axes.set_ylabel('Normalized rate')
            axes.set_title('Plasticity signal transformations', fontsize=mpl.rcParams['font.size'])
            axes.legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
            clean_axes(axes)
            fig.tight_layout()

            fig2, axes2 = plt.subplots(1, figsize=[8., 6.])
            axes2.plot(local_pot_filter_t / 1000., local_pot_filter / np.max(local_pot_filter), color='r',
                      label='Local potentiation signal filter')
            axes2.plot(local_depot_filter_t / 1000., local_depot_filter / np.max(local_depot_filter), color='c',
                      label='Local de-potentiation signal filter')
            axes2.plot(global_filter_t / 1000., global_filter / np.max(global_filter), color='k',
                      label='Global signal filter')
            axes2.set_xlabel('Time (s)')
            axes2.set_ylabel('Normalized\nfilter amplitude')
            axes2.set_title('Plasticity signal filters', fontsize=mpl.rcParams['font.size'])
            axes2.legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
            axes2.set_xlim(-0.5, max(5000., local_pot_filter_t[-1], global_filter_t[-1]) / 1000.)
            clean_axes(axes2)
            fig2.tight_layout()

            for induction_key in f[exported_data_key][cell_key]:
                for parameter in ramp_amp, ramp_width, peak_shift, ratio, min_val:
                    if induction_key not in parameter[cell_key]:
                        parameter[cell_key][induction_key] = {}
                i = int(float(induction_key)) - 1
                description = 'model_ramp_features'
                if i + 1 == 2 and self_consistent and \
                        'model_ramp_self_consistent_features' in f[exported_data_key][cell_key][induction_key]:
                    description = 'model_ramp_self_consistent_features'
                group = f[exported_data_key][cell_key][induction_key][description]

                ramp_amp[cell_key][induction_key]['target'] = group.attrs['target_ramp_amp']
                ramp_width[cell_key][induction_key]['target'] = group.attrs['target_ramp_width']
                peak_shift[cell_key][induction_key]['target'] = group.attrs['target_peak_shift']
                ratio[cell_key][induction_key]['target'] = group.attrs['target_ratio']
                min_val[cell_key][induction_key]['target'] = group.attrs['target_min_val']
                ramp_amp[cell_key][induction_key]['model'] = group.attrs['model_ramp_amp']
                ramp_width[cell_key][induction_key]['model'] = group.attrs['model_ramp_width']
                peak_shift[cell_key][induction_key]['model'] = group.attrs['model_peak_shift']
                ratio[cell_key][induction_key]['model'] = group.attrs['model_ratio']
                min_val[cell_key][induction_key]['model'] = group.attrs['model_min_val']
                
                induction_start_times = group.attrs['induction_start_times']
                induction_stop_times = group.attrs['induction_stop_times']
                down_t = group['down_t'][:]
                example_local_pot_signal = group['example_local_pot_signal'][:]
                example_local_depot_signal = group['example_local_depot_signal'][:]
                global_signal = group['global_signal'][:]
                example_weight_dynamics = group['example_weight_dynamics'][:]

                fig4, axes4 = plt.subplots(3, sharex=True, figsize=[8., 8.])
                ymax0 = max(np.max(example_local_pot_signal), np.max(global_signal))
                bar_loc0 = ymax0 * 1.05
                ymax1 = max(np.max(example_local_depot_signal), np.max(global_signal))
                bar_loc1 = ymax1 * 1.05
                axes4[0].plot(down_t / 1000., example_local_pot_signal, c='r', label='Local potentiation signal')
                axes4[0].plot(down_t / 1000., global_signal, c='k', label='Global signal')
                axes4[0].set_ylim([-0.1 * ymax0, 1.1 * ymax0])
                axes4[0].hlines([bar_loc0] * len(induction_start_times),
                               xmin=induction_start_times / 1000.,
                               xmax=induction_stop_times / 1000., linewidth=2)
                axes4[0].set_xlabel('Time (s)')
                axes4[0].set_ylabel('Plasticity\nsignal\namplitudes')
                axes4[0].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
                axes4[1].plot(down_t / 1000., example_local_depot_signal, c='c',
                             label='Local de-potentiation signal')
                axes4[1].plot(down_t / 1000., global_signal, c='k', label='Global signal')
                axes4[1].set_ylim([-0.1 * ymax1, 1.1 * ymax1])
                axes4[1].hlines([bar_loc1] * len(induction_start_times),
                               xmin=induction_start_times / 1000.,
                               xmax=induction_stop_times / 1000., linewidth=2)
                axes4[1].set_xlabel('Time (s)')
                axes4[1].set_ylabel('Plasticity\nsignal\namplitudes')
                axes4[1].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
                axes4[2].plot(down_t / 1000., example_weight_dynamics)
                axes4[2].set_ylim([0., peak_weight * 1.1])
                axes4[2].hlines([peak_weight * 1.05] * len(induction_start_times),
                               xmin=induction_start_times / 1000.,
                               xmax=induction_stop_times / 1000., linewidth=2)
                axes4[2].set_ylabel('Synaptic weight\n(example\nsingle input)')
                axes4[2].set_xlabel('Time (s)')
                fig4.suptitle('Cell: %i, Induction: %i' % (int(float(cell_key)), int(float(induction_key))),
                              fontsize=mpl.rcParams['font.size'])
                clean_axes(axes4)
                fig4.tight_layout()
                fig4.subplots_adjust(top=0.9, hspace=0.5)

            fig3, axes3 = plt.subplots(2, 3, figsize=[16, 6])
            ymin = -1.
            ymax = 10.
            for induction_key in f[exported_data_key][cell_key]:
                i = int(float(induction_key)) - 1
                description = 'model_ramp_features'
                if i + 1 == 2 and self_consistent and \
                        'model_ramp_self_consistent_features' in f[exported_data_key][cell_key][induction_key]:
                    description = 'model_ramp_self_consistent_features'
                group = f[exported_data_key][cell_key][induction_key][description]
                model_weights = group['model_weights'][:]
                initial_weights = group['initial_weights'][:]
                delta_weights = np.subtract(model_weights, initial_weights)
                initial_model_ramp = group['initial_model_ramp'][:]
                target_ramp = group['target_ramp'][:]
                model_ramp = group['model_ramp'][:]
                ymin = min(ymin, np.min(model_ramp) - 1., np.min(target_ramp) - 1.)
                ymax = max(ymax, np.max(model_ramp) + 1., np.max(target_ramp) + 1.)
                if 'initial_exp_ramp' in group:
                    initial_exp_ramp = group['initial_exp_ramp'][:]
                    axes3[i][0].plot(binned_x, initial_exp_ramp, label='Before', c='darkgrey')
                    ymin = min(ymin, np.min(initial_exp_ramp) - 1.)
                    ymax = max(ymax, np.max(initial_exp_ramp) + 1.)
                axes3[i][0].plot(binned_x, target_ramp, label='After', c='r')
                axes3[i][0].set_title('Induction %i\nExperiment Vm:' % (i + 1), fontsize=mpl.rcParams['font.size'])
                axes3[i][1].plot(binned_x, initial_model_ramp, label='Before', c='darkgrey')
                axes3[i][1].plot(binned_x, model_ramp, label='After', c='c')
                axes3[i][1].set_title('\nModel Vm', fontsize=mpl.rcParams['font.size'])
                axes3[i][2].plot(peak_locs, delta_weights, c='k')
                axes3[i][2].set_title('Change in\nSynaptic Weights', fontsize=mpl.rcParams['font.size'])

                axes3[i][0].set_xlabel('Location (cm)')
                axes3[i][0].set_ylabel('Depolarization\namplitude (mV)')
                axes3[i][1].set_xlabel('Location (cm)')
                axes3[i][1].set_ylabel('Depolarization\namplitude (mV)')
                axes3[i][2].set_xlabel('Location (cm)')
                axes3[i][2].set_ylabel('Change in synaptic\nweight (a.u.)')
                xmin, xmax = axes3[i][2].get_xlim()
                axes3[i][2].plot([xmin, xmax], [0., 0.], c='darkgrey', alpha=0.5, ls='--')

            for induction_key in f[exported_data_key][cell_key]:
                i = int(float(induction_key)) - 1
                description = 'model_ramp_features'
                if i + 1 == 2 and self_consistent and \
                        'model_ramp_self_consistent_features' in f[exported_data_key][cell_key][induction_key]:
                    description = 'model_ramp_self_consistent_features'
                group = f[exported_data_key][cell_key][induction_key][description]
                mean_induction_start_loc = group.attrs['mean_induction_start_loc']
                mean_induction_stop_loc = group.attrs['mean_induction_stop_loc']
                axes3[i][0].set_ylim([ymin, ymax * 1.05])
                axes3[i][1].set_ylim([ymin, ymax * 1.05])
                axes3[i][2].set_ylim([-peak_weight, peak_weight * 1.1])
                axes3[i][0].hlines(ymax, xmin=mean_induction_start_loc, xmax=mean_induction_stop_loc)
                axes3[i][1].hlines(ymax, xmin=mean_induction_start_loc, xmax=mean_induction_stop_loc)
                axes3[i][2].hlines(peak_weight * 1.05, xmin=mean_induction_start_loc, xmax=mean_induction_stop_loc)
                axes3[i][0].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
                axes3[i][1].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
            clean_axes(axes3)
            fig3.suptitle('Cell: %i' % int(float(cell_key)), fontsize=mpl.rcParams['font.size'])
            fig3.tight_layout()
            fig3.subplots_adjust(top=0.9, hspace=0.75)
            if savefig_dir is not None:
                if not os.path.isdir(savefig_dir):
                    raise IOError('plot_BTSP_model_ramp_traces: cannot save figure to invalid directory: %s' %
                                  savefig_dir)
                else:
                    fig3.savefig('%s/%s_BTSP2_cell_%s_ramp_shape_traces.svg' %
                                 (savefig_dir, date_stamp, cell_key), format='svg')
            if show:
                plt.show()
    plt.close('all')
    mpl.rcParams['font.size'] = orig_fontsize
    plot_BTSP_model_ramp_summary(ramp_amp, ramp_width, peak_shift, ratio, min_val, date_stamp=date_stamp,
                                 savefig_dir=savefig_dir, show=show)
    if full_output:
        return ramp_amp, ramp_width, peak_shift, ratio, min_val


def plot_BTSP_model_ramp_summary(ramp_amp, ramp_width, peak_shift, ratio, min_val, date_stamp=None, savefig_dir=None,
                                 show=True):
    """

    :param ramp_amp: dict
    :param ramp_width: dict
    :param peak_shift: dict
    :param ratio: dict
    :param min_val: dict
    :param date_stamp: str
    :param savefig_dir: str: dir
    :param show: bool
    """
    orig_fontsize = mpl.rcParams['font.size']
    # mpl.rcParams['font.size'] = 24.
    if date_stamp is None:
        date_stamp = datetime.datetime.today().strftime('%Y%m%d%H%M')
    fig, axes = plt.subplots(2, 3, figsize=[12., 8.])
    tick_locs = [np.arange(0., 13., 3.), np.arange(0., 181., 60.), np.arange(-90., 91., 45.), np.arange(0., 6., 1.),
                 np.arange(0., 6., 1.)]
    for i, (parameter, label) in enumerate(zip([ramp_amp, ramp_width, peak_shift, ratio, min_val],
                                ['Ramp peak\namplitude (mV)', 'Ramp width (cm)', 'Peak shift (cm)',
                                 'Ramp asymmetry\n(ratio)', 'Ramp minimum\namplitude (mV)'])):
        col = i / 3
        row = i % 3
        axes[col][row].scatter(*zip(*[(parameter[cell_key]['1']['target'], parameter[cell_key]['1']['model'])
                            for cell_key in (cell_key for cell_key in parameter if '1' in parameter[cell_key])]),
                     color='darkgrey', alpha=0.5, label='Induction 1')
        axes[col][row].scatter(*zip(*[(parameter[cell_key]['2']['target'], parameter[cell_key]['2']['model'])
                            for cell_key in (cell_key for cell_key in parameter if '2' in parameter[cell_key])]),
                     color='r', alpha=0.5, label='Induction 2')
        axes[col][row].set_title(label, fontsize=mpl.rcParams['font.size'])
        if i == 0:
            axes[col][row].legend(loc='best', frameon=False, framealpha=0.5)
        ymin, ymax = axes[col][row].get_ylim()
        xmin, xmax = axes[col][row].get_xlim()
        min_lim = min(0., ymin, xmin)
        max_lim = max(0., ymax, xmax)
        axes[col][row].set_xlim([min_lim, max_lim])
        axes[col][row].set_ylim([min_lim, max_lim])
        axes[col][row].set_xticks(tick_locs[i])
        axes[col][row].set_yticks(tick_locs[i])
        axes[col][row].plot([min_lim, max_lim], [min_lim, max_lim], c='darkgrey', alpha=0.5, ls='--')
        axes[col][row].set_aspect('equal')
        axes[col][row].set_xlabel('Experiment')
        axes[col][row].set_ylabel('Model')
    clean_axes(axes)
    fig.tight_layout(h_pad=0.2)
    fig.subplots_adjust(hspace=0.75)
    if savefig_dir is not None:
        if not os.path.isdir(savefig_dir):
            raise IOError('plot_BTSP_model_ramp_traces: cannot save figure to invalid directory: %s' %
                          savefig_dir)
        else:
            fig.savefig('%s/%s_BTSP2_ramp_shape_summary.svg' % (savefig_dir, date_stamp), format='svg')
    if show:
        plt.show()
    plt.close('all')
    mpl.rcParams['font.size'] = orig_fontsize


def plot_exported_BTSP_model_ramp_features8(export_file_path):
    """

    :param export_file_path: str (path)
    """
    orig_fontsize = mpl.rcParams['font.size']
    # mpl.rcParams['font.size'] = 20.
    with h5py.File(export_file_path, 'r') as f:
        description = 'shared_context'
        group = f[description]
        peak_locs = group['peak_locs'][:]
        binned_x = group['binned_x'][:]
        signal_xrange = group['signal_xrange'][:]
        depot_rate = group['depot_rate'][:]
        pot_rate = group['pot_rate'][:]
        peak_weight = group.attrs['peak_weight']
        fig, axes = plt.subplots(1)
        axes.plot(signal_xrange, pot_rate, label='Potentiation rate', c='k')
        axes.plot(signal_xrange, depot_rate, label='Depotentiation rate', c='r')
        axes.set_xlabel('Normalized plasticity signal amplitude (a.u.)')
        axes.set_ylabel('Normalized rate')
        axes.set_title('Plasticity signal transformations')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes)
        fig.tight_layout()
        # plt.show()
        # plt.close()
        exported_data_key = 'exported_data'
        for cell_key in f[exported_data_key]:
            fig2, axes = plt.subplots(2, 2)
            ymin = -1.
            ymax = 10.
            for induction_key in f[exported_data_key][cell_key]:
                i = int(float(induction_key)) - 1
                description = 'model_ramp_features'
                group = f[exported_data_key][cell_key][induction_key][description]
                model_weights = group['model_weights'][:]
                initial_weights = group['initial_weights'][:]
                delta_weights = np.subtract(model_weights, initial_weights)
                target_ramp = group['target_ramp'][:]
                model_ramp = group['model_ramp'][:]
                axes[i][1].plot(peak_locs, delta_weights, c='r')
                axes[i][0].plot(binned_x, target_ramp, label='Experiment', c='k')
                axes[i][0].plot(binned_x, model_ramp, label='Model', c='r')
                axes[i][1].set_ylabel('Change in synaptic weight (a.u.)')
                axes[i][1].set_xlabel('Location (cm)')
                axes[i][0].set_ylabel('Ramp amplitude (mV)')
                axes[i][0].set_xlabel('Location (cm)')
                axes[i][1].set_title('Induction: %i' % int(float(induction_key)), fontsize=mpl.rcParams['font.size'])
                ymin = min(ymin, np.min(model_ramp) - 1., np.min(target_ramp) - 1.)
                ymax = max(ymax, np.max(model_ramp) + 1., np.max(target_ramp) + 1.)
            description = 'model_ramp_self_consistent_features'
            if description in f[exported_data_key][cell_key][induction_key]:
                group = f[exported_data_key][cell_key][induction_key][description]
                self_consistent_model_ramp = group['model_ramp'][:]
                self_consistent_initial_weights = group['initial_weights'][:]
                self_consistent_model_weights = group['model_weights'][:]
                self_consistent_delta_weights = np.subtract(self_consistent_model_weights,
                                                            self_consistent_initial_weights)
                axes[1][1].plot(peak_locs, self_consistent_delta_weights, c='c')
                axes[1][0].plot(binned_x, self_consistent_model_ramp, label='Model (self-consistent)', c='c')
                ymin = min(ymin, np.min(self_consistent_model_ramp) - 1.)
                ymax = max(ymax, np.max(self_consistent_model_ramp) + 1.)
            for induction_key in f[exported_data_key][cell_key]:
                i = int(float(induction_key)) - 1
                axes[i][0].set_ylim([ymin, ymax])
                axes[i][1].set_ylim([-peak_weight, peak_weight])
                axes[i][0].legend(loc='best', frameon=False, framealpha=0.5)
            clean_axes(axes)
            fig2.suptitle('Cell_id: %i' % int(float(cell_key)), fontsize=mpl.rcParams['font.size'])
            fig2.tight_layout(h_pad=0.2)
            plt.show()
            plt.close()
    mpl.rcParams['font.size'] = orig_fontsize


def plot_exported_BTSP_model_ramp_features7(processed_export_file_path):
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
        axes.plot(signal_xrange, pot_rate, label='Potentiation rate', c='k')
        axes.plot(signal_xrange, depot_rate, label='Depotentiation rate', c='r')
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
                model_weights = group['model_weights'][:]
                initial_weights = group['initial_weights'][:]
                delta_weights = np.subtract(model_weights, initial_weights)
                target_ramp = group['target_ramp'][:]
                model_ramp = group['model_ramp'][:]
                axes[i][1].plot(peak_locs, delta_weights, c='r')
                axes[i][0].plot(binned_x, target_ramp, label='Experiment', c='k')
                axes[i][0].plot(binned_x, model_ramp, label='Model', c='r')
                axes[i][1].set_ylabel('Change in synaptic weight (a.u.)')
                axes[i][1].set_xlabel('Location (cm)')
                axes[i][0].set_ylabel('Ramp amplitude (mV)')
                axes[i][0].set_xlabel('Location (cm)')
                axes[i][1].set_title('Induction: %i' % int(float(induction_key)), fontsize=mpl.rcParams['font.size'])
                ymin = min(ymin, np.min(model_ramp) - 1., np.min(target_ramp) - 1.)
                ymax = max(ymax, np.max(model_ramp) + 1., np.max(target_ramp) + 1.)
            try:
                self_consistent_model_ramp = f['model_ramp_self_consistent_features'][cell_key]['2']['model_ramp'][:]
                self_consistent_initial_weights = \
                    f['model_ramp_features'][cell_key]['1']['model_weights'][:]
                self_consistent_model_weights = \
                    f['model_ramp_self_consistent_features'][cell_key]['2']['model_weights'][:]
                self_consistent_delta_weights = np.subtract(self_consistent_model_weights,
                                                            self_consistent_initial_weights)
                axes[1][1].plot(peak_locs, self_consistent_delta_weights, c='c')
                axes[1][0].plot(binned_x, self_consistent_model_ramp, label='Model (self-consistent)', c='c')
                ymin = min(ymin, np.min(self_consistent_model_ramp) - 1.)
                ymax = max(ymax, np.max(self_consistent_model_ramp) + 1.)
            except KeyError:
                pass
            for induction_key in f[description][cell_key]:
                i = int(float(induction_key)) - 1
                axes[i][0].set_ylim([ymin, ymax])
                axes[i][1].set_ylim([-peak_weight, peak_weight])
                axes[i][0].legend(loc='best', frameon=False, framealpha=0.5)
            clean_axes(axes)
            fig.suptitle('Cell_id: %i' % int(float(cell_key)), fontsize=mpl.rcParams['font.size'])
            fig.tight_layout()
            plt.show()
            plt.close()
    mpl.rcParams['font.size'] = orig_fontsize


def plot_exported_BTSP_model_ramp_features6(processed_export_file_path):
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
        axes.plot(signal_xrange, pot_rate, label='Potentiation rate', c='k')
        axes.plot(signal_xrange, depot_rate, label='Depotentiation rate', c='r')
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
                model_weights = group['model_weights'][:]
                initial_weights = group['initial_weights'][:]
                delta_weights = np.subtract(model_weights, initial_weights)
                target_ramp = group['target_ramp'][:]
                model_ramp = group['model_ramp'][:]
                axes[i][1].plot(peak_locs, delta_weights, c='r')
                axes[i][0].plot(binned_x, target_ramp, label='Experiment', c='k')
                axes[i][0].plot(binned_x, model_ramp, label='Model', c='r')
                axes[i][1].set_ylabel('Change in synaptic weight (a.u.)')
                axes[i][1].set_xlabel('Location (cm)')
                axes[i][0].set_ylabel('Ramp amplitude (mV)')
                axes[i][0].set_xlabel('Location (cm)')
                axes[i][1].set_title('Induction: %i' % int(float(induction_key)), fontsize=mpl.rcParams['font.size'])
                ymin = min(ymin, np.min(model_ramp) - 1., np.min(target_ramp) - 1.)
                ymax = max(ymax, np.max(model_ramp) + 1., np.max(target_ramp) + 1.)
            try:
                self_consistent_model_ramp = f['model_ramp_self_consistent_features'][cell_key]['2']['model_ramp'][:]
                self_consistent_initial_weights = \
                    f['model_ramp_features'][cell_key]['1']['model_weights'][:]
                self_consistent_model_weights = \
                    f['model_ramp_self_consistent_features'][cell_key]['2']['model_weights'][:]
                self_consistent_delta_weights = np.subtract(self_consistent_model_weights,
                                                            self_consistent_initial_weights)
                axes[1][1].plot(peak_locs, self_consistent_delta_weights, c='c')
                axes[1][0].plot(binned_x, self_consistent_model_ramp, label='Model (self-consistent)', c='c')
                ymin = min(ymin, np.min(self_consistent_model_ramp) - 1.)
                ymax = max(ymax, np.max(self_consistent_model_ramp) + 1.)
            except KeyError:
                pass
            for induction_key in f[description][cell_key]:
                i = int(float(induction_key)) - 1
                axes[i][0].set_ylim([ymin, ymax])
                axes[i][1].set_ylim([-peak_weight, peak_weight])
                axes[i][0].legend(loc='best', frameon=False, framealpha=0.5)
            clean_axes(axes)
            fig.suptitle('Cell_id: %i' % int(float(cell_key)), fontsize=mpl.rcParams['font.size'])
            fig.tight_layout()
            plt.show()
            plt.close()
    mpl.rcParams['font.size'] = orig_fontsize


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
        axes.plot(signal_xrange, pot_rate, label='Potentiation rate', c='k')
        axes.plot(signal_xrange, depot_rate, label='Depotentiation rate', c='r')
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
                axes[i][1].plot(peak_locs, delta_weights, c='r')
                axes[i][0].plot(binned_x, target_ramp, label='Experiment', c='k')
                axes[i][0].plot(binned_x, model_ramp, label='Model', c='r')
                axes[i][1].set_ylabel('Change in synaptic weight (a.u.)')
                axes[i][1].set_xlabel('Location (cm)')
                axes[i][0].set_ylabel('Ramp amplitude (mV)')
                axes[i][0].set_xlabel('Location (cm)')
                axes[i][1].set_title('Induction: %i' % int(float(induction_key)), fontsize=mpl.rcParams['font.size'])
                ymin = min(ymin, np.min(model_ramp) - 1., np.min(target_ramp) - 1.)
                ymax = max(ymax, np.max(model_ramp) + 1., np.max(target_ramp) + 1.)
            try:
                self_consistent_model_ramp = f['model_ramp_self_consistent_features'][cell_key]['2']['model_ramp'][:]
                self_consistent_initial_weights = \
                    f['model_ramp_self_consistent_features'][cell_key]['2']['initial_weights'][:]
                self_consistent_weights = f['model_ramp_self_consistent_features'][cell_key]['2']['weights'][:]
                self_consistent_delta_weights = np.subtract(self_consistent_weights, self_consistent_initial_weights)
                axes[1][1].plot(peak_locs, self_consistent_delta_weights, c='c')
                axes[1][0].plot(binned_x, self_consistent_model_ramp, label='Model (self-consistent)', c='c')
                ymin = min(ymin, np.min(self_consistent_model_ramp) - 1.)
                ymax = max(ymax, np.max(self_consistent_model_ramp) + 1.)
            except KeyError:
                pass
            for induction_key in f[description][cell_key]:
                i = int(float(induction_key)) - 1
                axes[i][0].set_ylim([ymin, ymax])
                axes[i][1].set_ylim([-peak_weight, peak_weight])
                axes[i][0].legend(loc='best', frameon=False, framealpha=0.5)
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
