from plot_results import *
import random

svg_title = '050916 - compare density distributions'
#svg_title = None

saved_parameters, ampa_forces, bin_centers, density, weights, intra_peaks, intra_phases = {}, {}, {}, {}, {}, {}, {}
filenames = {'uniform': '050916 - uniform density - predict precession parameters.pkl',
             'biased': '050916 - biased density - predict precession parameters.pkl'}
for condition in filenames:
    saved_parameters[condition] = read_from_pkl(data_dir+filenames[condition])
    stim_t, ampa_forces[condition], bin_centers[condition], density[condition], weights[condition], \
        intra_peaks[condition], intra_phases[condition] = saved_parameters[condition]

if svg_title is not None:
    remember_font_size = mpl.rcParams['font.size']
    mpl.rcParams['font.size'] = 8

"""
colors = ['r', 'k']
fig, axes = plt.subplots(1)
for i, (condition, title) in enumerate(zip(['biased', 'uniform'], ['Biased input density', 'Uniform input density'])):
    axes.plot(bin_centers[condition], density[condition], label=title, color=colors[i])
axes.set_ylim(0., 700.)
axes.set_ylabel('CA3 Input Density (/s)')
axes.set_xlim(0., 7500.)
axes.set_xlabel('Time (s)')
axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    fig.set_size_inches(2.05, 1.40)
    fig.savefig(data_dir+svg_title+' - Input Density.svg', format='svg', transparent=True)
plt.show()
plt.close()

fig, axes = plt.subplots(1)
for i, (condition, title) in enumerate(zip(['biased', 'uniform'], ['Biased input density', 'Uniform input density'])):
    axes.plot(intra_peaks[condition], intra_phases[condition], label=title, color=colors[i])
axes.set_ylim(0., 360.)
axes.set_yticks([0., 90., 180., 270., 360.])
axes.set_ylabel('Theta phase (o)')
axes.set_xlim(0., 7500.)
axes.set_xlabel('Time (s)')
axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
#axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    fig.set_size_inches(2.05, 1.40)
    fig.savefig(data_dir+svg_title+' - Predicted Precession.svg', format='svg', transparent=True)
plt.show()
plt.close()
"""


local_random = random.Random()
peak_times, peak_phases, binned_times, binned_phases = {}, {}, {}, {}
"""
for parameter in peak_times, peak_phases, binned_times, binned_phases:
    parameter['signal'] = []
    parameter['signal_all'] = []
    parameter['noise'] = []
    #parameter['inh'] = []
    parameter['inh_all'] = []
    parameter['exc'] = []
with h5py.File(data_dir+'output050316 - cell107 - e3200-i500-type_A_0-inh0 - 10 trials.hdf5', 'r') as f:
    trial = f.itervalues().next()
    time_offset = trial.attrs['phase_offset']
    for key in trial['train']:
        if trial['train'][key].attrs['group'] == 'CA3':
            peak_time_array, peak_phase_array = get_waveform_phase_vs_time(trial['successes'][key],
                                                                           time_offset=time_offset)
            peak_loc = trial['train'][key].attrs['peak_loc']
            if 4250. <= peak_loc <= 4750.:
                condition = 'signal'
                peak_times['signal_all'].extend(peak_time_array)
                peak_phases['signal_all'].extend(peak_phase_array)
            else:
                condition = 'noise'
            if local_random.uniform(0., 1.) <= 0.25:
                peak_times[condition].extend(peak_time_array)
                peak_phases[condition].extend(peak_phase_array)
            peak_times['exc'].extend(peak_time_array)
            peak_phases['exc'].extend(peak_phase_array)
    for key in trial['inh_train']:
        group = trial['inh_train'][key].attrs['group']
        #group = 'inh'
        #if group in ['perisomatic', 'apical dendritic', 'distal apical dendritic']:
        peak_time_array, peak_phase_array = get_waveform_phase_vs_time(trial['inh_train'][key],
                                                                       time_offset=time_offset)
        if local_random.uniform(0., 1.) <= 0.05:
            if group not in peak_times:
                for parameter in peak_times, peak_phases, binned_times, binned_phases:
                    parameter[group] = []
            peak_times[group].extend(peak_time_array)
            peak_phases[group].extend(peak_phase_array)
        peak_times['inh_all'].extend(peak_time_array)
        peak_phases['inh_all'].extend(peak_phase_array)

for group in ['exc', 'signal_all', 'inh_all']:
    binned_times[group], binned_phases[group] = plot_phase_precession([np.array(peak_times[group])],
                                        [np.array(peak_phases[group])], group, plot=False, adjust=False, num_bins=10)

fig, axes = plt.subplots(1)
axes.scatter(peak_times['signal'], peak_phases['signal'], s=0.1, color='r', label='Signal', zorder=2)
axes.scatter(peak_times['noise'], peak_phases['noise'], s=0.1, color='lightgrey', label='Noise', zorder=0)
axes.plot(binned_times['exc'], binned_phases['exc'], color='k', zorder=1)
axes.plot(binned_times['signal_all'], binned_phases['signal_all'], color='r', zorder=3)
axes.set_ylim(0., 360.)
axes.set_yticks([0., 90., 180., 270., 360.])
axes.set_ylabel('Theta phase (o)')
axes.set_xlim(0., 7500.)
axes.set_xlabel('Time (s)')
axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'], scatterpoints=1)
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    #fig.set_size_inches(2.05, 1.40)
    fig.set_size_inches(1.54, 1.40)
    fig.savefig(data_dir+svg_title+' - CA3 Population Input Precession.svg', format='svg', transparent=True)
plt.show()
plt.close()

fig, axes = plt.subplots(1)
#axes.scatter(peak_times['inh'], peak_phases['inh'], s=0.01, color='purple', label='Inhibitory inputs', alpha=0.4)
for group in [group for group in ['perisomatic', 'apical dendritic', 'distal apical dendritic', 'axo-axonic',
                                  'tuft feedforward', 'tuft feedback'] if group in peak_times]:
#group = 'inh'
    axes.scatter(peak_times[group], peak_phases[group], s=0.1, color='purple', label='Inhibitory inputs', alpha=0.2)
axes.plot(binned_times['inh_all'], binned_phases['inh_all'], color='purple')
axes.set_ylim(0., 360.)
axes.set_yticks([0., 90., 180., 270., 360.])
axes.set_ylabel('Theta phase (o)')
axes.set_xlim(0., 7500.)
axes.set_xlabel('Time (s)')
axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'], scatterpoints=1)
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    #fig.set_size_inches(2.05, 1.40)
    fig.set_size_inches(1.54, 1.40)
    fig.savefig(data_dir+svg_title+' - Inhibitory Input Precession.svg', format='svg', transparent=True)
plt.show()
plt.close()

hist, edges = {}, {}
colors = ['k', 'purple']
for i, condition in enumerate(['exc', 'inh_all']):
    hist[condition], edges[condition] = np.histogram(peak_phases[condition], 24, density=True)
    fig, axes = plt.subplots(1)
    axes.plot(hist[condition]*15.*100., edges[condition][1:], color=colors[i])
    axes.set_ylim(0., 360.)
    #axes.set_yticks([0., 90., 180., 270., 360.])
    #axes.set_ylabel('Theta phase (o)')
    axes.set_xlim(3., 6.)
    axes.set_xlabel('Prob. (%)')
    axes.set_xticks([3., 4.5, 6.])
    axes.set_xticklabels([3, 4.5, 6])
    #axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'], scatterpoints=1)
    clean_axes(axes)
    axes.spines['left'].set_visible(False)
    axes.tick_params(direction='out')
    if svg_title is not None:
        # fig.set_size_inches(2.05, 1.40)
        fig.set_size_inches(0.41, 1.40)
        fig.savefig(data_dir + svg_title + ' - ' + condition + ' - Phase Histogram.svg', format='svg',
                    transparent=True)
    plt.show()
    plt.close()
"""
if svg_title is not None:
    mpl.rcParams['font.size'] = 20

for condition in ['train', 'successes']:
    peak_times[condition] = []
    peak_phases[condition] = []

with h5py.File(data_dir+'output050316 - cell107 - e3200-i500-type_A_0-inh0 - 10 trials.hdf5', 'r') as f:
    trial = f.itervalues().next()
    for key in trial['train']:
        peak_loc = trial['train'][key].attrs['peak_loc']
        if 1495. <= peak_loc <= 1505.:
            single_train = trial['train'][key][:]
            break
    """
    for trial in f.itervalues():
        time_offset = trial.attrs['phase_offset']
        for condition in ['train', 'successes']:
            peak_time_array, peak_phase_array = get_waveform_phase_vs_time(trial[condition][key][:],
                                                                           time_offset=time_offset)
            peak_times[condition].extend(peak_time_array)
            peak_phases[condition].extend(peak_phase_array)
fig, axes = plt.subplots(1)
axes.scatter(peak_times['train'], peak_phases['train'], label='Failures', color='grey', zorder=0)
axes.scatter(peak_times['successes'], peak_phases['successes'], label='Successes', color='r', zorder=1)
axes.set_ylim(0., 360.)
axes.set_yticks([0., 90., 180., 270., 360.])
axes.set_ylabel('Theta phase (o)')
axes.set_xlim(0., 3000.)
axes.set_xlabel('Time (s)')
axes.set_xticks([0., 500., 1000., 1500., 2000., 2500., 3000.])
axes.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5, 3])
axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'], scatterpoints=1)
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    #fig.set_size_inches(2.05, 1.40)
    fig.set_size_inches(4.42, 2.98)
    fig.savefig(data_dir+svg_title+' - CA3 Spike Failures.svg', format='svg', transparent=True)
plt.show()
plt.close()
"""

class Pr(object):
    """
    This object contains internal variables to track the evolution in time of parameters governing synaptic release
    probability, used during optimization, and then exported to pr.mod for use during patterned input simulations.
    """
    def __init__(self, P0=0.2, f=1.769, tau_F=67.351, d=0.878, tau_D=92.918, dt=0.01):
        self.P0 = P0
        self.f = f
        self.tau_F = tau_F
        self.d = d
        self.tau_D = tau_D
        self.P = P0
        self.dt = dt
        self.tlast = None
        self.F0 = 1.
        self.D0 = 1.
        self.F = None
        self.D = None
        self.t_array = None
        self.P_array = None

    def stim(self, stim_time):
        """
        Evolve the dynamics until the current stim_time, report the current P, and update the internal parameters.
        :param stim_time: float
        :return: float
        """
        if self.tlast is not None:
            t_segment = np.arange(self.t_array[-1], stim_time, self.dt)
            self.F = 1. + (self.F0 - 1.) * np.exp(-(t_segment - self.t_array[-1]) / self.tau_F)
            self.D = 1. - (1. - self.D0) * np.exp(-(t_segment - self.t_array[-1]) / self.tau_D)
            P_segment = np.minimum(np.ones_like(t_segment), np.multiply(self.F, self.D) * self.P0)
            self.P_array = np.append(self.P_array, P_segment)
            self.t_array = np.append(self.t_array, t_segment)
            self.F0 = self.F[-1] + self.f
            self.D0 = self.D[-1] * self.d
            self.P = self.P_array[-1]
        else:
            self.t_array = np.arange(0., stim_time, self.dt)
            self.P_array = np.ones_like(self.t_array) * self.P0
            self.F0 += self.f
            self.D0 *= self.d
        self.tlast = stim_time

thisPr = Pr()
for spike in single_train:
    thisPr.stim(spike)
fig, axes = plt.subplots(1)
axes.plot(thisPr.t_array, thisPr.P_array, color='k')
axes.set_ylim(0., 1.)
axes.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.])
axes.set_ylabel('Release Probability')
axes.set_xlim(0., 3000.)
axes.set_xlabel('Time (s)')
axes.set_xticks([0., 500., 1000., 1500., 2000., 2500., 3000.])
axes.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5, 3])
#axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    #fig.set_size_inches(2.05, 1.40)
    fig.set_size_inches(4.42, 2.98)
    fig.savefig(data_dir+svg_title+' - CA3 Example Pr.svg', format='svg', transparent=True)
plt.show()
plt.close()

if svg_title is not None:
    mpl.rcParams['font.size'] = remember_font_size
