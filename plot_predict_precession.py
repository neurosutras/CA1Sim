from plot_results import *
import random

svg_title = '072516'  #  - compare density distributions'
# svg_title = None

if svg_title is not None:
    remember_font_size = mpl.rcParams['font.size']
    mpl.rcParams['font.size'] = 8

"""
saved_parameters, ampa_forces, bin_centers, density, weights, intra_peaks, intra_phases = {}, {}, {}, {}, {}, {}, {}
filenames = {'uniform': '050916 - uniform density - predict precession parameters.pkl',
             'biased': '050916 - biased density - predict precession parameters.pkl'}
for condition in filenames:
    saved_parameters[condition] = read_from_pkl(data_dir+filenames[condition])
    stim_t, ampa_forces[condition], bin_centers[condition], density[condition], weights[condition], \
        intra_peaks[condition], intra_phases[condition] = saved_parameters[condition]

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

# rec_filename = 'output050316 - cell107 - e3200-i500-type_A_0-inh0 - 10 trials'
rec_filename = 'output060716 - cell141 - e3200-i600-subtr_inh_0-inh0 - 10 trials'
local_random = random.Random()
local_random.seed(6)
peak_times, peak_phases, binned_times, binned_phases = {}, {}, {}, {}

for parameter in peak_times, peak_phases, binned_times, binned_phases:
    for condition in ['signal', 'noise', 'inh', 'inh_all', 'exc_all']:
        parameter[condition] = []
with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
    trial = f.itervalues().next()
    time_offset = trial.attrs['phase_offset']
    for key in trial['train']:
        if trial['train'][key].attrs['group'] == 'CA3':
            peak_time_array, peak_phase_array = get_waveform_phase_vs_time(trial['successes'][key],
                                                                           time_offset=time_offset)
            peak_loc = trial['train'][key].attrs['peak_loc']
            if 4250. <= peak_loc <= 4750.:
                if local_random.uniform(0., 1.) <= 0.25:
                    condition = 'signal'
                    peak_times[condition].extend(peak_time_array)
                    peak_phases[condition].extend(peak_phase_array)
            else:
                condition = 'noise'
                if local_random.uniform(0., 1.) <= 0.2:
                    peak_times[condition].extend(peak_time_array)
                    peak_phases[condition].extend(peak_phase_array)
            condition = 'exc_all'
            peak_times[condition].extend(peak_time_array)
            peak_phases[condition].extend(peak_phase_array)
    inh_keys = {}
    for key in trial['inh_train']:
        group = trial['inh_train'][key].attrs['group']
        if group not in inh_keys:
            inh_keys[group] = []
        inh_keys[group].append(key)
        condition = 'inh_all'
        peak_time_array, peak_phase_array = get_waveform_phase_vs_time(trial['inh_train'][key],
                                                                       time_offset=time_offset)
        peak_times[condition].extend(peak_time_array)
        peak_phases[condition].extend(peak_phase_array)
    for group in inh_keys:
        for key in local_random.sample(inh_keys[group], int(0.05*len(inh_keys[group]))):
            peak_time_array, peak_phase_array = get_waveform_phase_vs_time(trial['inh_train'][key],
                                                                           time_offset=time_offset)
            condition = 'inh'
            peak_times[condition].extend(peak_time_array)
            peak_phases[condition].extend(peak_phase_array)
"""
hist, edges = {}, {}
colors = ['k', 'purple']
for i, condition in enumerate(['exc_all', 'inh_all']):
    hist[condition], edges[condition] = np.histogram(peak_phases[condition], np.arange(0., 375., 15.), density=True)
    fig, axes = plt.subplots(1)
    axes.plot(hist[condition]*15.*100., np.add(edges[condition][:-1], 7.5), color=colors[i])
    axes.set_ylim(0., 360.)
    #axes.set_yticks([0., 90., 180., 270., 360.])
    #axes.set_ylabel('Theta phase (o)')
    axes.set_xlim(2.5, 5.5)
    axes.set_xlabel('Prob. (%)')
    axes.set_xticks([3., 4., 5.])
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
for condition in ['exc_all', 'signal', 'inh_all']:
    binned_times[condition], binned_phases[condition] = plot_phase_precession([np.array(peak_times[condition])],
                                                                              [np.array(peak_phases[condition])],
                                                                              condition, plot=False, adjust=True,
                                                                              num_bins=10)
    indexes = np.where((np.array(binned_times[condition]) > 2000.) & (np.array(binned_times[condition]) < 7000.))[0]
    binned_times[condition] = np.array(binned_times[condition][indexes])
    binned_phases[condition] = np.array(binned_phases[condition][indexes])

fig, axes = plt.subplots(1)

for condition in ['signal', 'noise', 'inh']:
    indexes = np.where((np.array(peak_times[condition]) > 2000.) & (np.array(peak_times[condition]) < 7000.))[0]
    peak_times[condition] = np.array(peak_times[condition])[indexes]
    peak_phases[condition] = np.array(peak_phases[condition])[indexes]
axes.scatter(peak_times['signal'], peak_phases['signal'], s=0.1, color='r', label='Signal', zorder=2)
axes.scatter(peak_times['noise'], peak_phases['noise'], s=0.1, color='lightgrey', label='Noise', zorder=0)
axes.plot(binned_times['exc_all'], binned_phases['exc_all'], color='k', zorder=1)
axes.plot(binned_times['signal'], binned_phases['signal'], color='r', zorder=3)
axes.set_ylim(0., 360.)
axes.set_yticks([0., 90., 180., 270., 360.])
axes.set_ylabel('Theta phase (o)')
axes.set_xlim(2000., 7000.)
axes.set_xlabel('Time (s)')
axes.set_xticks([2500., 3500., 4500., 5500., 6500.])
axes.set_xticklabels([-2, -1, 'center', 1, 2])
#axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'], scatterpoints=1)
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    #fig.set_size_inches(2.05, 1.40)
    fig.set_size_inches(1.54, 1.40)
    fig.savefig(data_dir+svg_title+' - CA3 Population Input Precession.svg', format='svg', transparent=True)
plt.show()
plt.close()


fig, axes = plt.subplots(1)
axes.scatter(peak_times['inh'], peak_phases['inh'], s=0.1, color='purple', label='Inhibitory inputs')
axes.plot(binned_times['inh_all'], binned_phases['inh_all'], color='purple')
axes.set_ylim(0., 360.)
axes.set_yticks([0., 90., 180., 270., 360.])
axes.set_ylabel('Theta phase (o)')
axes.set_xlim(2000., 7000.)
axes.set_xlabel('Time (s)')
axes.set_xticks([2500., 3500., 4500., 5500., 6500.])
axes.set_xticklabels([-2, -1, 'center', 1, 2])
#axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'], scatterpoints=1)
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    #fig.set_size_inches(2.05, 1.40)
    fig.set_size_inches(1.54, 1.40)
    fig.savefig(data_dir+svg_title+' - Inhibitory Input Precession.svg', format='svg', transparent=True)
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
"""
for condition in ['CA3', 'peri-som']:
    peak_times[condition] = []
    peak_phases[condition] = []

with h5py.File(data_dir+'output060716 - cell141 - e3200-i600-subtr_inh_0-inh0 - 10 trials.hdf5', 'r') as f:
    trial = f.itervalues().next()
    for exc_key in trial['train']:
        peak_loc = trial['train'][exc_key].attrs['peak_loc']
        if 1495. <= peak_loc <= 1505.:
            single_train = trial['train'][exc_key][:]
            break
    for inh_key in trial['inh_train']:
        if trial['inh_train'][inh_key].attrs['group'] == 'perisomatic':
            break
    for trial in f.itervalues():
        time_offset = trial.attrs['phase_offset']
        peak_time_array, peak_phase_array = get_waveform_phase_vs_time(trial['train'][exc_key][:],
                                                                       time_offset=time_offset)
        peak_times['CA3'].extend(peak_time_array)
        peak_phases['CA3'].extend(peak_phase_array)
        peak_time_array, peak_phase_array = get_waveform_phase_vs_time(trial['inh_train'][inh_key][:],
                                                                       time_offset=time_offset)
        peak_times['peri-som'].extend(peak_time_array)
        peak_phases['peri-som'].extend(peak_phase_array)
fig, axes = plt.subplots(1)
axes.scatter(peak_times['CA3'], peak_phases['CA3'], label='CA3 excitatory', color='c', zorder=1, s=0.1)
axes.scatter(peak_times['peri-som'], peak_phases['peri-som'], label='Peri-som inhibitory', color='purple', zorder=0,
             s=0.1)
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
    #fig.set_size_inches(1.368, 0.907)
    fig.set_size_inches(1.3198, 1.2169)
    fig.savefig(data_dir+svg_title+' - example spike trains.svg', format='svg', transparent=True)
plt.show()
plt.close()


filename = ''
saved_parameters = read_from_pkl(data_dir+filename)
stim_t, ampa_forces, peak_locs, intra_peaks, intra_phases = saved_parameters

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

"""
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
"""

if svg_title is not None:
    mpl.rcParams['font.size'] = remember_font_size
