__author__ = 'Aaron D. Milstein'
import matplotlib as mpl
import matplotlib.lines as mlines
import scipy.stats as stats
import matplotlib.gridspec as gridspec
from matplotlib import cm
from function_lib import *

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.size'] = 12.
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['text.usetex'] = False


def plot_superimpose_conditions(rec_filename, legend=False):
    """
    File contains simulation results from iterating through some changes in parameters or stimulation conditions.
    This function produces one plot per recorded vector. Each plot superimposes the recordings from each of the
    simulation iterations.
    :param rec_filename: str
    :param legend: bool
    """
    with h5py.File(data_dir+rec_filename, 'r') as f:
        rec_ids = []
        sim_ids = []
        for sim in f.values():
            if 'description' in sim.attrs and not sim.attrs['description'] in sim_ids:
                sim_ids.append(sim.attrs['description'])
            for rec in sim['rec'].values():
                if 'description' in rec.attrs:
                    rec_id = rec.attrs['description']
                else:
                    rec_id = rec.attrs['type']+str(rec.attrs['index'])
                if not rec_id in (id['id'] for id in rec_ids):
                    rec_ids.append({'id': rec_id, 'ylabel': rec.attrs['ylabel']+' ('+rec.attrs['units']+')'})
        if len(rec_ids) == 1:
            fig, axis = plt.subplots(1)
            axes = [axis]
        else:
            fig, axes = plt.subplots(1, len(rec_ids))
        for i in range(len(rec_ids)):
            axes[i].set_xlabel('Time (ms)')
            axes[i].set_ylabel(rec_ids[i]['ylabel'])
            axes[i].set_title(rec_ids[i]['id'])
        for sim in f.values():
            if 'description' in sim.attrs:
                sim_id = sim.attrs['description']
            else:
                sim_id = ''
            tvec = sim['time']
            for rec in sim['rec'].values():
                if ('description' in rec.attrs):
                    rec_id = rec.attrs['description']
                else:
                    rec_id = rec.attrs['type']+str(rec.attrs['index'])
                i = [index for index, id in enumerate(rec_ids) if id['id'] == rec_id][0]
                axes[i].plot(tvec[:], rec[:], label=sim_id)
        if legend:
            for i in range(len(rec_ids)):
                axes[i].legend(loc='best', framealpha=0.5, frameon=False)
        fig.tight_layout()
        plt.show()


def plot_synaptic_param_distribution(cell, syn_type, param_name, scale_factor=1., param_label=None,
                                 ylabel='Peak conductance', yunits='uS', svg_title=None):
    """
    Takes a cell as input rather than a file. No simulation is required, this method just takes a fully specified cell
    and plots the relationship between distance and the specified synaptic parameter for all spines. Used while
    debugging specification of synaptic parameters.
    :param cell: :class:'HocCell'
    :param syn_type: str
    :param param_name: str
    :param scale_factor: float
    :param param_label: str
    :param ylabel: str
    :param yunits: str
    :param svg_title: str
    """
    colors = ['k', 'r', 'c', 'y', 'm', 'g', 'b']
    dend_types = ['basal', 'trunk', 'apical', 'tuft']

    if svg_title is not None:
        remember_font_size = mpl.rcParams['font.size']
        mpl.rcParams['font.size'] = 20
    fig, axes = plt.subplots(1)
    maxval, minval = None, None
    for i, sec_type in enumerate(dend_types):
        syn_list = []
        distances = []
        param_vals = []

        for branch in cell.get_nodes_of_subtype(sec_type):
            for spine in branch.spines:
                syn_list.extend(spine.synapses)
            syn_list.extend(branch.synapses)
        for syn in [syn for syn in syn_list if syn_type in syn._syn]:
            if syn.node.type == 'spine_head':
                this_distance = cell.get_distance_to_node(cell.tree.root, syn.node.parent.parent, syn.loc)
            else:
                this_distance = cell.get_distance_to_node(cell.tree.root, syn.node, syn.loc)
            distances.append(this_distance)
            if sec_type == 'basal':
                    distances[-1] *= -1
            param_vals.append(getattr(syn.target(syn_type), param_name) * scale_factor)
        if param_vals:
            axes.scatter(distances, param_vals, color=colors[i], label=sec_type)
            if maxval is None:
                maxval = max(param_vals)
            else:
                maxval = max(maxval, max(param_vals))
            if minval is None:
                minval = min(param_vals)
            else:
                minval = min(minval, min(param_vals))

    axes.set_ylabel(ylabel + ' (' + yunits + ')')
    if (maxval is not None) and (minval is not None):
        buffer = 0.1 * (maxval - minval)
        if buffer == 0:
            buffer = 0.5 * maxval
        axes.set_ylim(minval - buffer, maxval + buffer)
    axes.set_xlabel('Distance to soma (um)')
    axes.set_xlim(-200., 525.)
    axes.set_xticks([-150., 0., 150., 300., 450.])
    plt.legend(loc='best', scatterpoints=1, frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    if param_label is not None:
        plt.title(param_label, fontsize=mpl.rcParams['font.size'])
    clean_axes(axes)
    axes.tick_params(direction='out')
    fig.tight_layout()
    if not svg_title is None:
        if param_label is not None:
            svg_title = svg_title+' - '+param_label+'.svg'
        else:
            svg_title = svg_title+' - '+syn_type+'_'+param_name+' distribution.svg'
        fig.set_size_inches(5.27, 4.37)
        fig.savefig(data_dir+svg_title, format='svg', transparent=True)
    fig.show()
    plt.show()
    if svg_title is not None:
        mpl.rcParams['font.size'] = remember_font_size


def plot_mech_param_distribution(cell, mech_name, param_name, export=None, overwrite=False, scale_factor=10000.,
                        param_label=None, ylabel='Conductance density', yunits='pS/um2', svg_title=None):
    """
    Takes a cell as input rather than a file. No simulation is required, this method just takes a fully specified cell
    and plots the relationship between distance and the specified mechanism parameter for all dendritic segments. Used
    while debugging specification of mechanism parameters.
    :param cell: :class:'HocCell'
    :param mech_name: str
    :param param_name: str
    :param export: str (name of hdf5 file for export)
    :param overwrite: bool (whether to overwrite or append to potentially existing hdf5 file)
    :param scale_factor: float
    :param param_label: str
    :param ylabel: str
    :param yunits: str
    :param svg_title: str
    """
    if svg_title is not None:
        remember_font_size = mpl.rcParams['font.size']
        mpl.rcParams['font.size'] = 20
    dend_types = ['basal', 'trunk', 'apical', 'tuft']
    fig, axes = plt.subplots(1)
    maxval, minval = 1., 0.
    distances = {}
    param_vals = {}
    sec_types_list = [sec_type for sec_type in dend_types if sec_type in cell._node_dict.keys()]
    num_colors = 10
    color_x = np.linspace(0., 1., num_colors)
    colors = [cm.Set1(x) for x in color_x]
    for i, sec_type in enumerate(sec_types_list):
        distances[sec_type] = []
        param_vals[sec_type] = []
        for branch in cell.get_nodes_of_subtype(sec_type):
            for seg in [seg for seg in branch.sec if hasattr(seg, mech_name)]:
                distances[sec_type].append(cell.get_distance_to_node(cell.tree.root, branch, seg.x))
                if sec_type == 'basal':
                    distances[sec_type][-1] *= -1
                param_vals[sec_type].append(getattr(getattr(seg, mech_name), param_name) * scale_factor)
        if param_vals[sec_type]:
            axes.scatter(distances[sec_type], param_vals[sec_type],
                         color=colors[i], label=sec_type, alpha=0.5)
            if maxval is None:
                maxval = max(param_vals[sec_type])
            else:
                maxval = max(maxval, max(param_vals[sec_type]))
            if minval is None:
                minval = min(param_vals[sec_type])
            else:
                minval = min(minval, min(param_vals[sec_type]))
    axes.set_xlabel('Distance to soma (um)')
    distances_list = []
    for dist_list in distances.values():
        distances_list.extend(dist_list)
    xmax0 = max(0.1, max(distances_list))
    xmin0 = min(0, min(distances_list))
    xmin = xmin0 - 0.01 * (xmax0 - xmin0)
    xmax = xmax0 + 0.01 * (xmax0 - xmin0)
    axes.set_xlim(xmin, xmax)
    axes.set_ylabel(ylabel + ' (' + yunits + ')')
    if (maxval is not None) and (minval is not None):
        buffer = 0.01 * (maxval - minval)
        if buffer == 0:
            buffer = 0.5 * maxval
        axes.set_ylim(minval - buffer, maxval + buffer)
    if param_label is not None:
        axes.set_title(param_label, fontsize=mpl.rcParams['font.size'])
    axes.legend(loc='best', scatterpoints=1, frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    clean_axes(axes)
    axes.tick_params(direction='out')
    if not svg_title is None:
        if param_label is not None:
            svg_title = svg_title + ' - ' + param_label + '.svg'
        else:
            svg_title = svg_title + ' - ' + mech_name + '_' + param_name + ' distribution.svg'
        fig.set_size_inches(5.27, 4.37)
        fig.savefig(data_dir + svg_title, format='svg', transparent=True)
    plt.show()
    plt.close()
    if svg_title is not None:
        mpl.rcParams['font.size'] = remember_font_size

    if export is not None:
        if overwrite:
            f = h5py.File(data_dir + export + '.hdf5', 'w')
        else:
            f = h5py.File(data_dir + export + '.hdf5', 'a')
        if 'mech_filename' in f.attrs.keys():
            if not (f.attrs['mech_filename'] == '{}'.format(cell.mech_filename)):
                raise Exception('Specified mechanism filename {} does not match the mechanism filename '
                                'of the cell {}'.format(f.attrs['mech_filename'], cell.mech_filename))
        else:
            f.attrs['mech_filename'] = '{}'.format(cell.mech_filename)
        if mech_name in f.keys():
            if param_name in f[mech_name].keys():
                return
            else:
                f[mech_name].create_group(param_name)
        else:
            f.create_group(mech_name)
            f[mech_name].create_group(param_name)
        for sec_type in param_vals.keys():
            f[mech_name][param_name].create_group(sec_type)
            f[mech_name][param_name][sec_type].create_dataset('values', data=param_vals[sec_type])
        if not 'distances' in f.keys():
            f.create_group('distances')
            for sec_type in distances.keys():
                f['distances'].create_group(sec_type)
                f['distances'][sec_type].create_dataset('values', data=distances[sec_type])
        f.close()


def process_patterned_input_simulation(rec_filename, title, dt=0.02):
    """

    :param rec_file_name: str
    :param title: str
    :param dt: float
    :return: list of array
    # remember .attrs['phase_offset'] could be inside ['train'] for old files
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
        sim = next(iter(f.values()))
        equilibrate = sim.attrs['equilibrate']
        track_equilibrate = sim.attrs['track_equilibrate']
        track_length = sim.attrs['track_length']
        if 'input_field_duration' in sim.attrs:
            input_field_duration = sim.attrs['input_field_duration']
        elif 'input_field_width' in sim.attrs and 'run_vel' in sim.attrs:
            input_field_duration = sim.attrs['input_field_width'] / sim.attrs['run_vel'] * 1000.
        duration = sim.attrs['duration']
        stim_dt = sim.attrs['stim_dt']
        bins = int((1.5 + track_length) * input_field_duration / 20.)
        track_duration = duration - equilibrate - track_equilibrate
        stim_t = np.arange(-track_equilibrate, track_duration, stim_dt)
        start = int(track_equilibrate/stim_dt)
        spatial_bin = input_field_duration/50.
        intervals = []
        pop_input = []
        output = []
        for sim in f.values():
            exc_input_sum = None
            for key, train in sim['train'].items():
                this_train = np.array(train)
                intervals.extend(np.diff(this_train))
                this_exc_rate = get_binned_firing_rate(this_train, stim_t)
                if exc_input_sum is None:
                    exc_input_sum = np.array(this_exc_rate)
                else:
                    exc_input_sum = np.add(exc_input_sum, this_exc_rate)
            pop_input.append(exc_input_sum)
            this_output = get_smoothed_firing_rate(np.array(sim['output']), stim_t[start:], bin_dur=3.*spatial_bin,
                                           bin_step=spatial_bin, dt=stim_dt)
            output.append(this_output)
        pop_psd = []
        for this_pop_input in pop_input:
            pop_freq, this_pop_psd = signal.periodogram(this_pop_input, 1000./stim_dt)
            pop_psd.append(this_pop_psd)
        pop_psd = np.mean(pop_psd, axis=0)
        left = np.where(pop_freq >= 4.)[0][0]
        right = np.where(pop_freq >= 11.)[0][0]
        pop_psd /= np.max(pop_psd[left:right])
        mean_output = np.mean(output, axis=0)
        plt.hist(intervals, bins=int((max(intervals) - min(intervals)) / 3.), normed=True)
        plt.xlim(0., 200.)
        plt.ylabel('Probability')
        plt.xlabel('Inter-Spike Interval (ms)')
        plt.title('Distribution of Input Inter-Spike Intervals - '+title)
        plt.show()
        plt.close()
        peak_locs = [sim.attrs['peak_loc'] for sim in next(iter(f.values()))['train'].itervalues()]
        plt.hist(peak_locs, bins=bins)
        plt.xlabel('Time (ms)')
        plt.ylabel('Count (20 ms Bins)')
        plt.title('Distribution of Input Peak Locations - '+title)
        plt.xlim((np.min(peak_locs), np.max(peak_locs)))
        plt.show()
        plt.close()
        for sim in f.values():
            t = np.arange(0., duration, dt)
            vm = np.interp(t, sim['time'], sim['rec']['0'])
            start = int((equilibrate + track_equilibrate)/dt)
            plt.plot(np.subtract(t[start:], equilibrate + track_equilibrate), vm[start:])
            plt.xlabel('Time (ms)')
            plt.ylabel('Voltage (mV)')
            plt.title('Somatic Vm - '+title)
            plt.ylim((-70., -50.))
        plt.show()
        plt.close()
    rec_t = np.arange(0., track_duration, dt)
    #spikes_removed = get_removed_spikes_alt(rec_filename, plot=0)
    spikes_removed = get_removed_spikes(rec_filename, plot=0)
    # down_sample traces to 2 kHz after clipping spikes for theta and ramp filtering
    down_dt = 0.5
    down_t = np.arange(0., track_duration, down_dt)
    # 2000 ms Hamming window, ~2 Hz low-pass for ramp, ~5 - 10 Hz bandpass for theta, ~0.2 Hz low-pass for residuals
    window_len = int(2000. / down_dt)
    pad_len = int(window_len / 2.)
    theta_filter = signal.firwin(window_len, [5., 10.], nyq=1000. / 2. / down_dt, pass_zero=False)
    ramp_filter = signal.firwin(window_len, 2., nyq=1000. / 2. / down_dt)
    slow_vm_filter = signal.firwin(window_len, .2, nyq=1000. / 2. / down_dt)
    theta_traces = []
    theta_removed = []
    ramp_traces = []
    slow_vm_traces = []
    residuals = []
    intra_psd = []
    theta_envelopes = []
    for trace in spikes_removed:
        intra_freq, this_intra_psd = signal.periodogram(trace, 1000. / dt)
        intra_psd.append(this_intra_psd)
        down_sampled = np.interp(down_t, rec_t, trace)
        padded_trace = np.zeros(len(down_sampled) + window_len)
        padded_trace[pad_len:-pad_len] = down_sampled
        padded_trace[:pad_len] = down_sampled[::-1][-pad_len:]
        padded_trace[-pad_len:] = down_sampled[::-1][:pad_len]
        filtered = signal.filtfilt(theta_filter, [1.], padded_trace, padlen=pad_len)
        this_theta_envelope = np.abs(signal.hilbert(filtered))
        filtered = filtered[pad_len:-pad_len]
        up_sampled = np.interp(rec_t, down_t, filtered)
        theta_traces.append(up_sampled)
        this_theta_removed = trace - up_sampled
        theta_removed.append(this_theta_removed)
        this_theta_envelope = this_theta_envelope[pad_len:-pad_len]
        up_sampled = np.interp(rec_t, down_t, this_theta_envelope)
        theta_envelopes.append(up_sampled)
        filtered = signal.filtfilt(ramp_filter, [1.], padded_trace, padlen=pad_len)
        filtered = filtered[pad_len:-pad_len]
        up_sampled = np.interp(rec_t, down_t, filtered)
        ramp_traces.append(up_sampled)
        filtered = signal.filtfilt(slow_vm_filter, [1.], padded_trace, padlen=pad_len)
        filtered = filtered[pad_len:-pad_len]
        up_sampled = np.interp(rec_t, down_t, filtered)
        slow_vm_traces.append(up_sampled)
        this_residual = this_theta_removed - up_sampled
        residuals.append(this_residual)
    intra_psd = np.mean(intra_psd, axis=0)
    left = np.where(intra_freq >= 4.)[0][0]
    right = np.where(intra_freq >= 11.)[0][0]
    intra_psd /= np.max(intra_psd[left:right])
    # mean_across_trials = np.mean(theta_removed, axis=0)
    # variance_across_trials = np.var(theta_removed, axis=0)
    binned_mean = [[] for i in range(len(residuals))]
    binned_variance = [[] for i in range(len(residuals))]
    binned_t = []
    bin_duration = 3. * spatial_bin
    interval = int(bin_duration / dt)
    for j in range(0, int(track_duration / bin_duration)):
        binned_t.append(j * bin_duration + bin_duration / 2.)
        for i, residual in enumerate(residuals):
            binned_variance[i].append(np.var(residual[j * interval:(j + 1) * interval]))
            binned_mean[i].append(np.mean(theta_removed[i][j * interval:(j + 1) * interval]))
    mean_theta_envelope = np.mean(theta_envelopes, axis=0)
    mean_ramp = np.mean(ramp_traces, axis=0)
    mean_binned_vm = np.mean(binned_mean, axis=0)
    mean_binned_var = np.mean(binned_variance, axis=0)
    scatter_vm_mean = np.array(binned_mean).flatten()
    scatter_vm_var = np.array(binned_variance).flatten()
    print('Mean Theta Envelope for %s: %.2f' % (title, np.mean(mean_theta_envelope)))
    plt.plot(binned_t, mean_binned_vm)
    plt.xlabel('Time - 180 ms bins')
    plt.ylabel('Voltage (mV)')
    plt.title('Somatic Vm Mean - Across Trials - ' + title)
    plt.show()
    plt.close()
    plt.plot(binned_t, mean_binned_var)
    plt.xlabel('Time (ms)')
    plt.ylabel('Vm Variance (mV' + r'$^2$' + ')')
    plt.title('Somatic Vm Variance - Across Trials - ' + title)
    plt.show()
    plt.close()
    plt.scatter(scatter_vm_mean, scatter_vm_var)
    plt.xlabel('Mean Vm (mV)')
    plt.ylabel('Vm Variance (mV' + r'$^2$' + ')')
    plt.title('Mean - Variance Analysis - ' + title)
    plt.show()
    plt.close()
    plt.plot(pop_freq, pop_psd, label='Total Population Input Spikes')
    plt.plot(intra_freq, intra_psd, label='Single Cell Intracellular Vm')
    plt.xlim(4., 11.)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Normalized Power Density')
    plt.title('Power Spectral Density - ' + title)
    plt.legend(loc='best')
    plt.show()
    plt.close()
    return rec_t, residuals, mean_theta_envelope, scatter_vm_mean, scatter_vm_var, binned_t, mean_binned_vm, \
           mean_binned_var, mean_ramp, mean_output


def process_patterned_input_simulation_theta_freq(rec_filenames, conditions=None, theta_dur=None, field_center=4500.,
                                                  win_dur=1800., dt=0.02):
    """
    :param rec_file_names: dict of str
    :param conditions: list of str
    :param theta_dur: dict of {str: float}
    :param field_center: float
    :param win_dur: float
    :param dt: float
    :return: list of dict of array
    """
    if conditions is None:
        conditions = ['modinh0', 'modinh3']
    if theta_dur is None:
        theta_dur = {'orig': 150., 'modinh': 145.}
    peaks, phases, IPI, binned_peaks, binned_phases, binned_t, binned_IPI, theta_env = {}, {}, {}, {}, {}, {}, \
                                                                                            {}, {}
    for parameter in peaks, phases, IPI, binned_peaks, binned_phases, binned_t, binned_IPI:
        for group in ['exc', 'successes', 'inh', 'intra']:
            parameter[group] = {}
    for condition in conditions:
        rec_filename = rec_filenames[condition]
        with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
            sim = next(iter(f.values()))
            equilibrate = sim.attrs['equilibrate']
            track_equilibrate = sim.attrs['track_equilibrate']
            duration = sim.attrs['duration']
            stim_dt = sim.attrs['stim_dt']
            track_duration = duration - equilibrate - track_equilibrate
            stim_t = np.arange(-track_equilibrate, track_duration, stim_dt)
            stim_t_start = int(track_equilibrate/stim_dt)
            t = stim_t[stim_t_start:] - stim_t[stim_t_start]
            rec_t = np.arange(0., track_duration, dt)
            down_dt = 0.5
            # ~5 - 10 Hz bandpass for theta
            window_len = int(2000. / down_dt)
            theta_filter = signal.firwin(window_len, [5., 10.], nyq=1000. / 2. / down_dt, pass_zero=False)
            time_offset = {'orig': [trial.attrs['phase_offset'] for trial in f.values()]}
            if 'mod_inh_time_offset' in sim.attrs:
                time_offset['modinh'] = [trial.attrs['mod_inh_time_offset'] for trial in f.values()]
            for i, trial in enumerate(f.values()):
                for group, key in zip(['exc', 'successes', 'inh'], ['train', 'successes', 'inh_train']):
                    if key in trial:
                        input_sum = None
                        for train in trial[key].values():
                            this_rate = get_binned_firing_rate(train[:], stim_t)
                            if input_sum is None:
                                input_sum = np.array(this_rate)
                            else:
                                input_sum = np.add(input_sum, this_rate)
                        if condition not in peaks[group]:
                            peaks[group][condition] = []
                            phases[group][condition] = {}
                            IPI[group][condition] = []
                        theta_trace = general_filter_trace(t, input_sum[stim_t_start:], filter=theta_filter,
                                                          duration=track_duration, dt=stim_dt)
                        for LFP_type in time_offset:
                            this_peaks, this_phases = get_waveform_phase_vs_time(t, theta_trace,
                                                                                 cycle_duration=theta_dur[LFP_type],
                                                                                 time_offset=time_offset[LFP_type][i])
                            if LFP_type not in phases[group][condition]:
                                phases[group][condition][LFP_type] = []
                            if LFP_type == 'orig':
                                peaks[group][condition].append(this_peaks)
                                this_IPI = np.diff(this_peaks)
                                IPI[group][condition].append(this_IPI)
                            phases[group][condition][LFP_type].append(this_phases)
        spikes_removed = get_removed_spikes(rec_filename, dt=dt, plot=0)
        group = 'intra'
        for i, trace in enumerate(spikes_removed):
            if condition not in peaks[group]:
                peaks[group][condition] = []
                phases[group][condition] = {}
                IPI[group][condition] = []
                theta_env[condition] = []
            theta_trace = general_filter_trace(rec_t, trace, filter=theta_filter,
                                               duration=track_duration, dt=dt)
            this_theta_env = np.abs(signal.hilbert(theta_trace))
            theta_env[condition].append(this_theta_env)
            for LFP_type in time_offset:
                this_peaks, this_phases = get_waveform_phase_vs_time(rec_t, theta_trace,
                                                                     cycle_duration=theta_dur[LFP_type],
                                                                     time_offset=time_offset[LFP_type][i])
                if LFP_type not in phases[group][condition]:
                    phases[group][condition][LFP_type] = []
                if LFP_type == 'orig':
                    peaks[group][condition].append(this_peaks)
                    this_IPI = np.diff(this_peaks)
                    IPI[group][condition].append(this_IPI)
                phases[group][condition][LFP_type].append(this_phases)
    start = field_center - win_dur / 2.
    end = start + win_dur
    for group in peaks:
        for condition in peaks[group]:
            binned_phases[group][condition] = {}
            for LFP_type in phases[group][condition]:
                binned_peaks[group][condition], binned_phases[group][condition][LFP_type] = \
                    plot_phase_precession(peaks[group][condition], phases[group][condition][LFP_type],
                                          group+'_'+condition+'; LFP: '+LFP_type, fit_start=start,
                                          fit_end=end)
            binned_t[group][condition], binned_IPI[group][condition] = plot_IPI(peaks[group][condition],
                                                                                IPI[group][condition],
                                                                                group+'_'+condition)
    return t, rec_t, peaks, phases, IPI, binned_peaks, binned_phases, binned_t, binned_IPI, theta_env


def plot_patterned_input_sim_summary(rec_t, mean_theta_envelope, binned_t,  mean_binned_var, mean_ramp, mean_output,
                                     key_list=None, titles=None, baseline_range=[0., 600.], dt=0.02, svg_title=None):
    """
    Expects the output of process_patterned_input_simulation.
    Produces summary plots for ramp, variance, theta, and firing rate.
    :param rec_t: array
    :param mean_theta_envelope: array
    :param binned_t: array
    :param mean_binned_var: array
    :param mean_ramp: array
    :param mean_output: array
    :param key_list: list of str
    :param titles: list of str
    :param baseline_range: list of float
    :param dt: float
    :param svg_title: str
    """
    if svg_title is not None:
        remember_font_size = mpl.rcParams['font.size']
        mpl.rcParams['font.size'] = 20
    if key_list is None:
        key_list = ['modinh0', 'modinh1', 'modinh2']
    if titles is None:
        titles = ['Control', 'Reduced inhibition - In field', 'Reduced inhibition - Out of field']
    colors = ['k', 'y', 'orange']
    fig, axes = plt.subplots(1)
    baseline = np.mean(mean_ramp[key_list[0]][int(baseline_range[0]/dt):int(baseline_range[1]/dt)])
    for i, (condition, title) in enumerate(zip([key_list[0], key_list[2], key_list[1]], titles)):
        axes.plot(rec_t, np.subtract(mean_ramp[condition], baseline), color=colors[i], label=title)
    clean_axes(axes)
    axes.set_xlabel('Time (s)')
    axes.set_ylabel('DVm (mV)')
    axes.set_ylim(-0.8, 9.)
    axes.set_xlim(0., 7500.)
    axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
    axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
    axes.tick_params(direction='out')
    plt.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    if svg_title is not None:
        fig.set_size_inches(4.403, 3.631)
        fig.savefig(data_dir+svg_title+' - Summary - Ramp.svg', format='svg', transparent=True)
    plt.show()
    plt.close()
    fig, axes = plt.subplots(1)
    for i, (condition, title) in enumerate(zip([key_list[0], key_list[2], key_list[1]], titles)):
        axes.plot(rec_t, mean_output[condition], color=colors[i], label=title)
    clean_axes(axes)
    axes.set_xlabel('Time (s)')
    axes.set_ylabel('Firing rate (Hz)')
    axes.set_ylim(0., 45.)
    axes.set_xlim(0., 7500.)
    axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
    axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
    # plt.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    axes.tick_params(direction='out')
    if svg_title is not None:
        fig.set_size_inches(4.403, 3.631)
        fig.savefig(data_dir + svg_title + ' - Summary - Rate.svg', format='svg', transparent=True)
    plt.show()
    plt.close()
    fig, axes = plt.subplots(1)
    for i, (condition, title) in enumerate(zip([key_list[0], key_list[2], key_list[1]], titles)):
        axes.plot(rec_t, mean_theta_envelope[condition], color=colors[i], label=title)
    clean_axes(axes)
    axes.set_xlabel('Time (s)')
    axes.set_ylabel('Thetaintra (mV)')
    axes.set_ylim(0., 2.5)
    axes.set_xlim(0., 7500.)
    axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
    axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
    # plt.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    axes.tick_params(direction='out')
    if svg_title is not None:
        fig.set_size_inches(4.403, 3.631)
        fig.savefig(data_dir + svg_title + ' - Summary - Theta.svg', format='svg', transparent=True)
    plt.show()
    plt.close()
    mpl.rcParams['font.size'] = 8
    fig, axes = plt.subplots(1)
    for i, (condition, title) in enumerate(zip([key_list[0], key_list[2], key_list[1]], titles)):
        axes.plot(binned_t[condition], mean_binned_var[condition], color=colors[i], label=title)
    clean_axes(axes)
    axes.set_xlabel('Time (s)', fontsize=8)
    axes.set_ylabel('Variance (mV2)', fontsize=8)
    axes.set_ylim(0., 7.)
    axes.set_xlim(0., 7500.)
    axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
    axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
    plt.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    axes.tick_params(direction='out')
    if svg_title is not None:
        fig.set_size_inches(1.95, 1.16)
        fig.savefig(data_dir + svg_title + ' - Summary - Variance.svg', format='svg', transparent=True)
    plt.show()
    plt.close()
    if svg_title is not None:
        mpl.rcParams['font.size'] = remember_font_size
