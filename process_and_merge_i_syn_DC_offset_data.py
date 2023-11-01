import numpy as np
import sys, os
import h5py
import click
import scipy.signal as signal
from nested.utils import read_from_yaml
from function_lib import get_binned_firing_rate, get_smoothed_firing_rate


def get_removed_spikes(source_file_path, before=1.6, after=6., dt=0.02, th=10.):
    """

    :param source_file_path: str
    :param before: float : time to remove before spike
    :param after: float : time to remove after spike in case where trough or voltage recovery cannot be used
    :param dt: float : temporal resolution for interpolation and dvdt
    :param th: float : slope threshold
    :return: list of :class:'numpy.array'
    """
    removed = []
    with h5py.File(source_file_path, 'r') as f:
        sim = next(iter(f.values()))
        equilibrate = sim.attrs['equilibrate']
        duration = sim.attrs['duration']
        track_equilibrate = sim.attrs['track_equilibrate']
        for rec in f.values():
            t = np.arange(0., duration, dt)
            vm = np.interp(t, rec['time'], rec['rec']['0'])
            start = int((equilibrate + track_equilibrate) / dt)
            t = np.subtract(t[start:], equilibrate + track_equilibrate)
            vm = vm[start:]
            dvdt = np.gradient(vm, dt)
            crossings = np.where(dvdt >= th)[0]
            if not np.any(crossings):
                removed.append(vm)
            else:
                start = 0.
                i = 0
                left = max(0, crossings[i] - int(before/dt))
                previous_vm = vm[left]
                while i < len(crossings) and start < len(vm):
                    start = crossings[i]
                    left = max(0, crossings[i] - int(before/dt))
                    if not np.isnan(vm[left]):
                        previous_vm = vm[left]
                    recovers = np.where(vm[crossings[i]:] < previous_vm)[0]
                    if np.any(recovers):
                        recovers = crossings[i] + recovers[0]
                    falling = np.where(dvdt[crossings[i]:] < 0.)[0]
                    if np.any(falling):
                        falling = crossings[i] + falling[0]
                        rising = np.where(dvdt[falling:] >= 0.)[0]
                        if np.any(rising):
                            rising = falling + rising[0]
                    else:
                        rising = []
                    if np.any(recovers):
                        if np.any(rising):
                            right = min(recovers, rising)
                        else:
                            right = recovers
                    elif np.any(rising):
                        right = rising
                    else:
                        right = min(crossings[i] + int(after/dt), len(vm)-1)
                    # added to remove majority of complex spike:
                    if vm[right] >= -45. and np.any(recovers):
                        right = recovers
                    for j in range(left, right):
                        vm[j] = np.nan
                    i += 1
                    while i < len(crossings) and crossings[i] < right:
                        i += 1
                not_blank = np.where(~np.isnan(vm))[0]
                vm = np.interp(t, t[not_blank], vm[not_blank])
                removed.append(vm)
    return removed


def get_patterned_input_component_traces(source_file_path):
    """

    :param rec_file_name: str
    # remember .attrs['phase_offset'] could be inside ['train'] for old files
    """
    with h5py.File(source_file_path, 'r') as f:
        sim = next(iter(f.values()))
        equilibrate = sim.attrs['equilibrate']
        track_equilibrate = sim.attrs['track_equilibrate']
        duration = sim.attrs['duration']
        if 'dt' in sim.attrs:
            dt = sim.attrs['dt']
        else:
            dt = sim['time'][1] - sim['time'][0]
        track_duration = duration - equilibrate - track_equilibrate
        start = int((equilibrate + track_equilibrate)/dt)
        vm_array = []
        for sim in f.values():
            t = np.arange(0., duration, dt)
            vm = np.interp(t, sim['time'], sim['rec']['0'])
            vm = vm[start:]
            vm_array.append(vm)
    rec_t = np.arange(0., track_duration, dt)
    spikes_removed = get_removed_spikes(source_file_path, dt=dt)
    # down_sample traces to 2 kHz after clipping spikes for theta and ramp filtering
    down_dt = 0.5
    down_t = np.arange(0., track_duration, down_dt)
    # 2000 ms Hamming window, ~2 Hz low-pass for ramp, ~5 - 10 Hz bandpass for theta
    window_len = int(2000./down_dt)
    pad_len = int(window_len/2.)
    theta_filter = signal.firwin(window_len, [5., 10.], fs=1000./down_dt, pass_zero=False)
    ramp_filter = signal.firwin(window_len, 2., fs=1000./down_dt)
    theta_traces = []
    ramp_traces = []
    for trace in spikes_removed:
        down_sampled = np.interp(down_t, rec_t, trace)
        padded_trace = np.zeros(len(down_sampled)+window_len)
        padded_trace[pad_len:-pad_len] = down_sampled
        padded_trace[:pad_len] = down_sampled[::-1][-pad_len:]
        padded_trace[-pad_len:] = down_sampled[::-1][:pad_len]
        filtered = signal.filtfilt(theta_filter, [1.], padded_trace, padlen=pad_len)
        filtered = filtered[pad_len:-pad_len]
        up_sampled = np.interp(rec_t, down_t, filtered)
        theta_traces.append(up_sampled)
        filtered = signal.filtfilt(ramp_filter, [1.], padded_trace, padlen=pad_len)
        filtered = filtered[pad_len:-pad_len]
        up_sampled = np.interp(rec_t, down_t, filtered)
        ramp_traces.append(up_sampled)
    return rec_t, vm_array, theta_traces, ramp_traces


def get_patterned_input_filtered_synaptic_currents(source_file_path, syn_types=['AMPA', 'NMDA', 'GABA']):
    """

    :param source_file_path: str
    :param syn_types: list of str
    """
    with h5py.File(source_file_path, 'r') as f:
        sim = next(iter(f.values()))
        equilibrate = sim.attrs['equilibrate']
        track_equilibrate = sim.attrs['track_equilibrate']
        duration = sim.attrs['duration']
        track_duration = duration - equilibrate - track_equilibrate
        if 'dt' in sim.attrs:
            dt = sim.attrs['dt']
        else:
            dt = sim['time'][1] - sim['time'][0]
        start = int((equilibrate + track_equilibrate) / dt)
        t = np.arange(0., duration, dt)
        i_syn_list_dict = {}
        for syn_type in syn_types:
            i_syn_list_dict[syn_type] = []
        for sim in f.values():
            i_syn_dict = {}
            for syn_type in syn_types:
                i_syn_dict[syn_type] = np.zeros_like(t)
            for rec in sim['rec'].values():
                for syn_type in syn_types:
                    if 'i_%s' % syn_type in rec.attrs['description']:
                        interp_rec = np.interp(t, sim['time'], rec[:])
                        i_syn_dict[syn_type] += interp_rec
            for syn_type in syn_types:
                i_syn_list_dict[syn_type].append(i_syn_dict[syn_type][start:])
    
    rec_t = np.arange(0., track_duration, dt)
    # down_sample traces to 2 kHz for low-pass filtering
    down_dt = 0.5
    down_t = np.arange(0., track_duration, down_dt)
    # 2000 ms Hamming window, ~2 Hz low-pass for ramp, ~5 - 10 Hz bandpass for theta
    window_len = int(2000. / down_dt)
    pad_len = int(window_len / 2.)
    ramp_filter = signal.firwin(window_len, 2., fs=1000. / down_dt)
    filtered_i_syn_list_dict = {}
    for syn_type in syn_types:
        filtered_i_syn_list_dict[syn_type] = []
        for trace in i_syn_list_dict[syn_type]:
            down_sampled = np.interp(down_t, rec_t, trace)
            padded_trace = np.zeros(len(down_sampled) + window_len)
            padded_trace[pad_len:-pad_len] = down_sampled
            padded_trace[:pad_len] = down_sampled[::-1][-pad_len:]
            padded_trace[-pad_len:] = down_sampled[::-1][:pad_len]
            filtered = signal.filtfilt(ramp_filter, [1.], padded_trace, padlen=pad_len)
            filtered = filtered[pad_len:-pad_len]
            up_sampled = np.interp(rec_t, down_t, filtered)
            filtered_i_syn_list_dict[syn_type].append(up_sampled)
    return rec_t, filtered_i_syn_list_dict


def process_patterned_input_simulation_input_output(source_file_path):
    """

    :param source_file_path: str
    """
    with h5py.File(source_file_path, 'r') as f:
        sim = next(iter(f.values()))
        equilibrate = sim.attrs['equilibrate']
        track_equilibrate = sim.attrs['track_equilibrate']
        input_field_duration = sim.attrs['input_field_duration']
        duration = sim.attrs['duration']
        stim_dt = sim.attrs['stim_dt']
        track_duration = duration - equilibrate - track_equilibrate
        stim_t = np.arange(-track_equilibrate, track_duration, stim_dt)
        start = int(track_equilibrate/stim_dt)
        spatial_bin = input_field_duration/50.
        intervals = []
        pop_input = []
        successes = []
        inh_input = []
        output = []
        if 'successes' in sim:
            stochastic = True
        else:
            stochastic = False
        for sim in f.values():
            exc_input_sum = None
            successes_sum = None
            inh_input_sum = None
            for key, train in sim['train'].items():
                this_train = np.array(train)
                if len(this_train) > 0:
                    for i in range(len(this_train) - 1):
                        intervals.append(this_train[i+1] - this_train[i])
                this_exc_rate = get_binned_firing_rate(this_train, stim_t)
                if exc_input_sum is None:
                    exc_input_sum = np.array(this_exc_rate)
                else:
                    exc_input_sum = np.add(exc_input_sum, this_exc_rate)
                if stochastic:
                    this_success_rate = get_binned_firing_rate(np.array(sim['successes'][key]), stim_t)
                    if successes_sum is None:
                        successes_sum = np.array(this_success_rate)
                    else:
                        successes_sum = np.add(successes_sum, this_success_rate)
            pop_input.append(exc_input_sum)
            if stochastic:
                successes.append(successes_sum)
            for train in sim['inh_train'].values():
                this_inh_rate = get_binned_firing_rate(np.array(train), stim_t)
                if inh_input_sum is None:
                    inh_input_sum = np.array(this_inh_rate)
                else:
                    inh_input_sum = np.add(inh_input_sum, this_inh_rate)
            inh_input.append(inh_input_sum)
            this_output = get_smoothed_firing_rate(np.array(sim['output']), stim_t[start:], bin_dur=3.*spatial_bin,
                                           bin_step=spatial_bin, dt=stim_dt)
            output.append(this_output)
    
    pop_input_list = [pop_input[i][start:] for i in range(len(pop_input))]
    inh_input_list = [inh_input[i][start:] for i in range(len(pop_input))]
    successes_list = [successes[i][start:] for i in range(len(pop_input))]
    if stochastic:
        return stim_t[start:], pop_input_list, successes_list, inh_input_list, output
    else:
        return stim_t[start:], pop_input_list, inh_input_list, output


@click.command()
@click.option('--source-data-config-file-path', type=click.Path(exists=True, file_okay=True, dir_okay=False),
              default=None)
@click.option('--source-data-dir', type=click.Path(exists=True, file_okay=False, dir_okay=True),
              default=None)
@click.option('--target-file-path', type=str, default=None)
def main(source_data_config_file_path, source_data_dir, target_file_path):
    """
    
    :param source_data_config_file_path:
    :param source_data_dir:
    :param target_file_path:
    """
    target_data = dict()
    source_data_file_name_dict = read_from_yaml(file_path=source_data_config_file_path)
    with h5py.File(target_file_path, 'a') as target_file:
        for model_name, source_file_name_list in source_data_file_name_dict.items():
            if model_name not in target_file:
                target_file.create_group(model_name)
            model_group = target_file[model_name]
            for source_file_name in source_file_name_list:
                source_file_path = source_data_dir + '/' + source_file_name
                if not os.path.isfile(source_file_path):
                    raise Exception('Invalid file path: %s' % source_file_path)
                with h5py.File(source_file_path, 'r') as source_file:
                    seed_key = next(iter(source_file.keys()))
                    if seed_key not in model_group:
                        model_group.create_group(seed_key)
                    seed_group = model_group[seed_key]
                    data_group = source_file[seed_key]
                    seed_group.attrs.update(data_group.attrs.items())
                    seed_group.create_dataset('spike_times', data=data_group['output'][:], compression='gzip')
                
                rec_t, vm_array, theta_traces, ramp_traces = get_patterned_input_component_traces(source_file_path)
                seed_group.create_dataset('rec_t', data=rec_t, compression='gzip')
                seed_group.create_dataset('vm', data=vm_array[0], compression='gzip')
                seed_group.create_dataset('theta', data=theta_traces[0], compression='gzip')
                seed_group.create_dataset('ramp', data=ramp_traces[0], compression='gzip')
                
                stim_t, exc_input_list, successes_list, inh_input_list, output_list = (
                    process_patterned_input_simulation_input_output(source_file_path))
                seed_group.create_dataset('stim_t', data=stim_t, compression='gzip')
                seed_group.create_dataset('firing_rate', data=output_list[0], compression='gzip')
                
                _, filtered_i_syn_list_dict = get_patterned_input_filtered_synaptic_currents(source_file_path)
                for syn_type in filtered_i_syn_list_dict:
                    i_syn_key = 'i_%s' % syn_type
                    seed_group.create_dataset(i_syn_key, data=filtered_i_syn_list_dict[syn_type][0], compression='gzip')
                print('Processed and exported data from %s' % source_file_path)
                sys.stdout.flush()
    print('Finished exporting to %s' % target_file_path)
    sys.stdout.flush()


if __name__ == '__main__':
    main(standalone_mode=False)