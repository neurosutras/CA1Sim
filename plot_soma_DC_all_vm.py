from plot_results import *
from collections import defaultdict

rec_file_paths = {0: 'data/output012920202140-pid41451-seed0-e3200-i600-condition0-seed_0.hdf5',
                  1: 'data/output012920202140-pid41470-seed0-e3200-i600-condition1-seed_1.hdf5',
                  2: 'data/output012920202140-pid41508-seed0-e3200-i600-condition2-seed_2.hdf5',
                  3: 'data/output030920201109-pid71331-seed0-e3200-i600-condition3-seed_3.hdf5' }

t = None
vm_dict = defaultdict(lambda: defaultdict(dict))  # node_type: distance: condition
labels = {0: 'Silent', 1: 'Silent + depo', 2: 'Place', 3: 'Place + hyper'}
colors = {0: 'k', 1: 'purple', 2: 'c', 3: 'r'}

for condition, rec_file_path in rec_file_paths.items():
    with h5py.File(rec_file_path, 'r') as f:
        for sim in f.values():
            if t is None:
                equilibrate = sim.attrs['equilibrate']
                track_equilibrate = sim.attrs['track_equilibrate']
                duration = sim.attrs['duration'] - (equilibrate + track_equilibrate)
                dt = sim['time'].attrs['dt']
                t = sim['time'][:]
                start_index = np.where(t >= (equilibrate + track_equilibrate))[0][0]
                t = t[start_index:] - (equilibrate + track_equilibrate)
            for rec in sim['rec'].values():
                description = rec.attrs['description']
                sec_type = rec.attrs['type']
                distance = rec.attrs['soma_distance']
                if sec_type == 'soma':
                    node_type = 'soma'
                elif sec_type in ['apical', 'basal']:
                    node_type = '%s_dend' % sec_type
                    if sec_type == 'basal':
                        distance *= -1.
                elif sec_type == 'spine_head':
                    if 'apical' in description:
                        node_type = 'apical_spine'
                    elif 'basal' in description:
                        node_type = 'basal_spine'
                        distance *= -1.
                else:
                    continue
                vm_dict[node_type][distance][condition] = rec[:][start_index:]


down_dt = 0.5
down_t = np.arange(0., duration, down_dt)
# 2000 ms Hamming window, ~3 Hz low-pass filter
window_len = int(2000./down_dt)
pad_len = int(window_len / 2.)
ramp_filter = signal.firwin(window_len, 2., nyq=1000. / 2. / down_dt)
filtered_vm_dict = defaultdict(lambda: defaultdict(dict))  # node_type: distance: condition
mean_filtered_vm_dict = defaultdict(lambda: defaultdict(dict))  # node_type: distance: condition

for node_type in vm_dict:
    for distance in vm_dict[node_type]:
        for condition in vm_dict[node_type][distance]:
            vm = vm_dict[node_type][distance][condition]
            down_sampled = np.interp(down_t, t, vm)
            padded_trace = np.zeros(len(down_sampled) + window_len)
            padded_trace[pad_len:-pad_len] = down_sampled
            padded_trace[:pad_len] = down_sampled[::-1][-pad_len:]
            padded_trace[-pad_len:] = down_sampled[::-1][:pad_len]
            filtered = signal.filtfilt(ramp_filter, [1.], padded_trace, padlen=pad_len)
            filtered = filtered[pad_len:-pad_len]
            filtered_vm_dict[node_type][distance][condition] = filtered
            mean_filtered_vm_dict[node_type][distance][condition] = np.mean(filtered)

mean_filtered_vm_array_dict = defaultdict(lambda: defaultdict(lambda: {'distance': [], 'value': []}))
for node_type in mean_filtered_vm_dict:
    for distance in mean_filtered_vm_dict[node_type]:
        for condition in mean_filtered_vm_dict[node_type][distance]:
            val = mean_filtered_vm_dict[node_type][distance][condition]
            if node_type == 'soma':
                mean_filtered_vm_array_dict['Dendrites'][condition]['distance'].append(distance)
                mean_filtered_vm_array_dict['Dendrites'][condition]['value'].append(val)
                mean_filtered_vm_array_dict['Spines'][condition]['distance'].append(distance)
                mean_filtered_vm_array_dict['Spines'][condition]['value'].append(val)
            elif 'dend' in node_type:
                mean_filtered_vm_array_dict['Dendrites'][condition]['distance'].append(distance)
                mean_filtered_vm_array_dict['Dendrites'][condition]['value'].append(val)
            elif 'spine' in node_type:
                mean_filtered_vm_array_dict['Spines'][condition]['distance'].append(distance)
                mean_filtered_vm_array_dict['Spines'][condition]['value'].append(val)

for node_type in mean_filtered_vm_array_dict:
    for condition in mean_filtered_vm_array_dict[node_type]:
        distances = np.array(mean_filtered_vm_array_dict[node_type][condition]['distance'])
        values = np.array(mean_filtered_vm_array_dict[node_type][condition]['value'])
        indexes = np.argsort(distances)
        mean_filtered_vm_array_dict[node_type][condition]['distance'] = distances[indexes]
        mean_filtered_vm_array_dict[node_type][condition]['value'] = values[indexes]

fig, axes = plt.subplots(2, 2, figsize=(12., 7.))
for i, node_type in enumerate(mean_filtered_vm_array_dict):
    for condition in mean_filtered_vm_array_dict[node_type]:
        axes[0][i].plot(mean_filtered_vm_array_dict[node_type][condition]['distance'],
                        mean_filtered_vm_array_dict[node_type][condition]['value'], c=colors[condition], marker='o')
    axes[0][i].set_ylabel('Mean Filtered Vm (mV)')
    axes[0][i].set_ylim(-90., -45.)
    axes[0][i].set_yticks(np.arange(-85., -40., 10.))
    axes[0][i].set_xticks(np.arange(-150., 300., 50.))
    axes[0][i].set_xlabel('Distance from soma (um)')
    axes[0][i].set_title(node_type)
clean_axes(axes)
fig.tight_layout()
fig.show()


for plot_filtered in [False, True]:
    cols = 3
    rows = 2
    fig, axes = plt.subplots(rows, cols, figsize=(12., 7.))
    count = 0
    for node_type in vm_dict:
        for distance in vm_dict[node_type]:
            this_axis = axes[count // cols][count % cols]
            this_axis.set_ylim(-90., -20.)
            this_axis.set_ylabel('Vm (mV)')
            this_axis.set_xlabel('Time (ms)')
            this_axis.set_title('%s (%i um)' % (node_type, int(distance)))
            for condition in vm_dict[node_type][distance]:
                if plot_filtered:
                    this_axis.plot(down_t, filtered_vm_dict[node_type][distance][condition], label=labels[condition],
                                   c=colors[condition], linewidth=0.75)
                else:
                    this_axis.plot(t, vm_dict[node_type][distance][condition], label=labels[condition],
                                   c=colors[condition], linewidth=0.75)
            if count == 0:
                this_axis.legend(loc='best', frameon=False)
            if count == (cols * rows - 1):
                clean_axes(axes)
                fig.tight_layout()
                fig.show()
                fig, axes = plt.subplots(rows, cols, figsize=(12., 7.))
                count = 0
            else:
                count += 1
    axes[count // cols][count % cols].legend(loc='best', frameon=False)
    clean_axes(axes)
    fig.tight_layout()
    fig.show()