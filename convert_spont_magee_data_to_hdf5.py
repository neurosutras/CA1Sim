from plot_results import *
import pandas as pd
import sys

"""
Magee lab plateau-induced place field data is in a horrible ad-hoc series of text files. This organizes them into
a single .hdf5 file with a standard format:

cell_id:
    attrs:
        sampling_rate
        track_length
    induction_num:
        position:
            trial_num:
                data (raw)
                attrs:
                    induction_loc
                    induction_dur

        ramp:
            data (mean across laps)
"""


if len(sys.argv) > 1:
    cell_id = int(sys.argv[1])
else:
    raise Exception('No cell_id provided.')

output_file = '120216 magee lab spont'

cell_keys = {1: '130925s.txt',
             2: '131210s.txt',
             3: '140330s.txt',
             4: '140411s.txt',
             5: '140730s.txt',
             6: '141221s.txt',
             7: '160204s.txt'
             }

cell_legend = { # ramp had 97, not 100 spatial bins. Added 3 in the beginning
                1: [{'position': [0], 'ramp': 2, 'current': [1], 'induction_loc': [],
                    'induction_dur': []}],
                2: [{'position': [0], 'ramp': 2, 'current': [1], 'induction_loc': [],
                    'induction_dur': []}],
                3: [{'position': [1], 'ramp': 3, 'current': [2], 'induction_loc': [],
                    'induction_dur': []}],
                4: [{'position': [0], 'ramp': 2, 'current': [1], 'induction_loc': [],
                    'induction_dur': []}],
                5: [{'position': [1], 'ramp': 3, 'current': [2], 'induction_loc': [],
                    'induction_dur': []}],
                6: [{'position': [1], 'ramp': 3, 'current': [2], 'induction_loc': [],
                    'induction_dur': []}],
                7: [{'position': [0], 'ramp': 2, 'current': [1], 'induction_loc': [],
                    'induction_dur': []}]
               }


filename = cell_keys[cell_id]
df = pd.read_table(data_dir+'katie/'+filename, sep='\t', header=0)

all_data = {}
for c in range(len(df.values[0,:])):
    all_data[c] = df.values[:,c][~np.isnan(df.values[:,c])]

sampling_rate = 20000.  # Hz
# sampling_rate = 1000.  # Hz
track_length = 187.  # cm
V_scale = 1000.

for induction in cell_legend[cell_id]:
    for p in induction['position']:
        this_position = all_data[p]
        this_t = np.arange(0., len(this_position)*1000./sampling_rate, 1000./sampling_rate)
        plt.plot(this_t, this_position)
    plt.show()
    plt.close()

if 'current' in cell_legend[cell_id][0]:
    for induction in cell_legend[cell_id]:
        fig, axes = plt.subplots(2)
        for p, c in zip(induction['position'], induction['current']):
            this_position = all_data[p]
            all_data[c] = np.abs(all_data[c] - np.mean(all_data[c][:100]))
            this_current = all_data[c]
            this_t = np.arange(0., len(this_position) * 1000. / sampling_rate, 1000. / sampling_rate)
            axes[0].plot(this_position / np.max(this_position) * track_length, this_current)
            axes[1].plot(this_position / np.max(this_position) * track_length, this_t)
        plt.show()
        plt.close()
"""
if 'current' in cell_legend[cell_id][0]:
    for induction in cell_legend[cell_id]:
        for p, c in zip(induction['position'], induction['current']):
            this_position = all_data[p]
            this_current = all_data[c]
            start_index = np.where(this_current >= 0.5 * np.max(this_current))[0][0]
            this_induction_loc = this_position[start_index] / np.max(this_position) * track_length
            induction['induction_loc'].append(this_induction_loc)
            end_index = np.where(this_current[start_index:] <= 0.2 * np.max(this_current))[0][0]
            this_induction_dur = end_index * 1000./sampling_rate
            induction['induction_dur'].append(this_induction_dur)
            fig, axes = plt.subplots(2)
            axes[0].plot(this_position / np.max(this_position) * track_length, this_current)
            axes[1].plot(np.arange(0., track_length, track_length/100.),
                         all_data[cell_legend[cell_id][0]['ramp']][:100])
            plt.title('%.1f, %.1f' % (this_induction_loc, this_induction_dur))
            plt.show()
            plt.close()
"""
"""
with h5py.File(data_dir+output_file+'.hdf5', 'a') as f:
    f.create_group(str(cell_id))
    f[str(cell_id)].attrs['sampling_rate'] = sampling_rate
    f[str(cell_id)].attrs['track_length'] = track_length
    i = 1
    for induction in cell_legend[cell_id]:
        f[str(cell_id)].create_group(str(i))
        f[str(cell_id)][str(i)].create_group('position')
        for j, column in enumerate(induction['position']):
            f[str(cell_id)][str(i)]['position'].create_dataset(str(j), compression='gzip', compression_opts=9,
                                                               data=all_data[column]*10.)
            if type(induction['induction_loc']) != list:
                f[str(cell_id)][str(i)]['position'][str(j)].attrs['induction_loc'] = induction['induction_loc']
                f[str(cell_id)][str(i)]['position'][str(j)].attrs['induction_dur'] = induction['induction_dur']
            else:
                f[str(cell_id)][str(i)]['position'][str(j)].attrs['induction_loc'] = induction['induction_loc'][j]
                f[str(cell_id)][str(i)]['position'][str(j)].attrs['induction_dur'] = induction['induction_dur'][j]
        column = induction['ramp']
        f[str(cell_id)][str(i)].create_dataset('ramp', compression='gzip', compression_opts=9,
                                               data=all_data[column][:100] * V_scale)
        i += 1
"""

"""
with h5py.File(data_dir+output_file+'.hdf5', 'r') as f:
    for cell in f.itervalues():
        sampling_rate = cell.attrs['sampling_rate']
        track_length = cell.attrs['track_length']
        for induction in cell.itervalues():
            fig = plt.figure(1)
            for position in induction['position'].itervalues():
                plt.plot(np.arange(0., len(position)/sampling_rate, 1./sampling_rate), position)
            fig = plt.figure(2)
            plt.plot(np.arange(0., track_length, track_length/100.), induction['ramp'])
        plt.show()
        plt.close()
"""