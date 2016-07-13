__author__ = 'milsteina'
from plot_results import *
import sys
"""
The first 3 recs in the output file are from soma, proximal trunk, and distal trunk. The recordings then follow the
order [spine_vm, Pr, branch_vm] for 10 spines each in and out of field. Then i_AMPA and i_NMDA recordings have been
appended to the end of the list of the recording array. This script generates a grid of traces to compare in and out of
field activity.
"""

rec_filename = 'output051120161734-pid45981-seed3-e3200-i500-mod_inh0-track_spines_2.5_90'

spine_num = {}

if len(sys.argv) > 1:
    spine_num['out'] = int(sys.argv[1])
else:
    spine_num['out'] = 2
if len(sys.argv) > 2:
    spine_num['in'] = int(sys.argv[2])
else:
    spine_num['in'] = 4
if len(sys.argv) > 3:
    svg_title = str(sys.argv[3])
else:
    svg_title = None

f = h5py.File(data_dir+rec_filename+'.hdf5', 'r')
trial = f.itervalues().next()
equilibrate = trial.attrs['equilibrate']
track_equilibrate = trial.attrs['track_equilibrate']
duration = trial.attrs['duration']
track_duration = duration - equilibrate - track_equilibrate
dt = 0.02
t = np.arange(0., duration, 0.02)
start = int((equilibrate+track_equilibrate)/dt)

index_dict = {}
rec_key_int = 3

while rec_key_int < len(trial['rec']):
    rec_key = str(rec_key_int)
    description = trial['rec'][rec_key].attrs['description']
    index = trial['rec'][rec_key].attrs['index']
    if description == 'spine_vm':
        index_dict[index] = {}
        index_dict[index][description] = rec_key
        rec_key_int += 2
        rec_key = str(rec_key_int)
        description = trial['rec'][rec_key].attrs['description']
        index_dict[index][description] = rec_key  # branch_vm
        rec_key_int += 1
    else:
        index_dict[index][description] = rec_key  # i_AMPA or i_NMDA
        rec_key_int += 1
for train in trial['train'].itervalues():
    index = train.attrs['index']
    if index in index_dict:
        peak_loc = train.attrs['peak_loc']
        if 850. <= peak_loc <= 950.:
            index_dict[index]['peak_loc'] = 'out'
        elif 4450. <= peak_loc <= 4550.:
            index_dict[index]['peak_loc'] = 'in'
rec_key_dict = {}
rec_key_dict['in'] = [index for index in index_dict.itervalues() if index['peak_loc'] == 'in']
rec_key_dict['out'] = [index for index in index_dict.itervalues() if index['peak_loc'] == 'out']

xlim_dict = {}
xlim_dict['in'] = [4000., 5000.]
xlim_dict['out'] = [400., 1400.]

traces = {}

for loc in ['in', 'out']:
    traces[loc] = {}
    for trace in ['spine_vm', 'branch_vm', 'i_AMPA', 'i_NMDA']:
        rec_key = rec_key_dict[loc][spine_num[loc]][trace]
        this_trace = np.interp(t, trial['time'], trial['rec'][rec_key])
        this_trace = this_trace[start:]
        traces[loc][trace] = this_trace

soma_vm = np.interp(t, trial['time'], trial['rec']['0'])
soma_vm = soma_vm[start:]
t = t[start:]
t -= t[0]

if svg_title is not None:
    remember_font_size = mpl.rcParams['font.size']
    mpl.rcParams['font.size'] = 8

for loc in ['out', 'in']:
    fig, axes = plt.subplots(3, sharex=True)
    axes[0].plot(t, soma_vm, color='k', label='Soma Vm')
    clean_axes(axes[0])
    axes[0].tick_params(direction='out')
    axes[0].spines['bottom'].set_visible(False)
    axes[0].set_ylim(-70., -25.)
    axes[0].legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    axes[1].plot(t, traces[loc]['spine_vm'], color='c', label='Spine Vm')
    axes[1].plot(t, traces[loc]['branch_vm'], color='r', label='Branch Vm')
    clean_axes(axes[1])
    axes[1].tick_params(direction='out')
    axes[1].spines['bottom'].set_visible(False)
    axes[1].set_ylim(-70., -25.)
    axes[1].legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    axes[2].plot(t, traces[loc]['i_AMPA']*1000., color='g', label='i_AMPA', zorder=1)
    axes[2].plot(t, traces[loc]['i_NMDA']*1000., color='orange', label='i_NMDA', zorder=0)
    clean_axes(axes[2])
    axes[2].tick_params(direction='out')
    axes[2].set_ylim(-30., 5.)
    axes[2].legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    axes[2].set_xlim(xlim_dict[loc])
    if svg_title is not None:
        fig.set_size_inches(3.55, 3.2)
        fig.savefig(data_dir + svg_title + ' - ' + loc + ' field spine traces.svg', format='svg', transparent=True)
    plt.show()
    plt.close()

if svg_title is not None:
    mpl.rcParams['font.size'] = remember_font_size
f.close()

"""
rec_filename = 'output070416 - cell141 - e3200-i0-subtr0_inh0 - 10 trials'
f = h5py.File(data_dir+rec_filename+'.hdf5', 'r')
trial = f['3']
equilibrate = trial.attrs['equilibrate']
track_equilibrate = trial.attrs['track_equilibrate']
duration = trial.attrs['duration']
track_duration = duration - equilibrate - track_equilibrate
dt = 0.02
t = np.arange(0., duration, 0.02)
start = int((equilibrate+track_equilibrate)/dt)

vm = np.interp(t, trial['time'], trial['rec']['0'])
vm = vm[start:]
t = t[start:]
t -= t[0]

if svg_title is not None:
    remember_font_size = mpl.rcParams['font.size']
    mpl.rcParams['font.size'] = 8
fig, axes = plt.subplots(1)
axes.plot(t, vm, color='k')
clean_axes(axes)
axes.tick_params(direction='out')
#axes.spines['bottom'].set_visible(False)
axes.set_ylim(-70., 40.)
axes.set_ylabel('Vm (mV)')
axes.set_xlim(0., 7500.)
axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
axes.set_xlabel('Time (s)')
if svg_title is not None:
    fig.set_size_inches(3.55, 1.39)
    fig.savefig(data_dir + svg_title + ' - traces.svg', format='svg', transparent=True)
plt.show()
plt.close()
if svg_title is not None:
    mpl.rcParams['font.size'] = remember_font_size
f.close()
"""