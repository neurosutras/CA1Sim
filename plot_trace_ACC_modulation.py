from plot_results import *


rec_file_paths = {0: 'data/20200522_212514-pid55383-AAC_mod_inh1.00_10.hdf5',
                  1: 'data/20200522_212528-pid55404-AAC_mod_inh0.00_20.hdf5'}

t = None
vm_dict = dict()  # condition
labels = {0: 'Control', 1: 'No AAC Inhibition'}
colors = {0: 'k', 1: 'r'}

for condition, rec_file_path in rec_file_paths.items():
    with h5py.File(rec_file_path, 'r') as f:
        for sim in f.values():
            if t is None:
                equilibrate = sim.attrs['equilibrate']
                track_equilibrate = sim.attrs['track_equilibrate']
                dt = sim['time'].attrs['dt']
                t = sim['time'][:]
                start_index = np.where(t >= (equilibrate + track_equilibrate))[0][0]
                t = t[start_index:] - (equilibrate + track_equilibrate)
            for rec in sim['rec'].values():
                sec_type = rec.attrs['type']
                if sec_type == 'soma':
                    node_type = 'soma'
                    vm_dict[condition] = rec[:][start_index:]
                    break


fig, axes = plt.subplots()
for condition in vm_dict:
    axes.plot(t / 1000., vm_dict[condition], label=labels[condition], c=colors[condition], linewidth=0.75)
axes.legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
axes.set_ylim(-68., -45.)
axes.set_ylabel('Vm (mV)')
axes.set_xlabel('Time (s)')
clean_axes(axes)
fig.tight_layout()
fig.show()