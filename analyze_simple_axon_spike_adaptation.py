__author__ = 'milsteina'
from function_lib import *

rec_filename = '102815 test simple_axon_model_spike_height'

dt = 0.01
th_dvdt = 10.

with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
    equilibrate = f.values()[0].attrs['equilibrate']
    duration = f.values()[0].attrs['duration']
    offset = int((equilibrate+0.4)/dt)
    t = np.arange(0., duration, dt)
    soma_vm = np.interp(t, f['0']['time'], f['0']['rec']['0'])
    soma_vm = soma_vm[offset:]
    dvdt = np.gradient(soma_vm, [dt])
    th_x = np.where(dvdt > th_dvdt)[0]
    spike_peak_locs = []
    distances = []
    spike_amp_array = []
    if not th_x.any():
        print 'Does not cross threshold.'
    else:
        i = 0
        while i < len(th_x):
            start = th_x[i]
            stop = th_x[i] + int(4./dt)
            peak = np.max(soma_vm[start:stop])
            peak_loc = np.where(soma_vm[start:stop] == peak)[0][0] + start
            spike_peak_locs.append(peak_loc + offset)
            i = np.where(th_x > peak_loc)[0]
            if i.any():
                i = i[0]
            else:
                break
        for i, peak_loc in enumerate(spike_peak_locs):
            amp_list = []
            for rec in f['0']['rec'].values():
                if i == 0:
                    distances.append(rec.attrs['soma_distance'])
                vm = np.interp(t, f['0']['time'], rec)
                vm_rest = np.mean(vm[int((equilibrate-3.)/dt):int((equilibrate-1.)/dt)])
                start = peak_loc - int(1.5/dt)
                stop = peak_loc + int(2.5/dt)
                peak = np.max(vm[start:stop])
                amp = peak - vm_rest
                amp_list.append(amp)
            spike_amp_array.append(amp_list)

plt.scatter(distances, np.divide(spike_amp_array[3], spike_amp_array[0]), c='k')
plt.ylabel('Spike Height Attenuation (4th Spike / 1st Spike)')
plt.xlabel('Distance From Soma (um)')
plt.show()