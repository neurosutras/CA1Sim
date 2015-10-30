__author__ = 'milsteina'
from plot_results import *


rec_filenames = ['102715 test simple_axon_model_no_na',
                 '102715 test simple_axon_model_no_na_reduced_k',
                 '102715 test simple_axon_model']

figure_title = '102815 simple axon model'

titles = ['0% gNa, 100% gK', '0% gNa, 30% gK', '100% gNa, 100% gK']
propagation = {}
dt = 0.01

for i, title in enumerate(titles):
    distances, propagation[title] = process_simple_axon_model_output(rec_filenames[i])

colors = ['k', 'r', 'b']

target_voltages = [-5.0, 15.0, 40.0]

for target_voltage in target_voltages:
    for i, title in enumerate(titles):
        plt.scatter(distances, propagation[title][target_voltage], color=colors[i], label=title)
    plt.legend(loc='best', frameon=False, framealpha=0.5)
    plt.xlabel('Distance From Soma (um)')
    plt.ylabel('Voltage Propagation')
    plt.title('Somatic Voltage Amplitude ~ %i mV' % target_voltage)
    #plt.savefig(data_dir+figure_title+' - voltage prop amp %i.svg' % target_voltage, format='svg')
    plt.show()
    plt.close()

for i, title in enumerate(titles[:-1]):
    with h5py.File(data_dir+rec_filenames[i]+'.hdf5', 'r') as f:
        sim = f['2']
        equilibrate = sim.attrs['equilibrate']
        duration = sim.attrs['duration']
        stim_dur = sim.attrs['stim_dur']
        t = np.arange(0., duration, dt)
        soma_vm = np.interp(t, sim['time'], sim['rec']['0'])
        soma_rest = np.mean(soma_vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
        soma_vm -= soma_rest
        axon_vm = np.interp(t, sim['time'], sim['rec']['6'])  # 175 um from soma
        axon_rest = np.mean(axon_vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
        axon_vm -= axon_rest
        left = int((equilibrate - 20.)/dt)
        right = int((equilibrate + stim_dur + 40.)/dt)
        plt.plot(t[left:right], soma_vm[left:right], color='k', label='Soma')
        plt.plot(t[left:right], axon_vm[left:right], color='silver', label='Axon (175 um)')
        plt.ylim(-10., 120.)
        plt.title(title)
        plt.legend(loc='best', frameon=False, framealpha=0.5)
        plt.savefig(data_dir+figure_title+' - traces - %s.svg' % title, format='svg')
        plt.show()
        plt.close()

"""
rec_filename = '102815 test simple_axon_model_spike_traces'
with '100% gNa, 100% gK' as title:
    with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
        sim = f['0']
        equilibrate = sim.attrs['equilibrate']
        duration = sim.attrs['duration']
        stim_dur = sim.attrs['stim_dur']
        t = np.arange(0., duration, dt)
        soma_vm = np.interp(t, sim['time'], sim['rec']['0'])
        soma_rest = np.mean(soma_vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
        soma_vm -= soma_rest
        for rec, label in zip(['1', '3', '6'], ['AIS (25 um)', 'Axon (85 um)', 'Axon (175 um)']):
            axon_vm = np.interp(t, sim['time'], sim['rec'][rec])
            axon_rest = np.mean(axon_vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            axon_vm -= axon_rest
            left = int((equilibrate - 20.)/dt)
            right = int((equilibrate + stim_dur + 40.)/dt)
            plt.plot(t[left:right], soma_vm[left:right], color='k', label='Soma')
            plt.plot(t[left:right], axon_vm[left:right], color='silver', label=label)
            plt.ylim(-10., 120.)
            plt.title(title)
            plt.legend(loc='best', frameon=False, framealpha=0.5)
            plt.savefig(data_dir+figure_title+' - traces - %s.svg' % label, format='svg')
            plt.show()
            plt.close()
"""