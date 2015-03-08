__author__ = 'Aaron D. Milstein'
from IPython.parallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
import parallel_optimize_EPSP_amp_engine
"""
This simulation uses scipy.optimize to iterate through AMPA_KIN.gmax values to fit target EPSP amplitude at the soma.
Parallel version dynamically submits jobs to available cores.
"""
num_cores = parallel_optimize_EPSP_amp_engine.num_cores
c = Client()
v = c[:]
start_time = time.time()
v.execute('from parallel_optimize_EPSP_amp_engine import *')
v = c.load_balanced_view()
map_result = v.map_async(parallel_optimize_EPSP_amp_engine.optimize_single_synapse,
                         range(len(parallel_optimize_EPSP_amp_engine.syn_list)))
while not map_result.ready():
    clear_output()
    for stdout in [stdout for stdout in map_result.stdout if stdout][-num_cores:]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
    time.sleep(60)
results = map_result.get()
for stdout in [stdout for stdout in map_result.stdout if stdout][-num_cores:]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
print 'Parallel execution took:', time.time()-start_time, 's'
distances = {}
param_vals = {}
for result in results:
    if not result['sec_type'] in distances:
        distances[result['sec_type']] = []
        param_vals[result['sec_type']] = {}
for param_name in parallel_optimize_EPSP_amp_engine.param_names:
    for sec_type in param_vals:
        param_vals[sec_type][param_name] = []
for result in results:
    distances[result['sec_type']].append(result['distance'])
    for i, param_name in enumerate(parallel_optimize_EPSP_amp_engine.param_names):
        param_vals[result['sec_type']][param_name].append(result['result'][i])
new_rec_filename = '030815 kap_kad_ih_ampar_scale kd no_na optimize_EPSP_amp - EB1'
with h5py.File(data_dir+new_rec_filename+'.hdf5', 'w') as f:
    for sec_type in distances:
        f.create_group(sec_type)
        f.attrs['syn_type'] = parallel_optimize_EPSP_amp_engine.syn_type
        f[sec_type].create_dataset('distances', compression='gzip', compression_opts=9, data=distances[sec_type])
        for param_name in param_vals[sec_type]:
            f[sec_type].create_dataset(param_name, compression='gzip', compression_opts=9,
                                       data=param_vals[sec_type][param_name])
fig, axes = plt.subplots(1, max(2, len(distances)))
for i, sec_type in enumerate(distances):
    for param_name in param_vals[sec_type]:
        axes[i].scatter(distances[sec_type], param_vals[sec_type][param_name], label=param_name)
        axes[i].set_title(sec_type+' spines')
        axes[i].set_xlabel('Distance from Soma (um)')
        axes[i].set_ylabel('Peak Synaptic Conductance (uS)')
        axes[i].legend(loc='best')
plt.subplots_adjust(hspace=0.4, wspace=0.3, left=0.05, right=0.98, top=0.95, bottom=0.05)
plt.show()
plt.close()