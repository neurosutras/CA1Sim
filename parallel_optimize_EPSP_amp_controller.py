__author__ = 'Aaron D. Milstein'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
import parallel_optimize_EPSP_amp_engine
"""
This simulation uses scipy.optimize to iterate through AMPA_KIN.gmax values to fit target EPSP amplitude at the soma.
Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
filename_suffix = ' - AMPAR_scaling'
new_rec_filename = parallel_optimize_EPSP_amp_engine.mech_filename+filename_suffix

c = Client()
dv = c[:]
dv.clear()
dv.block = True
start_time = time.time()
dv.execute('from parallel_optimize_EPSP_amp_engine import *')
v = c.load_balanced_view()
map_result = v.map_async(parallel_optimize_EPSP_amp_engine.optimize_single_synapse,
                         range(len(parallel_optimize_EPSP_amp_engine.syn_list)))
while not map_result.ready():
    time.sleep(30)
    for stdout in [stdout for stdout in map_result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
    clear_output()
results = map_result.get()
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
with h5py.File(data_dir+new_rec_filename+'.hdf5', 'w') as f:
    f.attrs['syn_type'] = parallel_optimize_EPSP_amp_engine.syn_type
    for i, param_name in enumerate(parallel_optimize_EPSP_amp_engine.param_names):
        f.attrs[param_name] = parallel_optimize_EPSP_amp_engine.param_ylabels[i]
    for sec_type in distances:
        f.create_group(sec_type)
        f[sec_type].create_dataset('distances', compression='gzip', compression_opts=9, data=distances[sec_type])
        for param_name in param_vals[sec_type]:
            f[sec_type].create_dataset(param_name, compression='gzip', compression_opts=9,
                                       data=param_vals[sec_type][param_name])
plot_synaptic_parameter(new_rec_filename)