__author__ = 'Aaron D. Milstein' and 'Grace Ng'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
import parallel_optimize_EPSP_amp_engine_GC
"""
This simulation uses scipy.optimize to iterate through AMPA_KIN.gmax values to fit target EPSP amplitude at the soma.
Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
filename_suffix = ' - AMPAR_scaling'
new_rec_filename = parallel_optimize_EPSP_amp_engine_GC.mech_filename+filename_suffix

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    # mech_filename = '042117 GC optimizing spike stability'
    mech_filename = '042617 GC optimizing spike stability'
if len(sys.argv) > 3:
    cluster_id = sys.argv[3]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

dv = c[:]
dv.block = True
start_time = time.time()
dv.execute('run parallel_optimize_EPSP_amp_engine_GC %i \"%s\"' % (int(spines), mech_filename))
v = c.load_balanced_view()
result = v.map_async(parallel_optimize_EPSP_amp_engine_GC.optimize_single_synapse,
                         range(len(parallel_optimize_EPSP_amp_engine_GC.syn_list)))
while not result.ready():
    time.sleep(30)
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
    clear_output()
results = result.get()
print 'Parallel execution took:', time.time()-start_time, 's'
distances = {}
param_vals = {}
for result in results:
    if not result['sec_type'] in distances:
        distances[result['sec_type']] = []
        param_vals[result['sec_type']] = {}
for param_name in parallel_optimize_EPSP_amp_engine_GC.param_names:
    for sec_type in param_vals:
        param_vals[sec_type][param_name] = []
for result in results:
    distances[result['sec_type']].append(result['distance'])
    for i, param_name in enumerate(parallel_optimize_EPSP_amp_engine_GC.param_names):
        param_vals[result['sec_type']][param_name].append(result['result'][i])
with h5py.File(data_dir+new_rec_filename+'.hdf5', 'w') as f:
    f.attrs['syn_type'] = parallel_optimize_EPSP_amp_engine_GC.syn_type
    for i, param_name in enumerate(parallel_optimize_EPSP_amp_engine_GC.param_names):
        f.attrs[param_name] = parallel_optimize_EPSP_amp_engine_GC.param_ylabels[i]
    for sec_type in distances:
        f.create_group(sec_type)
        f[sec_type].create_dataset('distances', compression='gzip', compression_opts=9, data=distances[sec_type])
        for param_name in param_vals[sec_type]:
            f[sec_type].create_dataset(param_name, compression='gzip', compression_opts=9,
                                       data=param_vals[sec_type][param_name])
plot_synaptic_parameter(new_rec_filename)