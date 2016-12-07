__author__ = 'milsteina'
import parallel_rinp_engine
import sys
import os
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *

"""
Parallel version: Iterates through every section, injecting hyperpolarizing current and measuring input resistance.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    mech_filename = '110316 DG_GC pas spines'
if len(sys.argv) > 3:
    description = str(sys.argv[3])
    new_rec_filename = datetime.datetime.today().strftime('%m%d%Y%H%M') + '_Rinp_' + description
else:
    new_rec_filename = datetime.datetime.today().strftime('%m%d%Y%H%M') + '_Rinp'
if len(sys.argv) > 4:
    cluster_id = sys.argv[4]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()


dv = c[:]
dv.clear()
dv.block = True
start_time = time.time()
dv.execute('run parallel_rinp_engine %i \"%s\"' % (int(spines), mech_filename))
# time.sleep(120)
v = c.load_balanced_view()

num_secs = len(parallel_rinp_engine.nodes)

result = v.map_async(parallel_rinp_engine.test_single_section, range(num_secs))
#result = v.map_async(parallel_rinp_engine.test_single_section, range(40))
while not result.ready():
    time.sleep(30)
    clear_output()
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
print 'Parallel execution took: ', time.time()-start_time, ' s'
rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir+filename+'.hdf5')]

with h5py.File(data_dir+new_rec_filename+'.hdf5', 'w') as f:
    f.create_group("Rinp_data")
    f.create_group("avg_waves")
    rec_dict = {sec_type:[] for sec_type in ['axon', 'soma']}
    #create group for avg waves
    trial_ind = 0
    for filename in rec_file_list:
        with h5py.File(data_dir + filename + '.hdf5', 'r') as r:
            for trial in r.itervalues():
                f["Rinp_data"].create_group(str(trial_ind))
                for key, value in trial.attrs.iteritems():
                    f["Rinp_data"][str(trial_ind)].attrs[key] = value
                    #'Rinp_peak', 'decay_50', 'Rinp_baseline', 'Rinp_steady'
                rec = trial['rec'].itervalues().next()
                for key, value in rec.attrs.iteritems():
                    f["Rinp_data"][str(trial_ind)].attrs[key] = value
                    # adds 'cell', 'index', 'type', 'loc', 'soma_distance', 'branch_distance'
                sec_type = rec.attrs['type']
                if sec_type in ['ais', 'axon_hill']:
                    key = 'axon'
                elif sec_type in ['apical', 'trunk', 'tuft', 'basal']:
                    distance = rec.attrs['soma_distance']
                    if distance < 200.:
                        key = 'prox_'+sec_type
                    else:
                        key = 'dist_'+sec_type
                else:
                    key = sec_type
                if key not in rec_dict:
                    rec_dict[key] = []
                dt = 0.02
                interp_t, interp_vm = interpolate_tvec_vec(trial['time'], rec, parallel_rinp_engine.duration, dt)
                if trial_ind == 0:
                    f["avg_waves"].create_dataset('time', compression='gzip', compression_opts=9,
                                                  data=interp_t[int(parallel_rinp_engine.equilibrate/dt):] -
                                                  parallel_rinp_engine.equilibrate - parallel_rinp_engine.delay)
                #need to fix interpolation here
                rec_dict[key].append(interp_vm[int(parallel_rinp_engine.equilibrate/dt):])
                trial_ind += 1

    #average waveforms for soma, prox_dend, axon, dist_dend
    for key in rec_dict:
        if np.any(rec_dict[key]):
            f["avg_waves"].create_dataset(key, compression='gzip', compression_opts=9,
                                                            data=(np.mean(rec_dict[key], axis=0)))


for filename in rec_file_list:
    os.remove(data_dir+filename+'.hdf5')

plot_Rinp_general(new_rec_filename)
