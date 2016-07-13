__author__ = 'Aaron D. Milstein'
from ipyparallel import Client
from IPython.display import clear_output
from function_lib import *
import sys
import build_expected_EPSP_reference_engine
import os
import time
"""
This simulation steps through a list of spines, and saves output from stimulating each spine in isolation, including
location-specific changes in the weights of inputs during patterened input simulation. Can be used to compare expected
and actual somatic depolarization. Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

if len(sys.argv) > 1:
    cluster_id = sys.argv[1]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()
if len(sys.argv) > 2:
    synapses_seed = int(sys.argv[2])
else:
    synapses_seed = 1

num_exc_syns = 2900
num_inh_syns = 500

new_rec_filename = '021016 expected reference'+'-seed'+str(synapses_seed)+'-e'+str(num_exc_syns)+'-i'+str(num_inh_syns)

dv = c[:]
dv.clear()
dv.block = True
start_time = time.time()
result = dv.execute('from build_expected_EPSP_reference_engine import *')
while not result.ready():
    time.sleep(30)
result = dv.execute('local_container.distribute_synapses('+str(synapses_seed)+', '+str(num_exc_syns)+', '+
                    str(num_inh_syns)+')')
while not result.ready():
    time.sleep(30)

v = c.load_balanced_view()
num_exc_syns = dv['len(local_container.stim_exc_syn_list)'][0]

result = v.map_async(build_expected_EPSP_reference_engine.stim_single_exc_syn, range(num_exc_syns))
while not result.ready():
    time.sleep(30)
    clear_output()
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir+filename+'.hdf5')]
combine_output_files(rec_file_list, new_rec_filename)
for filename in rec_file_list:
    os.remove(data_dir+filename+'.hdf5')
print 'Parallel simulation took %i s to stimulate %i synapses' % (time.time() - start_time, num_exc_syns)