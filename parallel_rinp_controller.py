__author__ = 'milsteina'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import parallel_rinp_engine
import sys
import os
"""
Parallel version: Iterates through every section, injecting hyperpolarizing current and measuring input resistance.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

new_rec_filename = parallel_rinp_engine.mech_filename+' - Rinp'

if len(sys.argv) > 1:
    cluster_id = sys.argv[1]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

dv = c[:]
dv.clear()
dv.block = True
start_time = time.time()
result = dv.execute('from parallel_rinp_engine import *')
while not result.ready():
    time.sleep(30)
v = c.load_balanced_view()

num_secs = len(parallel_rinp_engine.nodes)

result = v.map_async(parallel_rinp_engine.test_single_section, range(num_secs))
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
combine_output_files(rec_file_list, new_rec_filename)
for filename in rec_file_list:
    os.remove(data_dir+filename+'.hdf5')
#plot_Rinp(new_rec_filename)
