__author__ = 'milsteina'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import parallel_spine_attenuation_ratio_engine
import sys
"""
Parallel version: Iterates through spines injecting an EPSC-shaped current and recording from the spine and branch to
calculate the amplitude attenuation ratio, spine neck, and branch impedance. Mechanism dictionary is specified in
parallel_spine_attenuation_ratio_engine.py

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
filename_suffix = ' - spine AR'
new_rec_filename = parallel_spine_attenuation_ratio_engine.mech_filename+filename_suffix

num_syns = len(parallel_spine_attenuation_ratio_engine.spine_syn_list)

if len(sys.argv) > 1:
    cluster_id = sys.argv[1]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

dv = c[:]
dv.clear()
dv.block = True
start_time = time.time()
dv.execute('from parallel_spine_attenuation_ratio_engine import *')
v = c.load_balanced_view()
result = v.map_async(parallel_spine_attenuation_ratio_engine.calculate_single_attenuation_ratio, range(num_syns))
while not result.ready():
    time.sleep(30)
    clear_output()
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
print 'Parallel execution took:', time.time()-start_time, 's'
rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir+filename+'.hdf5')]
combine_output_files(rec_file_list, new_rec_filename)
for filename in rec_file_list:
    os.remove(data_dir+filename+'.hdf5')
#plot_AR(new_rec_filename)