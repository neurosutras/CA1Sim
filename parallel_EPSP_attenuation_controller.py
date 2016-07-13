__author__ = 'milsteina'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import parallel_EPSP_attenuation_engine
import sys
"""
Parallel version: Iterates through every spine and activates AMPA_KIN synapses. Allows measurement of EPSP attenuation
and kinetics. Mechanism dictionary specifying gmax gradients is specified in parallel_EPSP_attenuation_engine.py

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
filename_suffix = ' - epsp attenuation'
#condition = ' nmda no_ih'
#new_rec_filename = parallel_EPSP_attenuation_engine.mech_filename+condition+filename_suffix
new_rec_filename = parallel_EPSP_attenuation_engine.mech_filename+filename_suffix

num_syns = len(parallel_EPSP_attenuation_engine.syn_list)

if len(sys.argv) > 1:
    cluster_id = sys.argv[1]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

dv = c[:]
dv.clear()
dv.block = True
start_time = time.time()
dv.execute('from parallel_EPSP_attenuation_engine import *')
v = c.load_balanced_view()
result = v.map_async(parallel_EPSP_attenuation_engine.stimulate_single_synapse, range(num_syns))
while not result.ready():
    time.sleep(30)
    clear_output()
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
print 'Parallel execution took: %.3f s' % (time.time()-start_time)
rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir+filename+'.hdf5')]
combine_output_files(rec_file_list, new_rec_filename)
for filename in rec_file_list:
    os.remove(data_dir+filename+'.hdf5')
#plot_EPSP_attenuation(new_rec_filename)
#plot_EPSP_kinetics(new_rec_filename)