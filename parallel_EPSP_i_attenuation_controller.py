__author__ = 'milsteina'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import parallel_EPSP_i_attenuation_engine
import sys
"""
Parallel version: Iterates through every node and injects EPSC-shaped currents. Allows measurement of EPSP_i attenuation
and kinetics. Mechanism dictionary specifying gmax gradients is specified in parallel_EPSP_i_attenuation_engine.py

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

filename_suffix = ' - epsp_i attenuation'
#new_rec_filename = '042015 soma_pas spine_adjusted - EB2 - epsp_i attenuation'
#new_rec_filename = '042015 soma_pas kdr ka_scale - adjusted - EB2 - epsp_i attenuation'
#new_rec_filename = parallel_EPSP_i_attenuation_engine.mech_filename + filename_suffix
new_rec_filename = '101415 test changes to epsp_i'

num_syns = len(parallel_EPSP_i_attenuation_engine.nodes)
c = Client()
dv = c[:]
dv.clear()
dv.block = True
start_time = time.time()
dv.execute('from parallel_EPSP_i_attenuation_engine import *')
v = c.load_balanced_view()
result = v.map_async(parallel_EPSP_i_attenuation_engine.stimulate_single_synapse, range(num_syns))
while not result.ready():
    time.sleep(30)
    clear_output()
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
print 'Parallel execution took:', time.time()-start_time, 's'
rec_file_list = dv['rec_filename']
combine_output_files(rec_file_list, new_rec_filename)
for filename in rec_file_list:
    os.remove(data_dir+filename+'.hdf5')
#plot_EPSP_i_attenuation(new_rec_filename)
plot_EPSP_i_kinetics(new_rec_filename)
#plot_EPSP_i_vm(new_rec_filename)