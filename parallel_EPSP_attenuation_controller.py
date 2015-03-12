__author__ = 'milsteina'
from IPython.parallel import Client
from IPython.display import clear_output
from plot_results import *
import parallel_EPSP_attenuation_engine
import sys
"""
Parallel version: Iterates through every spine, activating AMPA_KIN synapses and measuring EPSP attenuation.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
new_rec_filename = '031215 kap_kad_ih_ampar_scale kd no_na - EB2 - epsp_attenuation sample'

num_syns = len(parallel_EPSP_attenuation_engine.syn_list)
c = Client()
v = c[:]
start_time = time.time()
v.execute('from parallel_EPSP_attenuation_engine import *')
result = v.map_async(parallel_EPSP_attenuation_engine.stimulate_single_synapse, range(num_syns))
while not result.ready():
    clear_output()
    for stdout in result.stdout:
        if stdout:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
    sys.stdout.flush()
    time.sleep(60)
rec_file_list = result.get()
for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
    lines = stdout.split('\n')
    if lines[-2]:
        print lines[-2]
print 'Parallel execution took:', time.time()-start_time, 's'
combine_output_files(rec_file_list, new_rec_filename)
plot_EPSP_attenuation(new_rec_filename)
plot_EPSP_kinetics(new_rec_filename)