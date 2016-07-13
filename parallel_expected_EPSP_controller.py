__author__ = 'milsteina'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import parallel_expected_EPSP_engine
import sys
"""
Parallel version: Iterates through every spine and activates synapses containing AMPA_KIN and NMDA_KIN mechanisms, for
comparing Expected to Actual depolarization.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
filename_suffix = ' - expected EPSP'
new_rec_filename = parallel_expected_EPSP_engine.mech_filename+filename_suffix

num_syns = len(parallel_expected_EPSP_engine.syn_list)
c = Client()
dv = c[:]
dv.clear()
dv.block = True
start_time = time.time()
dv.execute('from parallel_expected_EPSP_engine import *')
v = c.load_balanced_view()
result = v.map_async(parallel_expected_EPSP_engine.stimulate_single_synapse, range(num_syns))
while not result.ready():
    clear_output()
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
    time.sleep(60)
print 'Parallel execution took: %.3f s' % (time.time()-start_time)
rec_file_list = dv['rec_filename']
combine_output_files(rec_file_list, new_rec_filename)