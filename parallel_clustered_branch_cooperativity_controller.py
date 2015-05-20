__author__ = 'milsteina'
from IPython.parallel import Client
from IPython.display import clear_output
from plot_results import *
import parallel_branch_cooperativity_engine
import sys
"""
Parallel version: Branches are organized into groups based on paths to the soma or trunk. Iterates through every
branch path, activating increasing numbers of spines containing AMPA_KIN and NMDA_KIN mechanisms. Allows comparison of
Expected and Actual EPSPs.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
filename_suffix = ' - actual EPSP - clustered basal sample'
new_rec_filename = parallel_branch_cooperativity_engine.mech_filename+filename_suffix
#new_rec_filename = '051915 test bash script - branch cooperativity'

num_paths = len(parallel_branch_cooperativity_engine.path_list)
#c = Client(profile='mpi')
c = Client()
dv = c[:]
dv.clear()
dv.block = True
start_time = time.time()
dv.execute('from parallel_branch_cooperativity_engine import *')
v = c.load_balanced_view()
instructions = []
for path_index in range(1):  # num_paths):
    #for num_syns in range(1, len(parallel_branch_cooperativity_engine.path_list[path_index]['spines'])):
    for num_syns in range(1, min(len(parallel_branch_cooperativity_engine.path_list[path_index]['spines'])+1,
                                 len(c)+1)):
        instructions.append((path_index, num_syns))
result = v.map_async(parallel_branch_cooperativity_engine.stimulate_synapse_group, instructions)
while not result.ready():
    clear_output()
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
    time.sleep(60)
print 'Parallel execution took: %.3f s' % (time.time()-start_time)
#print result.get()
rec_file_list = dv['rec_filename']
combine_output_files(rec_file_list, new_rec_filename)
#plot_superimpose_conditions(new_rec_filename)