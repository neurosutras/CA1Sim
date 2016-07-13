__author__ = 'Aaron D. Milstein'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
import parallel_clustered_branch_cooperativity_nmda_engine_110315
import os
"""
This simulation steps through a list of grouped_spines, and saves output from stimulating each spine (expected) and
group of spines (actual) to generate input-output plots. Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
new_rec_filename = '020616 clustered nmda cooperativity'

if len(sys.argv) > 1:
    cluster_id = sys.argv[1]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()
dv = c[:]
dv.clear()
dv.block = True
global_start_time = time.time()
dv.execute('from parallel_clustered_branch_cooperativity_nmda_engine_110315 import *')
time.sleep(300)
v = c.load_balanced_view()

start_time = time.time()

instructions = []
for i in range(len(parallel_clustered_branch_cooperativity_nmda_engine_110315.groups_to_stim)):
    for j in range(len(parallel_clustered_branch_cooperativity_nmda_engine_110315.groups_to_stim[i]['spines'])):
        instructions.append((i, j))
result = v.map_async(parallel_clustered_branch_cooperativity_nmda_engine_110315.stim_single_expected, instructions)
while not result.ready():
    time.sleep(30)
    clear_output()
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir+filename+'.hdf5')]
combine_output_files(rec_file_list, new_rec_filename+'_expected')
for filename in rec_file_list:
    os.remove(data_dir+filename+'.hdf5')

instructions = []
for i in range(len(parallel_clustered_branch_cooperativity_nmda_engine_110315.groups_to_stim)):
    for j in range(1, len(parallel_clustered_branch_cooperativity_nmda_engine_110315.groups_to_stim[i]['spines'])+1):
        instructions.append((i, j))
result = v.map_async(parallel_clustered_branch_cooperativity_nmda_engine_110315.stim_actual_group, instructions)
while not result.ready():
    time.sleep(30)
    clear_output()
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir+filename+'.hdf5')]
combine_output_files(rec_file_list, new_rec_filename+'_actual')
for filename in rec_file_list:
    os.remove(data_dir+filename+'.hdf5')
print 'Parallel simulation took %i s to stimulate %i groups of spines' % (time.time() - start_time,
                                        len(parallel_clustered_branch_cooperativity_nmda_engine_110315.groups_to_stim))