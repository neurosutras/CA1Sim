__author__ = 'Aaron D. Milstein and Grace Ng'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
import parallel_clustered_branch_cooperativity_nmda_engine_GC
import os
"""
This simulation steps through a list of grouped_synapses, and saves output from stimulating each synapse (expected) and
group of synapses (actual) to generate input-output plots. Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
filename_suffix = ' - clustered nmda cooperativity'
new_rec_filename = parallel_clustered_branch_cooperativity_nmda_engine_GC.mech_filename+filename_suffix

if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    # mech_filename = '042617 GC optimizing spike stability'
    mech_filename = '051917 GC optimizing EPSP'
if len(sys.argv) > 3:
    cluster_id = sys.argv[3]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

dv = c[:]
dv.block = True
global_start_time = time.time()
dv.execute('run parallel_clustered_branch_cooperativity_nmda_engine_GC %i \"%s\"' % (int(spines), mech_filename))
#time.sleep(300)
v = c.load_balanced_view()

start_time = time.time()

instructions = []
for i in range(len(parallel_clustered_branch_cooperativity_nmda_engine_GC.groups_to_stim)):
    for j in range(len(parallel_clustered_branch_cooperativity_nmda_engine_GC.groups_to_stim[i]['synapses'])):
        instructions.append((i, j))
result = v.map_async(parallel_clustered_branch_cooperativity_nmda_engine_GC.stim_single_expected, instructions)
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
for i in range(len(parallel_clustered_branch_cooperativity_nmda_engine_GC.groups_to_stim)):
    for j in range(1, len(parallel_clustered_branch_cooperativity_nmda_engine_GC.groups_to_stim[i]['synapses'])+1):
        instructions.append((i, j))
result = v.map_async(parallel_clustered_branch_cooperativity_nmda_engine_GC.stim_actual_group, instructions)
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
print 'Parallel simulation took %i s to stimulate %i groups of synapses' % (time.time() - start_time,
                                        len(parallel_clustered_branch_cooperativity_nmda_engine_GC.groups_to_stim))