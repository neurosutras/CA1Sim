from mpi4py import MPI
from specify_cells2 import *
from neurotrees.io import read_trees
from neurotrees.io import append_tree_attributes
# import mkl
import sys
import os
import gc
        
# mkl.set_num_threads(1)

comm = MPI.COMM_WORLD
rank = comm.rank

if rank == 0:
    print '%i ranks have been allocated' % comm.size
sys.stdout.flush()

neurotrees_dir = morph_dir
# forest_file = '122016_DGC_forest_test_copy.h5'
# neurotrees_dir = os.environ['PI_SCRATCH']+'/DGC_forest/hdf5/'
# neurotrees_dir = os.environ['PI_HOME']+'/'
# forest_file = 'DGC_forest_full.h5'
# forest_file = 'DGC_forest_syns_012717.h5'
forest_file = 'DGC_forest_syns_test_012717.h5'

# forest_file = 'DGC_forest_test.h5'

morph_dict = read_trees(MPI._addressof(comm), neurotrees_dir+forest_file, 'GC')[0]

GID = morph_dict.keys()
GID.sort()

print 'MPI rank %d received %i GCs: [%i:%i]' % (rank, len(GID), GID[0], GID[-1])
sys.stdout.flush()

synapse_dict = {}
mismatched_section_dict = {}

start_time = time.time()
#block_size = int(6000/comm.size)
block_size = 5
start_index = 0
end_index = start_index+block_size
count = 0

while start_index < len(GID):
# while start_index < block_size:
    synapse_dict = {}
    for gid in GID[start_index:end_index]:
        print 'Rank: %d, gid: %i' % (rank, gid)
        cell = DG_GC(neurotree_dict=morph_dict[gid], gid=gid, full_spines=False)
        # this_mismatched_sections = cell.get_mismatched_neurotree_sections()
        # if this_mismatched_sections is not None:
        #    mismatched_section_dict[gid] = this_mismatched_sections
        synapse_dict[gid] = cell.export_neurotree_synapse_attributes()
        del cell
        gc.collect()
        sys.stdout.flush()
        count += 1
    append_tree_attributes(MPI._addressof(comm), neurotrees_dir+forest_file, 'GC', synapse_dict,
                           namespace='Synapse_Attributes', value_chunk_size=48000)
    if end_index >= len(GID):
        last_index = len(GID)-1
    else:
        last_index = end_index-1
    print 'MPI rank %d wrote to file synapse locations for GCs: [%i:%i]' % (rank, GID[start_index], GID[last_index])
    del synapse_dict
    gc.collect()
    sys.stdout.flush()
    start_index += block_size
    end_index += block_size

# len_mismatched_section_dict_fragments = comm.gather(len(mismatched_section_dict), root=0)
# len_GID_fragments = comm.gather(len(GID), root=0)
count_fragments = comm.gather(count, root=0)
if rank == 0:
    print '%i ranks took %i s to compute synapse locations for %i morphologies' % (comm.size,
                                                                                   time.time() - start_time,
                                                                                   np.sum(count_fragments))
    # print '%i morphologies have mismatched section indexes' % np.sum(len_mismatched_section_dict_fragments)
