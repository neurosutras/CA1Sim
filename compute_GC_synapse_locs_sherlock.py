from mpi4py import MPI
from specify_cells2 import *
from neurotrees.io import read_trees
import h5py
# from neurotrees.io import append_tree_attributes
import mkl
import sys
import os
import gc
        
mkl.set_num_threads(1)

comm = MPI.COMM_WORLD
rank = comm.rank

if rank == 0:
    print '%i ranks have been allocated' % comm.size

# neurotrees_dir = morph_dir
# forest_file = '122016_DGC_forest_test_copy.h5'
neurotrees_dir = os.environ['PI_SCRATCH']+'/DGC_forest/hdf5/'
forest_file = 'DGC_forest_full_h5py.h5'
# forest_file = 'DGC_forest_test.h5'
# forest_file = 'DGC_forest_test_010917.h5'
synapse_file = 'DGC_forest_full_synapse_attributes_h5py.h5'
global_start_time = time.time()
start_time = time.time()
morph_dict = read_trees(MPI._addressof(comm), neurotrees_dir+forest_file, 'GC')

GID = morph_dict.keys()
GID.sort()
# GID = GID[:int(1000/comm.size)]

if rank == 0:
    print 'Rank %d: Neurotrees read took %.2f s' % (rank, time.time() - start_time)
print 'Rank %d: received %i GCs: [%i:%i]' % (rank, len(GID), GID[0], GID[-1])
sys.stdout.flush()

synapse_dict = {}
# mismatched_section_dict = {}

# start_time = time.time()
#block_size = int(6000/comm.size)
block_size = 5
if 'SYN_START_INDEX' in os.environ:
    start_index = int(os.environ['SYN_START_INDEX'])
else:
    start_index = 0
end_index = start_index+block_size
count = 0

while start_index < len(GID):
#while start_index < final_index:
    synapse_dict = {}
    for gid in GID[start_index:end_index]:
        print 'Rank %d: processing gid %i' % (rank, gid)
        # start_time = time.time()
        cell = DG_GC(neurotree_dict=morph_dict[gid], gid=gid, full_spines=False)
        # print 'Rank %d: building cell took %.2f s' % (rank, time.time() - start_time)
        """
        this_mismatched_sections = cell.get_mismatched_neurotree_sections()
        if this_mismatched_sections is not None:
            mismatched_section_dict[gid] = this_mismatched_sections
        start_time = time.time()
        print 'Rank %d: checking for mismatched section ids took %.2f s' % (rank, time.time() - start_time)
        start_time = time.time()
        """
        synapse_dict[gid] = cell.export_neurotree_synapse_attributes()
        del cell
        gc.collect()
        # print 'Rank %d: exporting synapse attributes to dict took %.2f s' % (rank, time.time() - start_time)
        sys.stdout.flush()
        count += 1
    start_time = time.time()
    gid_block = comm.allgather(synapse_dict.keys())
    flattened_gid_block = [item for sublist in gid_block for item in sublist]
    groups = synapse_dict.itervalues().next().keys()
    with h5py.File(neurotrees_dir+synapse_file, 'a', driver='mpio', comm=comm) as f:
        for gid in flattened_gid_block:
            f.create_group(str(gid))
            if gid in synapse_dict:
                this_rank = True
            else:
                this_rank = False
            for group in groups:
                if this_rank:
                    dataset_size = len(synapse_dict[gid][group])
                else:
                    dataset_size = 0
                dataset_size = max(comm.allgather(dataset_size))
                f[str(gid)].create_dataset(group, (dataset_size,))
                if this_rank:
                    f[str(gid)][group][0:] = synapse_dict[gid][group]
    if end_index >= len(GID):
        last_index = len(GID)-1
    else:
        last_index = end_index-1
    print 'Rank %d: wrote synapse attributes to file for GCs: [%i:%i]' % (rank, GID[start_index], GID[last_index])
    if rank == 0:
        print 'Writing synapse attributes to file took %.2f s' % (time.time() - start_time)
    sys.stdout.flush()
    del synapse_dict
    gc.collect()
    start_index += block_size
    end_index += block_size

# len_mismatched_section_dict_fragments = comm.gather(len(mismatched_section_dict), root=0)
# len_GID_fragments = comm.gather(len(GID), root=0)
count_fragments = comm.gather(count, root=0)
if rank == 0:
    print '%i ranks took %i s to process synapse attributes for %i morphologies' % (comm.size,
                                                                                   time.time() - global_start_time,
                                                                                   np.sum(count_fragments))
    # print '%i morphologies have mismatched section indexes' % np.sum(len_mismatched_section_dict_fragments)
sys.stdout.flush()
