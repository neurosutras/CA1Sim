from mpi4py import MPI
from specify_cells2 import *
from neurotrees.io import read_trees
from neurotrees.io import write_tree_attributes
import mkl
import sys
import os

mkl.set_num_threads(1)

comm = MPI.COMM_WORLD
rank = comm.rank  # The process ID (integer 0-3 for 4-process run)

# neurotrees_dir = morph_dir
# forest_file = '122016_DGC_forest_test_copy.h5'
neurotrees_dir = os.environ['PI_SCRATCH']+'/DGC_forest/hdf5/'
forest_file = 'DGC_forest.h5'
morph_dict = read_trees(MPI._addressof(comm), neurotrees_dir+forest_file, 'GC')

GID = morph_dict.keys()

print 'MPI rank %d received %i GCs: [%i:%i]' % (rank, len(GID), GID[0], GID[-1])

synapse_dict = {}
mismatched_section_dict = {}

start_time = time.time()
for gid in GID[:1]:
    print 'Rank: %d, gid: %i' % (rank, gid)
    cell = DG_GC(neurotree_dict=morph_dict[gid], gid=gid, full_spines=True)
    this_mismatched_sections = cell.get_mismatched_neurotree_sections()
    if this_mismatched_sections is not None:
        mismatched_section_dict[gid] = this_mismatched_sections
    synapse_dict[gid] = cell.export_neurotree_synapse_attributes()
    sys.stdout.flush()

write_tree_attributes(MPI._addressof(comm), neurotrees_dir+forest_file, 'GC', synapse_dict)

mismatched_section_dict_fragments = comm.gather(mismatched_section_dict, root=0)
if rank == 0:
    mismatched_section_dict = {key: value for piece in mismatched_section_dict_fragments for
                               key, value in piece.items()}
    print '%i nodes took %i s to compute synapse locations for %i morphologies' % (comm.size,
                                                                                   time.time() - start_time, len(GID))
    print '%i morphologies have mismatched section indexes' % len(mismatched_section_dict)