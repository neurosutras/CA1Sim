from mpi4py import MPI
from function_lib import *
from collections import Counter

comm = MPI.COMM_WORLD
rank = comm.rank

projections_dir = morph_dir
# projections_dir = os.environ['PI_HOME']+'/'
projections_file = '100716_dentate_MPPtoDGC.h5'
degree_file = 'DGC_forest_connectivity_degrees_orig.h5'

f = h5py.File(projections_dir+projections_file, 'r', driver='mpio', comm=comm)


def pulse_read_iterator(chunk_size, total_size):
    count = 0
    piece_size = int(np.ceil(float(chunk_size)/comm.size))
    global_start = 0
    global_end = min(total_size, chunk_size)
    # if rank == 0:
    #    print 'total_size: %i' % total_size
    while global_start < total_size:
        #if rank == 0:
        #    print 'global_start: %i' % global_start
        pieces_indexes = []
        pieces_sizes = []
        local_start = global_start
        if global_end > total_size:
            global_end = total_size
        while local_start < global_end:
            # if rank == 0:
            #    print 'local_start: %i' % local_start
            pieces_indexes.append(local_start)
            this_piece_size = min(piece_size, global_end-local_start)
            pieces_sizes.append(this_piece_size)
            local_start += piece_size
        while len(pieces_indexes) < comm.size:
            pieces_indexes.append(global_end-1)
            pieces_sizes.append(0)
        if rank == 0:
            print 'Scattering chunk %i' % count
        count += 1
        global_start += chunk_size
        global_end += chunk_size
        yield pieces_indexes, pieces_sizes

# chunk = pulse_read_iterator(100000000, len(f['Projections']['MPPtoGC']['Source']))
chunk = pulse_read_iterator(1000000, int(len(f['Projections']['MPPtoGC']['Source'])/100))

if rank == 0:
    in_degree = Counter()
    out_degree = Counter()
else:
    in_degree, out_degree = None, None

global_start_time = time.time()
for chunk_count, (pieces_indexes, pieces_sizes) in enumerate(chunk):
# for chunk_count, (pieces_indexes, pieces_sizes) in [enumerate(chunk).next()]:
    start_time = time.time()
    if rank == 0:
        print 'Chunk %i: Indexes: %s, Sizes: %s' % (chunk_count, str(pieces_indexes), str(pieces_sizes))

    this_piece_index = comm.scatter(pieces_indexes, root=0)
    this_piece_size = comm.scatter(pieces_sizes, root=0)
    """
    dset = f['Projections']['MPPtoGC']['Source']
    with dset.collective:
        sources = dset[this_piece_index:this_piece_index+this_piece_size]
    dset = f['Projections']['MPPtoGC']['Destination']
    with dset.collective:
        destinations = dset[this_piece_index:this_piece_index+this_piece_size]
    if this_piece_size > 0:
        print 'Rank %i received %i data points' % (rank, len(sources))
        print 'Index: %i, Destination: %i, Source: %i' % (this_piece_index, destinations[0], sources[0])
    """
    if this_piece_size > 0:
        sources = f['Projections']['MPPtoGC']['Source'][this_piece_index:this_piece_index+this_piece_size]
        destinations = f['Projections']['MPPtoGC']['Destination'][this_piece_index:this_piece_index + this_piece_size]
        print 'Rank %i received %i data points' % (rank, len(sources))
        print 'Index: %i, Destination: %i, Source: %i' % (this_piece_index, destinations[0], sources[0])
    else:
        sources = np.array([], dtype='uint32')
        destinations = np.array([], dtype='uint32')

    local_in_degree = Counter()
    local_out_degree = Counter()

    for i in range(len(destinations)):
        local_in_degree.update({destinations[i]: 1})
        local_out_degree.update({sources[i]: 1})
    local_in_degree = dict(local_in_degree)
    local_out_degree = dict(local_out_degree)

    in_degree_piece_list = comm.gather(local_in_degree, root=0)
    out_degree_piece_list = comm.gather(local_out_degree, root=0)

    if rank == 0:
        for in_degree_piece in in_degree_piece_list:
            in_degree.update(in_degree_piece)
        for out_degree_piece in out_degree_piece_list:
            out_degree.update(out_degree_piece)
        print 'Processed chunk %i in %i s (Elapsed time: %i s)' % (chunk_count, time.time() - start_time,
                                                                   time.time() - global_start_time)

comm.barrier()
f.close()

if rank == 0:
    print 'Processed %i edges in %i s' % (np.sum(in_degree.values()), time.time() - global_start_time)
    with h5py.File(projections_dir+degree_file, 'w') as f:
        f.create_group('Projections')
        f['Projections'].create_group('MPPtoGC')
        f['Projections']['MPPtoGC'].create_group('In degree')
        f['Projections']['MPPtoGC']['In degree'].create_dataset('gid', compression='gzip', compression_opts=9,
                                                   data=np.array(in_degree.keys(), dtype='uint32'))
        f['Projections']['MPPtoGC']['In degree'].create_dataset('count', compression='gzip', compression_opts=9,
                                                                data=np.array(in_degree.values(), dtype='uint32'))
        f['Projections']['MPPtoGC'].create_group('Out degree')
        f['Projections']['MPPtoGC']['Out degree'].create_dataset('gid', compression='gzip', compression_opts=9,
                                                                data=np.array(out_degree.keys(), dtype='uint32'))
        f['Projections']['MPPtoGC']['Out degree'].create_dataset('count', compression='gzip', compression_opts=9,
                                                                data=np.array(out_degree.values(), dtype='uint32'))
