from mpi4py import MPI
from function_lib import *
from collections import Counter

comm = MPI.COMM_WORLD
rank = comm.rank

# projections_dir = morph_dir
projections_dir = os.environ['PI_HOME']+'/'
projections_file = 'DGC_forest_connectivity.h5'
degree_file = 'DGC_forest_connectivity_degrees_new.h5'

f = h5py.File(projections_dir+projections_file, 'r', driver='mpio', comm=comm)


def pulse_read_iterator(dset, total_size):
    """

    :param dset: :class:'h5py.Dataset'
    :param total_size: int
    :return: tuple of list
    """
    count = 0
    piece_size = comm.size
    global_start = 0
    global_end = min(total_size, piece_size)
    # if rank == 0:
    #    print 'total_size: %i' % total_size
    while global_start < total_size:
        #if rank == 0:
        #    print 'global_start: %i' % global_start
        gid_list = []
        start_index_list = []
        end_index_list = []
        local_start = global_start
        if global_end > total_size:
            global_end = total_size
        while local_start < global_end:
            # if rank == 0:
            #    print 'local_start: %i' % local_start
            this_gid = dset['gid'][local_start]
            this_start_index = dset['ptr'][local_start]
            this_end_index = dset['ptr'][local_start + 1]
            gid_list.append(this_gid)
            start_index_list.append(this_start_index)
            end_index_list.append(this_end_index)
            local_start += 1
        while len(gid_list) < comm.size:
            gid_list.append(this_gid)
            start_index_list.append(global_end - 1)
            end_index_list.append(global_end - 1)
        if rank == 0:
            print 'Scattering chunk %i' % count
        count += 1
        global_start += piece_size
        global_end += piece_size
        yield gid_list, start_index_list, end_index_list


chunk = pulse_read_iterator(f['Populations']['GC']['Connectivity']['source_gid'],
                            len(f['Populations']['GC']['Connectivity']['source_gid']['gid']))

if rank == 0:
    in_degree = Counter()
    out_degree = Counter()
else:
    in_degree, out_degree = None, None

global_start_time = time.time()
for chunk_count, (gid_list, start_index_list, end_index_list) in enumerate(chunk):
# for chunk_count, (gid_list, start_index_list, end_index_list) in [enumerate(chunk).next()]:
    start_time = time.time()
    if rank == 0:
        print 'Chunk %i: gid_list: %s, num edges: %i' % (chunk_count, str(gid_list),
                                                         np.sum(np.subtract(end_index_list, start_index_list)))

    this_gid = comm.scatter(gid_list, root=0)
    this_start_index = comm.scatter(start_index_list, root=0)
    this_end_index = comm.scatter(end_index_list, root=0)
    this_piece_size = this_end_index - this_start_index

    if this_piece_size > 0:
        source_gids = f['Populations']['GC']['Connectivity']['source_gid']['value'][this_start_index:this_end_index]
        print 'Rank %i received %i edges' % (rank, len(source_gids))
        print 'target_gid: %i, source_gid: %i' % (this_gid, source_gids[0])
    else:
        source_gids = np.array([], dtype='uint32')

    local_in_degree = Counter()
    local_out_degree = Counter()

    local_in_degree.update({this_gid: this_piece_size})
    for source_gid in source_gids:
        local_out_degree.update({source_gid: 1})
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
