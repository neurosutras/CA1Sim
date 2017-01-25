from function_lib import *
from mpi4py import MPI
from neurotrees.io import append_tree_attributes
from neurotrees.io import read_tree_attributes
# import mkl
import sys
import os
import gc

"""
Determine MPP-DG_GC synaptic connectivity based on target convergences, divergences, and axonal distances.

Algorithm:
1. For each cell population:
    i. Load the soma locations of each cell population from an .hdf5 file. The (U, V) coordinates of each cell will be
        projected onto a plane in the middle of the granule cell layer (U', V') (L = -1), and used to calculate the
        orthogonal arc distances (S-T and M-L) between the projected soma locations.
2. For each cell, for each type of connection:
    i. Compute a probability of connection across all possible sources, based on the estimated arc distances between
        their projected soma locations.
    ii. Load from a neurotree file the synapses metadata, including layer, type, syn_loc, and unique indexes for each
        synapse.
    ii. Write to a neurotree file a the source_gids and synapse_indexes for all the connections that have been
        specified in this step. Iterating through connection types will keep appending to this data structure.
    TODO: Implement a parallel write method to the separate connection edge data structure, use that instead.

"""

# mkl.set_num_threads(2)

comm = MPI.COMM_WORLD
rank = comm.rank  # The process ID (integer 0-3 for 4-process run)

if rank == 0:
    print '%i ranks have been allocated' % comm.size
sys.stdout.flush()

neurotrees_dir = morph_dir
# neurotrees_dir = os.environ['PI_SCRATCH']+'/DGC_forest/hdf5/'
# neurotrees_dir = os.environ['PI_HOME']+'/'
# forest_file = '122016_DGC_forest_with_syn_locs.h5'
forest_file = 'DGC_forest_connectivity_test.h5'

# synapse_dict = read_from_pkl(neurotrees_dir+'010117_GC_test_synapse_attrs.pkl')
synapse_dict = read_tree_attributes(MPI._addressof(comm), neurotrees_dir+forest_file, 'GC',
                                    namespace='Synapse_Attributes')

coords_dir = morph_dir
# coords_dir = os.environ['PI_SCRATCH']+'/DG/'
# coords_dir = os.environ['PI_HOME']+'/'
coords_file = 'DG_Soma_Coordinates.h5'

target_GID = synapse_dict.keys()
target_GID.sort()

print 'MPI rank %i received %i GCs: [%i:%i]' % (rank, len(target_GID), target_GID[0], target_GID[-1])
sys.stdout.flush()

soma_coords = {}
with h5py.File(coords_dir+coords_file, 'r') as f:
    for import_pop_key, local_pop_key in zip (['Granule Cells', 'MEC Cells'], ['GC', 'MEC']):
        soma_coords[local_pop_key] = {'u': f['Coordinates'][import_pop_key]['U Coordinate'][:],
                                      'v': f['Coordinates'][import_pop_key]['V Coordinate'][:],
                                      'GID': f['Coordinates'][import_pop_key]['GID'][:]}


spatial_resolution = 1.  # um
max_u = 10826.
max_v = 3002.

du = (0.98*np.pi-0.01*np.pi)/max_u*spatial_resolution
dv = (1.425*np.pi-(-0.23*np.pi))/max_v*spatial_resolution
u = np.arange(0.01*np.pi, 0.98*np.pi+du, du)
v = np.arange(-0.23*np.pi, 1.425*np.pi+dv, dv)

U, V = np.meshgrid(u, v, indexing='ij')

# for the middle of the granule cell layer:
L = -1.
X = np.array(-500.* np.cos(U) * (5.3 - np.sin(U) + (1. + 0.138 * L) * np.cos(V)))
Y = np.array(750. * np.sin(U) * (5.5 - 2. * np.sin(U) + (0.9 + 0.114*L) * np.cos(V)))
Z = np.array(2500. * np.sin(U) + (663. + 114. * L) * np.sin(V - 0.13 * (np.pi-U)))

euc_coords = np.array([X.T, Y.T, Z.T]).T

delta_U = np.sqrt((np.diff(euc_coords, axis=0)**2.).sum(axis=2))
delta_V = np.sqrt((np.diff(euc_coords, axis=1)**2.).sum(axis=2))

distance_U = np.cumsum(np.insert(delta_U, 0, 0., axis=0), axis=0)
distance_V = np.cumsum(np.insert(delta_V, 0, 0., axis=1), axis=1)

width_U = np.mean(np.max(distance_U, axis=0))
width_V = np.mean(np.max(distance_V, axis=1))

# MEC layer II stellate cells: Gatome et al., Neuroscience, 2010
pop_size = {'GC': 1000000, 'MEC': 38000}

pop_density = {pop: pop_size[pop] / width_U / width_V for pop in pop_size}  # cell density per um^2
# Sloviter, Lomo, Frontiers, 2012. says GCs might only project 200 um S-T..., MEC less than 1500.
axon_width = {'GC': (900., 900.), 'MEC': (1500., 1500.)}  # full width in um (S-T, M-L)

# based on EM synapse density and dendritic length in MML
convergence = {'GC': {'MEC': 3000}}  # {target: {source: # of synapses from source cell population onto 1 target cell}}
# can't find good estimate of MEC divergence. circular reasoning: convergence onto GC * # of GCs / # of MEC cells
divergence = {'MEC': {'GC': 79000}}  # {source: {target: # of synapses from 1 source cell to target cell population}}
layers = {'GC': {'MEC': [2]}}
syn_types = {'GC': {'MEC': 0}}


get_array_index = np.vectorize(lambda val_array, this_val: np.where(val_array >= this_val)[0][0], excluded=[0])


class AxonProb(object):
    """
    An object of this class will instantiate customized vectorized functions describing the connection probabilities
    for each connection specified in 'divergence'. Those functions can then be used to get the distribution of
    connection probabilities across all possible source neurons, given the soma coordinates of a target neuron. Heavy on
    approximations, but is fast to compute, and independent for each synapse, so does not require any sampling without
    replacement, MPI communication, or produce any undesirable order or edge effects.
    """
    def __init__(self, divergence, convergence, axon_width, pop_density):
        """

        :param divergence: dict: {source: {target: int}}
        :param convergence: dict: {target: {source: int}}
        :param axon_width: dict: {source: (tuple of float)}
        :param pop_density: dict: {target: float}
        """
        self.p_dist = {}
        for source in divergence:
            self.p_dist[source] = {}
            width = {'u': axon_width[source][0], 'v': axon_width[source][1]}
            sigma = {axis: width[axis] / 3. / np.sqrt(2.) for axis in width}
            for target in divergence[source]:
                surface_area = np.pi * width['u'] * width['v'] / 4.
                total_p = divergence[source][target] / convergence[target][source] / surface_area / pop_density[target]
                this_volume = np.pi * sigma['u'] * sigma['v']
                self.p_dist[source][target] = np.vectorize(lambda distance_u, distance_v:
                                                   total_p * np.exp(-((distance_u / sigma['u'])**2. +
                                                                      (distance_v / sigma['v'])**2.)) /
                                                   this_volume)

    def get_approximate_arc_distances(self, target_index_u, target_index_v, source_indexes_u, source_indexes_v,
                                      distance_U, distance_V):
        """
        Arc distances along 2 basis dimensions are calculated as the average of the arc distances along parallel edges
        of a parallelogram with the soma locations of the pair of neurons as vertices.
        :param target_index_u: int
        :param target_index_v: int
        :param source_indexes_u: array of int
        :param source_indexes_v: array of int
        :param distance_U: array of float
        :param distance_V: array of float
        :return: tuple of array of float
        """
        distance_u0 = np.subtract(distance_U[source_indexes_u, target_index_v],
                                  distance_U[target_index_u, target_index_v])
        distance_u1 = np.subtract(distance_U[source_indexes_u, source_indexes_v],
                                  distance_U[target_index_u, source_indexes_v])
        distance_u = np.mean(np.array([distance_u0, distance_u1]), axis=0)
        distance_v0 = np.subtract(distance_V[target_index_u, source_indexes_v],
                                  distance_V[target_index_u, target_index_v])
        distance_v1 = np.subtract(distance_V[source_indexes_u, source_indexes_v],
                                  distance_V[source_indexes_u, target_index_v])
        distance_v = np.mean(np.array([distance_v0, distance_v1]), axis=0)

        return distance_u, distance_v

    def filter_soma_coords(self, target_index_u, target_index_v, source, soma_coords, axon_width, u, v, distance_U,
                           distance_V):
        """
        Given the coordinates of a target neuron, filter the set of source neurons, and return the arc_distances in two
        dimensions and the gids of source neurons whose axons potentially contact the target neuron.
        :param target_index_u: int
        :param target_index_v: int
        :param source: str
        :param soma_coords: nested dict of array
        :param axon_width: dict: {source: (tuple of float)}
        :param u: array of float
        :param v: array of float
        :param distance_U: array of float
        :param distance_V: array of float
        :return: tuple of array of int
        """
        source_indexes_u_all = get_array_index(u, soma_coords[source]['u'])
        source_indexes_v_all = get_array_index(v, soma_coords[source]['v'])
        source_distance_u_all, source_distance_v_all = self.get_approximate_arc_distances(target_index_u,
                                                                                          target_index_v,
                                                                                          source_indexes_u_all,
                                                                                          source_indexes_v_all,
                                                                                          distance_U, distance_V)
        width = {'u': axon_width[source][0], 'v': axon_width[source][1]}
        source_indexes = np.where((np.abs(source_distance_u_all) <= width['u'] / 2.) &
                                  (np.abs(source_distance_v_all) <= width['v'] / 2.))[0]
        return source_distance_u_all[source_indexes], source_distance_v_all[source_indexes], \
               soma_coords[source]['GID'][source_indexes]

    def choose_sources(self, target, source, target_gid, num_syns, soma_coords, axon_width, u, v, distance_U,
                       distance_V, local_random=None):
        """
        Given the soma coordinates of a target neuron, returns a list of length num_syns containing the gids of source
        neurons based on the specified connection probabilities.
        :param target: str
        :param source: str
        :param target_gid: int
        :param num_syns: int
        :param soma_coords: nested dict of array
        :param axon_width: dict: {source: (tuple of float)}
        :param u: array of float
        :param v: array of float
        :param distance_U: array of float
        :param distance_V: array of float
        :param local_random: :class:'np.random.RandomState'
        :return: list of int
        """
        if local_random is None:
            local_random = np.random.RandomState()
        target_gid_index = np.where(soma_coords[target]['GID'] == target_gid)[0][0]
        target_index_u = np.where(u >= soma_coords[target]['u'][target_gid_index])[0][0]
        target_index_v = np.where(v >= soma_coords[target]['v'][target_gid_index])[0][0]

        source_distance_u, source_distance_v, source_gid = self.filter_soma_coords(target_index_u, target_index_v,
                                                                                   source, soma_coords, axon_width, u,
                                                                                   v, distance_U, distance_V)
        source_p = self.p_dist[source][target](source_distance_u, source_distance_v)
        source_p /= np.sum(source_p)
        return local_random.choice(source_gid, num_syns, p=source_p)


p_connect = AxonProb(divergence, convergence, axon_width, pop_density)

local_np_random = np.random.RandomState()
connectivity_seed_offset = 100000000  # make sure random seeds are not being reused for various types of
                                      # stochastic sampling
start_time = time.time()

connection_dict = {}
target = 'GC'

# block_size = int(7000/comm.size)
block_size = 5
if 'SYN_START_INDEX' in os.environ:
    start_index = int(os.environ['SYN_START_INDEX'])
else:
    start_index = 0
end_index = start_index+block_size

count = 0
while start_index < block_size:
#while start_index < len(target_GID):
    connection_dict = {}
    for target_gid in target_GID[start_index:end_index]:
        this_synapse_dict = synapse_dict[target_gid]
        if 'syn_id' not in this_synapse_dict:
            this_synapse_dict['syn_id'] = np.array(range(len(this_synapse_dict['syn_locs'])))
        connection_dict[target_gid] = {}
        connection_dict[target_gid]['syn_id'] = np.array([], dtype='uint32')
        connection_dict[target_gid]['source_gid'] = np.array([], dtype='uint32')
        for source in convergence[target]:
            target_layers = layers[target][source]
            target_syn_type = syn_types[target][source]
            target_indexes = np.where((np.in1d(this_synapse_dict['layer'], target_layers)) &
                                      (this_synapse_dict['syn_type'] == target_syn_type))[0]
            connection_dict[target_gid]['syn_id'] = \
                np.append(connection_dict[target_gid]['syn_id'],
                          this_synapse_dict['syn_id'][target_indexes]).astype('uint32', copy=False)
            these_source_gids = p_connect.choose_sources(target, source, target_gid, len(target_indexes), soma_coords,
                                                         axon_width, u, v, distance_U, distance_V, local_np_random)
            connection_dict[target_gid]['source_gid'] = np.append(connection_dict[target_gid]['source_gid'],
                                                                  these_source_gids).astype('uint32', copy=False)
        count += 1
    append_tree_attributes(MPI._addressof(comm), neurotrees_dir + forest_file, 'GC', connection_dict,
                           namespace='Connectivity', value_chunk_size=48000)
    if end_index >= len(target_GID):
        last_index = len(target_GID) - 1
    else:
        last_index = end_index - 1
    print 'MPI rank %i wrote to file connectivity for cell gid: %i (count: %i)' % (rank, target_GID[last_index],
                                                                                         count)
    sys.stdout.flush()
    del connection_dict
    gc.collect()
    start_index += block_size
    end_index += block_size


count_fragments = comm.gather(count, root=0)
# connection_dict_fragments = comm.gather(connection_dict, root=0)
if rank == 0:
    print '%i ranks took took %.2f s to compute connectivity for %i cells' % (comm.size, time.time() - start_time,
                                                                              np.sum(count_fragments))
    # connection_dict = {key: value for piece in connection_dict_fragments for key, value in piece.items()}
    # write_to_pkl(neurotrees_dir+'010117_DG_GC_MEC_connectivity_test.pkl', connection_dict)

"""
connection_dict = read_from_pkl(neurotrees_dir+'010117_DG_GC_MEC_connectivity_test.pkl')
target_gid = 500
plot_source_soma_locs(target_gid, 'GC', 'MEC', connection_dict, soma_coords, u, v, U, V, distance_U, distance_V,
                     X, Y, Z)

start_time = time.time()
syn_in_degree = {target: {source: {i: 0 for i in range(len(pop_locs_X[target]))}
                          for source in convergence[target]} for target in convergence}

syn_out_degree = {source: {target: {i: 0 for i in range(len(pop_locs_X[source]))}
                           for target in divergence[source]} for source in divergence}

with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
    for target in f:
        for target_index in f[target]:
            for source in f[target][target_index]:
                syn_in_degree[target][source][int(target_index)] += len(f[target][target_index][source])
                for source_index in f[target][target_index][source]:
                    syn_out_degree[source][target][source_index] += 1

print 'Calculating out degree took %i s' % (time.time() - start_time)


for target in convergence:
    for source in convergence[target]:
        sc = plt.scatter(pop_locs_X[target], pop_locs_Y[target], c=np.array(syn_in_degree[target][source].values()),
                         linewidths=0)
        plt.xlabel('Location S-T (um)')
        plt.ylabel('Location M-L (um)')
        plt.title(source+' -> '+target+' Convergence (# of Synapses)')
        plt.colorbar(sc)
        plt.show()
        plt.close()

for source in divergence:
    for target in divergence[source]:
        sc = plt.scatter(pop_locs_X[source], pop_locs_Y[source], c=np.array(syn_out_degree[source][target].values()),
                         linewidths=0)
        plt.xlabel('Location S-T (um)')
        plt.ylabel('Location M-L (um)')
        plt.title(source + ' -> ' + target + ' Divergence (# of Synapses)')
        plt.colorbar(sc)
        plt.show()
        plt.close()

"""