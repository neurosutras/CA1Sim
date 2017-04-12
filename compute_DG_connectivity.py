from function_lib import *
from mpi4py import MPI
from neurotrees.io import append_cell_attributes
from neurotrees.io import NeurotreeGen
from neurotrees.io import NeurotreeAttrGen
from neurotrees.io import bcast_cell_attributes
from neurotrees.io import population_ranges
# import mkl
import sys
import os
import gc

"""
Determine synaptic connectivity onto DG GCs based on target convergences, divergences, and axonal distances.

Algorithm:
1. For each cell population:
    i. Load the soma locations of each cell population from an .hdf5 file. The (U, V) coordinates of each cell will be
        projected onto a plane in the middle of the granule cell layer (U', V') (L = -1), and used to calculate the
        orthogonal arc distances (S-T and M-L) between the projected soma locations.
2. For each cell, for each type of connection:
    i. Compute a probability of connection across all possible sources, based on the estimated arc distances between
        their projected soma locations.
    ii. Load from a neurotree file the synapses metadata, including layer, type, syn_loc, sec_type, and unique indexes 
        for each synapse.
    ii. Write to a neurotree file the source_gids and synapse_indexes for all the connections that have been
        specified in this step. Iterating through connection types will keep appending to this data structure.
    TODO: Implement a parallel write method to the separate connection edge data structure, use that instead.

swc_type_enumerator = {'soma': 1, 'axon': 2, 'basal': 3, 'apical': 4, 'trunk': 5, 'tuft': 6}
syn_type_enumerator = {'excitatory': 0, 'inhibitory': 1, 'neuromodulatory': 2}

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
forest_file = 'DGC_forest_connectivity_040617.h5'

# synapse_dict = read_from_pkl(neurotrees_dir+'010117_GC_test_synapse_attrs.pkl')
#synapse_dict = read_tree_attributes(MPI._addressof(comm), neurotrees_dir+forest_file, 'GC',
#                                    namespace='Synapse_Attributes')

coords_dir = morph_dir
# coords_dir = os.environ['PI_SCRATCH']+'/DG/'
# coords_dir = os.environ['PI_HOME']+'/Full_Scale_Control/'
coords_file = 'dentate_Full_Scale_Control_coords.h5'

start_time = time.time()
last_time = start_time

soma_coords = {}
for population in population_ranges(MPI._addressof(comm), coords_dir+coords_file).iterkeys():
    soma_coords[population] = bcast_cell_attributes(MPI._addressof(comm), 0, coords_dir+coords_file, population,
                                                    namespace='Coordinates')

spatial_resolution = 1.  # um
max_u = 10750.
max_v = 2957.

du = (0.98*np.pi-0.01*np.pi)/max_u*spatial_resolution
dv = (1.425*np.pi-(-0.23*np.pi))/max_v*spatial_resolution
u = np.arange(0.01*np.pi-10.*du, 0.98*np.pi+10.*du, du)
v = np.arange(-0.23*np.pi-10*dv, 1.425*np.pi+10.*dv, dv)

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

# full width in um (S-T, M-L)
axon_width = {'GC': (900., 900.), 'MPP': (1500., 3000.), 'LPP': (1500., 3000.), 'MC': (4000., 4000.)}
axon_offset = {'MC': (1000., 0.)}


layers = {'GC': {'MPP': [2], 'LPP': [3], 'MC': [1]}}
syn_types = {'GC': {'MPP': [0], 'LPP': [0], 'MC': [0]}}
# sec_types = {'GC': {'MPP': [4], 'LPP': [4], 'MC': [4]}}
# Original synapse attributes mistakenly used 5 for apical swc_type
swc_types = {'GC': {'MPP': [5], 'LPP': [5], 'MC': [5]}}
# fraction of synapses matching this layer, syn_type, sec_type, and source
proportions = {'GC': {'MPP': [1.], 'LPP': [1.], 'MC': [1.]}}

get_array_index = np.vectorize(lambda val_array, this_val: np.where(val_array >= this_val)[0][0], excluded=[0])

for population in soma_coords:
    for cell in soma_coords[population].itervalues():
        cell['u_index'] = get_array_index(u, cell['U Coordinate'][0])
        cell['v_index'] = get_array_index(v, cell['V Coordinate'][0])


def filter_sources(target, layer, swc_type, syn_type):
    """
    
    :param target: str 
    :param layer: int
    :param swc_type: int
    :param syn_type: int
    :return: list
    """
    source_list = []
    proportion_list = []
    for source in layers[target]:
        for i, this_layer in enumerate(layers[target][source]):
            if this_layer == layer:
                if swc_types[target][source][i] == swc_type:
                    if syn_types[target][source][i] == syn_type:
                        source_list.append(source)
                        proportion_list.append(proportions[target][source][i])
    return source_list, proportion_list


def filter_synapses(synapse_dict, layer, swc_type, syn_type):
    """
    
    :param synapse_dict: 
    :param layer: int 
    :param swc_type: int
    :param syn_type: int
    :return: 
    """
    return np.where((synapse_dict['layer'] == layer) & (synapse_dict['swc_type'] == swc_type)
                    & (synapse_dict['syn_type'] == syn_type))[0]


class AxonProb(object):
    """
    An object of this class will instantiate customized vectorized functions describing the connection probabilities
    for each presynaptic population. These functions can then be used to get the distribution of connection 
    probabilities across all possible source neurons, given the soma coordinates of a target neuron. Heavy on
    approximations, but is fast to compute, and independent for each synapse, so does not require any sampling without
    replacement, MPI communication, or produce any undesirable order or edge effects.
    """
    def __init__(self, axon_width, axon_offset):
        """
        Warning: This method does not produce an absolute probability. It must be normalized so that the total area
        (volume) under the distribution is 1 before sampling.
        :param axon_width: dict: {source: (tuple of float)}
        :param axon_offset: dict: {source: (tuple of float)}
        """
        self.p_dist = {}
        self.width = {}
        self.offset = {}
        self.sigma = {}
        for source in axon_width:
            self.width[source] = {'u': axon_width[source][0], 'v': axon_width[source][1]}
            self.sigma[source] = {axis: self.width[source][axis] / 3. / np.sqrt(2.) for axis in self.width[source]}
            if source in axon_offset:
                self.offset[source] = {'u': axon_offset[source][0], 'v': axon_offset[source][1]}
            else:
                self.offset[source] = {'u': 0., 'v': 0.}
            self.p_dist[source] = (lambda source: np.vectorize(lambda distance_u, distance_v:
                                               np.exp(-(((abs(distance_u) - self.offset[source]['u']) /
                                                         self.sigma[source]['u'])**2. +
                                                        ((abs(distance_v) - self.offset[source]['v']) /
                                                         self.sigma[source]['v'])**2.))))(source)

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

    def filter_by_soma_coords(self, target, source, target_gid, soma_coords, distance_U, distance_V):
        """
        Given the coordinates of a target neuron, filter the set of source neurons, and return the arc_distances in two
        dimensions and the gids of source neurons whose axons potentially contact the target neuron.
        :param target: str
        :param source: str
        :param target_gid: int
        :param soma_coords: nested dict of array
        :param distance_U: array of float
        :param distance_V: array of float
        :return: tuple of array of int
        """
        target_index_u = soma_coords[target][target_gid]['u_index']
        target_index_v = soma_coords[target][target_gid]['v_index']
        source_distance_u = []
        source_distance_v = []
        source_gid = []
        for this_source_gid in soma_coords[source]:
            this_source_distance_u, this_source_distance_v = \
                self.get_approximate_arc_distances(target_index_u, target_index_v,
                                                   soma_coords[source][this_source_gid]['u_index'],
                                                   soma_coords[source][this_source_gid]['v_index'], distance_U,
                                                   distance_V)
            if ((np.abs(this_source_distance_u) <= self.width[source]['u'] / 2. + self.offset[source]['u']) &
                    (np.abs(this_source_distance_v) <= self.width[source]['v'] / 2. + self.offset[source]['v'])):
                source_distance_u.append(this_source_distance_u)
                source_distance_v.append(this_source_distance_v)
                source_gid.append(this_source_gid)

        return np.array(source_distance_u), np.array(source_distance_v), np.array(source_gid)

    def get_p(self, target, source, target_gid, soma_coords, distance_U, distance_V, plot=False):
        """
        Given the soma coordinates of a target neuron and a population source, return an array of connection 
        probabilities and an array of corresponding source gids.
        :param target: str
        :param source: str
        :param target_gid: int
        :param soma_coords: nested dict of array
        :param distance_U: array of float
        :param distance_V: array of float
        :param plot: bool
        :return: array of float, array of int
        """
        source_distance_u, source_distance_v, source_gid = self.filter_by_soma_coords(target, source, target_gid,
                                                                                      soma_coords, distance_U,
                                                                                      distance_V)
        p = self.p_dist[source](source_distance_u, source_distance_v)
        p /= np.sum(p)
        if plot:
            plt.scatter(source_distance_u, source_distance_v, c=p)
            plt.title(source+' -> '+target)
            plt.xlabel('Septotemporal distance (um)')
            plt.ylabel('Tranverse distance (um)')
            plt.show()
            plt.close()
        return p, source_gid

print 'Rank %i took %i s to sort cells' % (rank, time.time() - last_time)
last_time = time.time()

p_connect = AxonProb(axon_width, axon_offset)
print 'Rank %i took %i s to initialize connection generator' % (rank, time.time() - last_time)
last_time = time.time()

local_np_random = np.random.RandomState()
connectivity_seed_offset = 100000000  # make sure random seeds are not being reused for various types of
                                      # stochastic sampling

connection_dict = {}
p_dict = {}
source_gid_dict = {}

target = 'GC'
count = 0
for target_gid, attributes_dict in NeurotreeAttrGen(MPI._addressof(comm), neurotrees_dir+forest_file, target,
                                                    io_size=comm.size, namespace='Synapse_Attributes'):
    count += 1
    synapse_dict = attributes_dict['Synapse_Attributes']
    local_np_random.seed(target_gid + connectivity_seed_offset)
    connection_dict['source_gid'] = np.array([], dtype='uint32')
    connection_dict['syn_id'] = np.array([], dtype='uint32')

    layer_set, swc_type_set, syn_type_set = set(), set(), set()
    for source in layers[target].keys():
        layer_set.update(layers[target][source])
        swc_type_set.update(swc_types[target][source])
        syn_type_set.update(syn_types[target][source])

    for layer in layer_set:
        for swc_type in swc_type_set:
            for syn_type in syn_type_set:
                sources, this_proportions = filter_sources(target, layer, swc_type, syn_type)
                if sources:
                    print 'Connections in layer %i (swc_type: %i, syn_type: %i): %s' % (layer, swc_type, syn_type,
                                                                                        '[' + ', '.join(
                                                                                            ['%s' % xi for xi in
                                                                                             sources]) + ']')
                    p, source_gid = np.array([]), np.array([])
                    for source, this_proportion in zip(sources, this_proportions):
                        if source not in source_gid_dict:
                            this_p, this_source_gid = p_connect.get_p(target, source, target_gid, soma_coords, distance_U,
                                                                  distance_V, plot=True)
                            source_gid_dict[source] = this_source_gid
                            p_dict[source] = this_p
                        else:
                            this_source_gid = source_gid_dict[source]
                            this_p = p_dict[source]
                        p = np.append(p, this_p * this_proportion)
                        source_gid = np.append(source_gid, this_source_gid)
                    syn_indexes = filter_synapses(synapse_dict, layer, swc_type, syn_type)
                    connection_dict['syn_id'] = np.append(connection_dict['syn_id'],
                                                          synapse_dict['syn_id'][syn_indexes]).astype('uint32',
                                                                                                      copy=False)
                    this_source_gid = local_np_random.choice(source_gid, len(syn_indexes), p=p)
                    connection_dict['source_gid'] = np.append(connection_dict['source_gid'],
                                                              this_source_gid).astype('uint32', copy=False)
    if count > 0:
        break

print 'Rank %i took %i s to compute connectivity for target: %s, gid: %i' % (rank, time.time() - last_time, target,
                                                                             target_gid)
"""

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