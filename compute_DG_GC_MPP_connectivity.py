from function_lib import *
from mpi4py import MPI
from neurotrees.io import read_trees
from neurotrees.io import write_tree_attributes

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
        specified in this step. Iterating through connection types will keep appending to this edge data structure.

"""

comm = MPI.COMM_WORLD
rank = comm.rank  # The process ID (integer 0-3 for 4-process run)

forest_file = '122016_DGC_forest_test_copy.h5'
synapse_dict = read_synapses(MPI._addressof(comm), morph_dir+forest_file, 'GC')
coords_file = 'dentate_Full_Scale_Control_coords_PP.h5'

target_GID = synapse_dict.keys()

print 'MPI rank %d received %i GCs: [%i:%i]' % (rank, len(target_GID), target_GID[0], target_GID[-1])

soma_coords = {}
with h5py.File(data_dir+coords_file, 'r') as f:
    for import_pop_key, local_pop_key in zip (['Granule Cells', 'MEC Cells'], ['GC', 'MEC']):
        soma_coords[local_pop_key] = {'u': f['Coordinates'][import_pop_key]['U Coordinate'][:],
                                      'v': f['Coordinates'][import_pop_key]['V Coordinate'][:],
                                      'GID': f['Coordinates'][import_pop_key]['GID'][:]}


spatial_resolution = 100.  # um
max_u = 10748.
max_v = 2955.
u = np.arange(0.01*np.pi, 0.98*np.pi, (0.98*np.pi-0.01*np.pi)/max_u*spatial_resolution)
v = np.arange(-0.23*np.pi, 1.425*np.pi, (1.425*np.pi-(-0.23*np.pi))/max_v*spatial_resolution)

U, V = np.meshgrid(u, v, indexing='ij')

# for the middle of the granule cell layer:
L =  -1.
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
convergence = {'GC': {'MEC': 3000}}  # {target: {source: # of source cells that contact 1 target cell}}
# can't find good estimate of MEC divergence. circular reasoning: convergence onto GC * # of GCs / # of MEC cells
divergence = {'MEC': {'GC': 79000}}  # {source: {target: # of target cells contacted by 1 source cell}}


class AxonProb(object):
    """
    An object of this class will instantiate customized vectorized functions describing the connection probabilities
    for each connection specified in 'divergence'. Those functions can then be used to get the distribution of
    connection probabilities across all possible source neurons, given the soma coordinates of a target neuron. Heavy on
    approximations, but is fast to compute, and independent for each synapse, so does not require any sampling without
    replacement, MPI communication, or produce any undesirable order or edge effects.
    """
    def __init__(self, divergence, axon_width, pop_density):
        """

        :param divergence: dict: {source: {target: int}}
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
                total_p = divergence[source][target] / surface_area / pop_density[target]
                this_volume = np.pi * sigma['u'] * sigma['v']
                self.p_dist[source][target] = np.vectorize(lambda distance_u, distance_v:
                                                   total_p * np.exp(-((distance_u / sigma['u'])**2. +
                                                                      (distance_v / sigma['v'])**2.)) /
                                                   this_volume)

    def get_approximate_arc_distances(self, target_u_index, target_v_index, source_u_indexes, source_v_indexes,
                                      distance_U, distance_V):
        """
        Arc distances along 2 basis dimensions are calculated as the average of the arc distances along parallel edges
        of a parallelogram with the soma locations of the pair of neurons as vertices.
        :param target_u_index: int
        :param target_v_index: int
        :param source_u_indexes: array of int
        :param source_v_indexes: array of int
        :param distance_U: array of float
        :param distance_V: array of float
        :return: tuple of array of float
        """
        distance_u0 = np.subtract(distance_U[source_u_indexes, target_v_index],
                                  distance_U[target_u_index, target_v_index])
        distance_u1 = np.subtract(distance_U[source_u_indexes, source_v_indexes],
                                  distance_U[target_u_index, source_v_indexes])
        distance_u = np.mean(np.array([distance_u0, distance_u1]), axis=0)
        distance_v0 = np.subtract(distance_V[target_u_index, source_v_indexes],
                                  distance_V[target_u_index, target_v_index])
        distance_v1 = np.subtract(distance_V[source_u_indexes, source_v_indexes],
                                  distance_V[source_u_indexes, target_v_index])
        distance_v = np.mean(np.array([distance_v0, distance_v1]), axis=0)

        return distance_u, distance_v

    def choose_sources(self, target, source, target_gid, num_syns, soma_coords, u, v, distance_U, distance_V,
                       local_random=None):
        """
        Given the soma coordinates of a target neuron, returns a list of length num_syns containing the gids of source
        neurons based on the specified connection probabilities.
        :param target: str
        :param source: str
        :param target_gid: int
        :param num_syns: int
        :param u: array of float
        :param v: array of float
        :param distance_U: array of float
        :param distance_V: array of float
        :param soma_coords: nested dict of array
        :param local_random: :class:'np.random.RandomState'
        :return: list of int
        """
        if local_random is None:
            local_random = np.random.RandomState()
        target_gid_index = np.searchsorted(soma_coords[target]['GID'], target_gid)
        target_u_index = np.where(u == soma_coords[target]['u'][target_gid_index])[0][0]
        target_v_index = np.where(v == soma_coords[target]['v'][target_gid_index])[0][0]
        source_u_indexes = np.searchsorted(u, soma_coords['MEC']['u'])
        source_v_indexes = np.searchsorted(u, soma_coords['MEC']['v'])
        distance_u, distance_v = self.get_approximate_arc_distances(target_u_index, target_v_index, source_u_indexes,
                                                                    source_v_indexes, distance_U, distance_V)
        source_p = self.p_dist[source][target](distance_u, distance_v)
        source_p /= np.sum(source_p)
        return local_random.choice(soma_coords[source]['GID'], num_syns, p=source_p)

p_connect = AxonProb(divergence, axon_width, pop_density)

local_np_random = np.random.RandomState()
connectivity_seed_offset = 100000000  # need to make sure random seeds are not being reused for various types of
                                      # stochastic sampling
start_time = time.time()
with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
    for target in convergence:
        f.create_group(target)
        for j in range(len(pop_locs_X[target])):
            f[target].create_group(str(j))
            target_x = pop_locs_X[target][j]
            target_y = pop_locs_Y[target][j]
            local_np_random.seed(j)
            for source in convergence[target]:
                num_syns = convergence[target][source]
                syn_list = p_connect.choose_sources(target, source, target_x, target_y, num_syns, pop_locs_X,
                                                    pop_locs_Y, local_np_random)
                f[target][str(j)].create_dataset(source, compression='gzip', compression_opts=9, data=syn_list)
                if j % 10000 == 0:
                    print target, j, source, num_syns, time.time() - start_time


# took 45 min
print 'Connectivity algorithm took %i s' % (time.time() - start_time)


start_time = time.time()
syn_in_degree = {target: {source: {i: 0 for i in range(len(pop_locs_X[target]))}
                          for source in convergence[target]} for target in convergence}

syn_out_degree = {source: {target: {i: 0 for i in range(len(pop_locs_X[source]))}
                           for target in divergence[source]} for source in divergence}
"""
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


target = 'GC'
source = 'MEC'
with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
    for target_index in f[target]:
        if (x_width / 2. - 20. <= pop_locs_X[target][int(target_index)] <= x_width / 2. + 20.) and (
                                y_width / 2. - 20. <= pop_locs_Y[target][int(target_index)] <= y_width / 2. + 20.):
            example_target_index = int(target_index)
            example_source_indexes = f[target][target_index][source][:]
            break

H, x_edges, y_edges = np.histogram2d(pop_locs_X[source][example_source_indexes],
                                     pop_locs_Y[source][example_source_indexes], 100)
X_edges, Y_edges = np.meshgrid(x_edges, y_edges)
pc = plt.pcolor(X_edges, Y_edges, H)
plt.colorbar(pc)
plt.show()
plt.close()

unique, counts = np.unique(example_source_indexes, return_counts=True)
sc = plt.scatter(pop_locs_X[source][unique], pop_locs_Y[source][unique], c=counts, linewidths=0)
plt.colorbar(sc)
plt.show()
plt.close()
"""