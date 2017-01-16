from function_lib import *
from mpi4py import MPI
from neurotrees.io import append_tree_attributes
from neurotrees.io import read_tree_attributes
# import mkl
import sys
from plot_results import *
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
                                      'l': f['Coordinates'][import_pop_key]['L Coordinate'][:],
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


def plot_source_soma_locs(target_gid, target, source, connection_dict, soma_coords, u, v, U, V, distance_U,
                          distance_V, X, Y, Z):
    """

    :param target_gid: int
    :param target: str
    :param source: str
    :param connection_dict: dict
    :param soma_coords: dict
    :param u: array
    :param v: array
    :param U: array
    :param V: array
    :param distance_U: array
    :param distance_V: array
    :param X: array
    :param Y: array
    :param Z: array
    """
    target_gid_index = np.where(soma_coords[target]['GID'] == target_gid)[0][0]
    target_index_u = np.where(u >= soma_coords[target]['u'][target_gid_index])[0][0]
    target_index_v = np.where(v >= soma_coords[target]['v'][target_gid_index])[0][0]
    source_gids = connection_dict[target_gid]['source_gid']
    source_indexes_gid = get_array_index(soma_coords[source]['GID'], source_gids)
    source_indexes_u = get_array_index(u, soma_coords[source]['u'][source_indexes_gid])
    source_indexes_v = get_array_index(v, soma_coords[source]['v'][source_indexes_gid])

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111, projection='3d')
    ax.scatter(X[::100, ::50], Y[::100, ::50], Z[::100, ::50], c='grey', alpha=0.1, zorder=0)
    ax.scatter(X[source_indexes_u, source_indexes_v], Y[source_indexes_u, source_indexes_v],
               Z[source_indexes_u, source_indexes_v], c='b', zorder=1)
    ax.scatter(X[target_index_u, target_index_v], Y[target_index_u, target_index_v],
               Z[target_index_u, target_index_v], c='r', s=10, zorder=2)

    fig2 = plt.figure(2)
    H, u_edges, v_edges = np.histogram2d(soma_coords[source]['u'][source_indexes_gid],
                                         soma_coords[source]['v'][source_indexes_gid], 100)
    # H, u_edges, v_edges = np.histogram2d(u[source_indexes_u], v[source_indexes_v], 100)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H == 0, H)
    plt.pcolormesh(u_edges, v_edges, Hmasked)
    plt.xlabel('U')
    plt.ylabel('V')
    plt.title(source+':'+target)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')

    fig3 = plt.figure(3)
    H, u_edges, v_edges = np.histogram2d(distance_U[source_indexes_u, source_indexes_v],
                                         distance_V[source_indexes_u, source_indexes_v], 100)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H == 0, H)
    plt.pcolormesh(u_edges, v_edges, Hmasked)
    plt.xlabel('Arc_distance U (um)')
    plt.ylabel('Arc_distance V (um)')
    plt.title(source + ':' + target)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')

    fig4 = plt.figure(4)
    unique_source_indexes_gid, counts = np.unique(source_indexes_gid, return_counts=True)
    unique_source_indexes_u = get_array_index(u, soma_coords[source]['u'][unique_source_indexes_gid])
    unique_source_indexes_v = get_array_index(v, soma_coords[source]['v'][unique_source_indexes_gid])
    plt.scatter(distance_U[unique_source_indexes_u, unique_source_indexes_v],
                     distance_V[unique_source_indexes_u, unique_source_indexes_v], c=counts, linewidths=0)
    plt.xlabel('Arc_distance U (um)')
    plt.ylabel('Arc_distance V (um)')
    plt.title(source + ':' + target)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    plt.show()
    plt.close()


def plot_population_density(population, soma_coords, u, v, U, V, distance_U, distance_V):
    """

    :param population: str
    :param soma_coords: dict of array
    :param u: array
    :param v: array
    :param U: array
    :param V: array
    :param distance_U: array
    :param distance_V: array
    :return:
    """
    """
    fig1 = plt.figure(1)
    H, u_edges, v_edges = np.histogram2d(soma_coords[population]['u'], soma_coords[population]['v'], 100)
    # H, u_edges, v_edges = np.histogram2d(u[population_indexes_u], v[population_indexes_v], 100)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H == 0, H)
    plt.pcolormesh(u_edges, v_edges, Hmasked)
    plt.xlabel('U')
    plt.ylabel('V')
    plt.title(population)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    """
    x = np.vectorize(lambda this_u, this_v, this_l: -500. * np.cos(this_u) * (5.3 - np.sin(this_u) +
                                                                              (1. + 0.138 * this_l) * np.cos(this_v)))
    y = np.vectorize(lambda this_u, this_v, this_l: 750. * np.sin(this_u) * (5.5 - 2. * np.sin(this_u) +
                                                                             (0.9 + 0.114 * this_l) * np.cos(this_v)))
    z = np.vectorize(lambda this_u, this_v, this_l: 2500. * np.sin(this_u) + (663. + 114. * this_l) *
                                                                             np.sin(this_v - 0.13 * (np.pi - this_u)))

    soma_coords[population]['x'] = x(soma_coords[population]['u'], soma_coords[population]['v'],
                                     soma_coords[population]['l'])
    soma_coords[population]['y'] = y(soma_coords[population]['u'], soma_coords[population]['v'],
                                     soma_coords[population]['l'])
    soma_coords[population]['z'] = z(soma_coords[population]['u'], soma_coords[population]['v'],
                                     soma_coords[population]['l'])

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111, projection='3d')
    piece_size = int(len(soma_coords[population]['x'])/10000)
    ax.scatter(soma_coords[population]['x'][::piece_size], soma_coords[population]['y'][::piece_size],
               soma_coords[population]['z'][::piece_size], alpha=0.1, linewidth=0)
    scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]] * 3)
    ax.set_xlabel('X (um)')
    ax.set_ylabel('Y (um)')
    ax.set_zlabel('Z (um)')

    fig2 = plt.figure(2, figsize=plt.figaspect(1.)*2.)
    population_indexes_u = get_array_index(u, soma_coords[population]['u'])
    population_indexes_v = get_array_index(v, soma_coords[population]['v'])
    H, u_edges, v_edges = np.histogram2d(distance_U[population_indexes_u, population_indexes_v],
                                         distance_V[population_indexes_u, population_indexes_v], 100)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H == 0, H)
    ax = plt.gca()
    pcm = ax.pcolormesh(u_edges, v_edges, Hmasked)
    ax.set_xlabel('Arc distance (Septal - Temporal) (um)')
    ax.set_ylabel('Arc distance (Suprapyramidal - Infrapyramidal)  (um)')
    ax.set_title(population)
    ax.set_aspect('equal', 'box')
    clean_axes(ax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar = plt.colorbar(pcm, cax=cax)
    cbar.ax.set_ylabel('Counts')

    plt.show()
    plt.close()


plot_population_density('GC', soma_coords, u, v, U, V, distance_U, distance_V)
plot_population_density('MEC', soma_coords, u, v, U, V, distance_U, distance_V)

"""
connection_dict = read_from_pkl(neurotrees_dir+'010117_DG_GC_MEC_connectivity_test.pkl')
target_gid = 500
plot_source_soma_locs(target_gid, 'GC', 'MEC', connection_dict, soma_coords, u, v, U, V, distance_U, distance_V,
                     X, Y, Z)


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


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(X[::250], Y[::250], Z[::250], c=distance_U[::250], linewidth=0)
plt.colorbar(sc)
plt.show()
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(X[::250], Y[::250], Z[::250], c=distance_V[::250], linewidth=0)
plt.colorbar(sc)
plt.show()
plt.close()

"""
