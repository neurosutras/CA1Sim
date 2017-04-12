from plot_results import *
import sys
from scipy import interpolate
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

neurotrees_dir = morph_dir
# neurotrees_dir = os.environ['PI_SCRATCH']+'/DGC_forest/hdf5/'
# neurotrees_dir = os.environ['PI_HOME']+'/'
# forest_file = '122016_DGC_forest_with_syn_locs.h5'
forest_file = 'DGC_forest_connectivity.h5'
degrees_file = 'DGC_forest_connectivity_degrees_new.h5'

coords_dir = morph_dir
# coords_dir = os.environ['PI_SCRATCH']+'/DG/'
# coords_dir = os.environ['PI_HOME']+'/'
coords_file = 'DG_Soma_Coordinates.h5'

# forest_file = '100716_dentate_MPPtoDGC.h5'
# degrees_file = 'DGC_forest_connectivity_degrees_orig.h5'


# target_GID = synapse_dict.keys()
# target_GID.sort()

# print 'MPI rank %i received %i GCs: [%i:%i]' % (rank, len(target_GID), target_GID[0], target_GID[-1])
# sys.stdout.flush()

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

get_array_index = np.vectorize(lambda val_array, this_val: np.where(val_array >= this_val)[0][0], excluded=[0])


def plot_in_degree_single_target(target_gid, target, source, forest_file, soma_coords, u, v, distance_U, distance_V):
    """

    :param target_gid: int
    :param target: str
    :param source: str
    :param forest_file: str
    :param soma_coords: dict
    :param u: array
    :param v: array
    :param distance_U: array
    :param distance_V: array
    """
    raw_source_gids = soma_coords[source]['GID'][:]
    sorted_source_indexes = range(len(raw_source_gids))
    sorted_source_indexes.sort(key=raw_source_gids.__getitem__)
    with h5py.File(neurotrees_dir+forest_file, 'r') as f:
        target_connection_index = np.where(f['Populations'][target]['Connectivity']['source_gid']['gid'][:] ==
                                           target_gid)[0][0]
        start_index = f['Populations'][target]['Connectivity']['source_gid']['ptr'][target_connection_index]
        end_index = f['Populations'][target]['Connectivity']['source_gid']['ptr'][target_connection_index + 1]
        source_indexes_gid = get_array_index(raw_source_gids[sorted_source_indexes],
                                             f['Populations'][target]['Connectivity']['source_gid']['value'][
                                             start_index:end_index])
    source_indexes_u = get_array_index(u, soma_coords[source]['u'][sorted_source_indexes][source_indexes_gid])
    source_indexes_v = get_array_index(v, soma_coords[source]['v'][sorted_source_indexes][source_indexes_gid])

    fig1 = plt.figure(1)
    H, u_edges, v_edges = np.histogram2d(distance_U[source_indexes_u, source_indexes_v],
                                         distance_V[source_indexes_u, source_indexes_v], 100)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H == 0, H)
    ax = plt.gca()
    pcm = ax.pcolormesh(u_edges, v_edges, Hmasked)
    ax.set_xlabel('Arc distance (Septal - Temporal) (um)')
    ax.set_ylabel('Arc distance (Suprapyramidal - Infrapyramidal)  (um)')
    ax.set_title(target + ' <-- ' + source + ' (in degree)')
    ax.set_aspect('equal', 'box')
    clean_axes(ax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar = plt.colorbar(pcm, cax=cax)
    cbar.ax.set_ylabel('Counts (100 bins)')

    fig2 = plt.figure(2)
    unique_source_indexes_gid, counts = np.unique(source_indexes_gid, return_counts=True)
    unique_source_indexes_u = \
        get_array_index(u, soma_coords[source]['u'][sorted_source_indexes][unique_source_indexes_gid])
    unique_source_indexes_v = \
        get_array_index(v, soma_coords[source]['v'][sorted_source_indexes][unique_source_indexes_gid])
    ax = plt.gca()
    ax.scatter(distance_U[unique_source_indexes_u, unique_source_indexes_v],
                     distance_V[unique_source_indexes_u, unique_source_indexes_v], c=counts, linewidths=0)
    ax.set_xlabel('Arc distance (Septal - Temporal) (um)')
    ax.set_ylabel('Arc distance (Suprapyramidal - Infrapyramidal)  (um)')
    ax.set_title(target + ' <-- ' + source + ' (multiple innervation)')
    ax.set_aspect('equal', 'box')
    clean_axes(ax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar = plt.colorbar(pcm, cax=cax)
    cbar.ax.set_ylabel('Counts')

    plt.show()
    plt.close()


def plot_in_degree_single_target_orig(target_gid, target, source, projection, forest_file, width_U):
    """
    TODO: Requires parallel read.
    :param target_gid: int
    :param target: str
    :param source: str
    :param projection: str
    :param forest_file: str
    :param width_U: array

    with h5py.File(neurotrees_dir + forest_file, 'r') as f:
        start_index = 0
        while start_index < len(f['Projections'][projection]['Destination']):
            if f['Projections'][projection]['Destination']['']
        raw_source_gids = soma_coords[source]['GID'][:]
        sorted_source_indexes = range(len(raw_source_gids))
        sorted_source_indexes.sort(key=raw_source_gids.__getitem__)

        source_distances_U = np.arange(0., width_U, width_U / len(raw_source_gids))


            target_connection_index = np.where(f['Populations'][target]['Connectivity']['source_gid']['gid'][:] ==
                                               target_gid)[0][0]
            start_index = f['Populations'][target]['Connectivity']['source_gid']['ptr'][target_connection_index]
            end_index = f['Populations'][target]['Connectivity']['source_gid']['ptr'][target_connection_index + 1]
            source_indexes_gid = get_array_index(raw_source_gids[sorted_source_indexes],
                                                 f['Populations'][target]['Connectivity']['source_gid']['value'][
                                                 start_index:end_index])
        source_indexes_u = get_array_index(u, soma_coords[source]['u'][sorted_source_indexes][source_indexes_gid])
        source_indexes_v = get_array_index(v, soma_coords[source]['v'][sorted_source_indexes][source_indexes_gid])

        fig1 = plt.figure(1)
        H, u_edges, v_edges = np.histogram2d(distance_U[source_indexes_u, source_indexes_v],
                                             distance_V[source_indexes_u, source_indexes_v], 100)
        H = np.rot90(H)
        H = np.flipud(H)
        Hmasked = np.ma.masked_where(H == 0, H)
        ax = plt.gca()
        pcm = ax.pcolormesh(u_edges, v_edges, Hmasked)
        ax.set_xlabel('Arc distance (Septal - Temporal) (um)')
        ax.set_ylabel('Arc distance (Suprapyramidal - Infrapyramidal)  (um)')
        ax.set_title(target + ' <-- ' + source + ' (in degree)')
        ax.set_aspect('equal', 'box')
        clean_axes(ax)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        cbar = plt.colorbar(pcm, cax=cax)
        cbar.ax.set_ylabel('Counts (100 bins)')

        fig2 = plt.figure(2)
        unique_source_indexes_gid, counts = np.unique(source_indexes_gid, return_counts=True)
        unique_source_indexes_u = \
            get_array_index(u, soma_coords[source]['u'][sorted_source_indexes][unique_source_indexes_gid])
        unique_source_indexes_v = \
            get_array_index(v, soma_coords[source]['v'][sorted_source_indexes][unique_source_indexes_gid])
        ax = plt.gca()
        ax.scatter(distance_U[unique_source_indexes_u, unique_source_indexes_v],
                         distance_V[unique_source_indexes_u, unique_source_indexes_v], c=counts, linewidths=0)
        ax.set_xlabel('Arc distance (Septal - Temporal) (um)')
        ax.set_ylabel('Arc distance (Suprapyramidal - Infrapyramidal)  (um)')
        ax.set_title(target + ' <-- ' + source + ' (multiple innervation)')
        ax.set_aspect('equal', 'box')
        clean_axes(ax)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        cbar = plt.colorbar(pcm, cax=cax)
        cbar.ax.set_ylabel('Counts')

        plt.show()
        plt.close()
    """

def plot_out_degree_single_source(source_gid, target, source, forest_file, soma_coords, u, v, distance_U, distance_V):
    """
    TODO: This requires parallel read. Should just write a script that inverts the connectivity.
    :param source_gid: int
    :param target: str
    :param source: str
    :param forest_file: str
    :param soma_coords: dict
    :param u: array
    :param v: array
    :param distance_U: array
    :param distance_V: array
    """
    raw_target_gids = soma_coords[target]['GID'][:]
    sorted_target_indexes = range(len(raw_target_gids))
    sorted_target_indexes.sort(key=raw_target_gids.__getitem__)
    with h5py.File(neurotrees_dir+forest_file, 'r') as f:
        target_connection_index = np.where(f['Populations'][target]['Connectivity']['source_gid']['gid'][:] ==
                                           target_gid)[0][0]
        start_index = f['Populations'][target]['Connectivity']['source_gid']['ptr'][target_connection_index]
        end_index = f['Populations'][target]['Connectivity']['source_gid']['ptr'][target_connection_index + 1]
        source_indexes_gid = get_array_index(raw_source_gids[sorted_source_indexes],
                                             f['Populations'][target]['Connectivity']['source_gid']['value'][
                                             start_index:end_index])
    source_indexes_u = get_array_index(u, soma_coords[source]['u'][sorted_source_indexes][source_indexes_gid])
    source_indexes_v = get_array_index(v, soma_coords[source]['v'][sorted_source_indexes][source_indexes_gid])

    fig1 = plt.figure(1)
    H, u_edges, v_edges = np.histogram2d(distance_U[source_indexes_u, source_indexes_v],
                                         distance_V[source_indexes_u, source_indexes_v], 100)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H == 0, H)
    ax = plt.gca()
    pcm = ax.pcolormesh(u_edges, v_edges, Hmasked)
    ax.set_xlabel('Arc distance (Septal - Temporal) (um)')
    ax.set_ylabel('Arc distance (Suprapyramidal - Infrapyramidal)  (um)')
    ax.set_title(target + ' <-- ' + source + ' (in degree)')
    ax.set_aspect('equal', 'box')
    clean_axes(ax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar = plt.colorbar(pcm, cax=cax)
    cbar.ax.set_ylabel('Counts (100 bins)')

    fig2 = plt.figure(2)
    unique_source_indexes_gid, counts = np.unique(source_indexes_gid, return_counts=True)
    unique_source_indexes_u = \
        get_array_index(u, soma_coords[source]['u'][sorted_source_indexes][unique_source_indexes_gid])
    unique_source_indexes_v = \
        get_array_index(v, soma_coords[source]['v'][sorted_source_indexes][unique_source_indexes_gid])
    ax = plt.gca()
    ax.scatter(distance_U[unique_source_indexes_u, unique_source_indexes_v],
                     distance_V[unique_source_indexes_u, unique_source_indexes_v], c=counts, linewidths=0)
    ax.set_xlabel('Arc distance (Septal - Temporal) (um)')
    ax.set_ylabel('Arc distance (Suprapyramidal - Infrapyramidal)  (um)')
    ax.set_title(target + ' <-- ' + source + ' (multiple innervation)')
    ax.set_aspect('equal', 'box')
    clean_axes(ax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar = plt.colorbar(pcm, cax=cax)
    cbar.ax.set_ylabel('Counts')

    plt.show()
    plt.close()


def plot_in_degree(target, source, projection, degrees_file, soma_coords, u, v, distance_U, distance_V, bin_size=50.):
    """

    :param target: str
    :param source: str
    :param projection: str
    :param degrees_file: :class:'h5py.Group'
    :param soma_coords: dict
    :param u: array
    :param v: array
    :param distance_U: array
    :param distance_V: array
    :param bin_size: float (um)
    """
    with h5py.File(neurotrees_dir+degrees_file, 'r') as f:
        in_degree_group = f['Projections'][projection]['In degree']
        raw_target_gids = soma_coords[target]['GID'][:]
        sorted_target_indexes = range(len(raw_target_gids))
        sorted_target_indexes.sort(key=raw_target_gids.__getitem__)

        target_indexes_gid = get_array_index(raw_target_gids[sorted_target_indexes], in_degree_group['gid'][:])
        target_indexes_u = get_array_index(u, soma_coords[target]['u'][sorted_target_indexes][target_indexes_gid])
        target_indexes_v = get_array_index(v, soma_coords[target]['v'][sorted_target_indexes][target_indexes_gid])

        this_step_size = int(bin_size / spatial_resolution)
        this_in_degree = interpolate.griddata((distance_U[target_indexes_u, target_indexes_v],
                                                distance_V[target_indexes_u, target_indexes_v]),
                                               in_degree_group['count'][:],
                                               (distance_U[::this_step_size, ::this_step_size],
                                                distance_V[::this_step_size, ::this_step_size]), method='nearest')
        fig1 = plt.figure(1, figsize=plt.figaspect(1.) * 2.)
        ax = plt.gca()
        pcm = plt.pcolormesh(distance_U[::this_step_size, ::this_step_size],
                             distance_V[::this_step_size, ::this_step_size], this_in_degree)
        ax.set_xlabel('Arc distance (septal - temporal) (um)')
        ax.set_ylabel('Arc distance (supra - infrapyramidal)  (um)')
        ax.set_title(target+' <-- '+source+' (in degree)')
        ax.set_aspect('equal', 'box')
        clean_axes(ax)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        cbar = plt.colorbar(pcm, cax=cax)
        cbar.ax.set_ylabel('Counts')

        plt.show()
        plt.close()


def plot_out_degree(target, source, projection, degrees_file, soma_coords, u, v, distance_U, distance_V, bin_size=50.):
    """

    :param target: str
    :param source: str
    :param projection: str
    :param degrees_file: :class:'h5py.Group'
    :param soma_coords: dict
    :param u: array
    :param v: array
    :param distance_U: array
    :param distance_V: array
    :param bin_size: float (um)
    """
    with h5py.File(neurotrees_dir+degrees_file, 'r') as f:
        out_degree_group = f['Projections'][projection]['Out degree']
        raw_source_gids = soma_coords[source]['GID'][:]
        sorted_source_indexes = range(len(raw_source_gids))
        sorted_source_indexes.sort(key=raw_source_gids.__getitem__)

        source_indexes_gid = get_array_index(raw_source_gids[sorted_source_indexes], out_degree_group['gid'][:])
        source_indexes_u = get_array_index(u, soma_coords[source]['u'][sorted_source_indexes][source_indexes_gid])
        source_indexes_v = get_array_index(v, soma_coords[source]['v'][sorted_source_indexes][source_indexes_gid])

        this_step_size = int(bin_size / spatial_resolution)
        this_out_degree = interpolate.griddata((distance_U[source_indexes_u, source_indexes_v],
                                               distance_V[source_indexes_u, source_indexes_v]),
                                               out_degree_group['count'][:],
                                               (distance_U[::this_step_size, ::this_step_size],
                                                distance_V[::this_step_size, ::this_step_size]), method='nearest')
        fig1 = plt.figure(1, figsize=plt.figaspect(1.) * 2.)
        ax = plt.gca()
        pcm = plt.pcolormesh(distance_U[::this_step_size, ::this_step_size],
                             distance_V[::this_step_size, ::this_step_size], this_out_degree)
        ax.set_xlabel('Arc distance (septal - temporal) (um)')
        ax.set_ylabel('Arc distance (supra - infrapyramidal)  (um)')
        ax.set_title(source+' --> '+target+' (out degree)')
        ax.set_aspect('equal', 'box')
        clean_axes(ax)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        cbar = plt.colorbar(pcm, cax=cax)
        cbar.ax.set_ylabel('Counts')

        plt.show()
        plt.close()


def plot_out_degree_orig(target, source, projection, degrees_file, width_U):
    """

    :param target: str
    :param source: str
    :param projection: str
    :param degrees_file: :class:'h5py.Group'
    :param width_U: array
    """
    with h5py.File(neurotrees_dir+degrees_file, 'r') as f:
        out_degree_group = f['Projections'][projection]['Out degree']
        raw_source_gids = out_degree_group['gid'][:]
        sorted_source_indexes = range(len(raw_source_gids))
        sorted_source_indexes.sort(key=raw_source_gids.__getitem__)
        source_distances_U = np.arange(0., width_U, width_U / len(raw_source_gids))

        fig1 = plt.figure(1, figsize=plt.figaspect(1.) * 2.)
        ax = plt.gca()
        ax.scatter(source_distances_U[sorted_source_indexes], out_degree_group['count'][:], linewidths=0)
        ax.set_xlabel('Arc distance (Septal - Temporal) (um)')
        ax.set_ylabel('Synapse count')
        ax.set_title(source+' --> '+target+' (out degree)')
        clean_axes(ax)

        plt.show()
        plt.close()


def plot_population_density(population, soma_coords, u, v, U, V, distance_U, distance_V, max_u, max_v, bin_size=50.):
    """

    :param population: str
    :param soma_coords: dict of array
    :param u: array
    :param v: array
    :param U: array
    :param V: array
    :param distance_U: array
    :param distance_V: array
    :param max_u: float: u_distance
    :param max_v: float: v_distance
    :param bin_size: float
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
    c
    scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]] * 3)
    ax.set_xlabel('X (um)')
    ax.set_ylabel('Y (um)')
    ax.set_zlabel('Z (um)')

    step_sizes = [int(max_u / bin_size), int(max_v / bin_size)]
    fig2 = plt.figure(2, figsize=plt.figaspect(1.)*2.)
    population_indexes_u = get_array_index(u, soma_coords[population]['u'])
    population_indexes_v = get_array_index(v, soma_coords[population]['v'])
    H, u_edges, v_edges = np.histogram2d(distance_U[population_indexes_u, population_indexes_v],
                                         distance_V[population_indexes_u, population_indexes_v], step_sizes)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H == 0, H)
    ax = plt.gca()
    pcm = ax.pcolormesh(u_edges, v_edges, Hmasked)
    ax.set_xlabel('Arc distance (septal - temporal) (um)')
    ax.set_ylabel('Arc distance (supra - infrapyramidal)  (um)')
    ax.set_title(population)
    ax.set_aspect('equal', 'box')
    clean_axes(ax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar = plt.colorbar(pcm, cax=cax)
    cbar.ax.set_ylabel('Counts')

    plt.show()
    plt.close()


# plot_population_density('GC', soma_coords, u, v, U, V, distance_U, distance_V, max_u, max_v)

plot_population_density('MEC', soma_coords, u, v, U, V, distance_U, distance_V, max_u, max_v)

# plot_in_degree('GC', 'MEC', 'MPPtoGC', degrees_file, soma_coords, u, v, distance_U, distance_V)


# plot_out_degree('GC', 'MEC', 'MPPtoGC', degrees_file, soma_coords, u, v, distance_U, distance_V)
# plot_out_degree_orig('GC', 'MEC', 'MPPtoGC', degrees_file, width_U)
"""
target_gid_index = np.where((soma_coords['GC']['u'] > u[0] + (u[-1]-u[0])/2.) &
                            (soma_coords['GC']['v'] > v[0] + (v[-1]-v[0])/2.))[0][0]
target_gid = soma_coords['GC']['GID'][target_gid_index]
# plot_in_degree_single_target(target_gid, 'GC', 'MEC', forest_file, soma_coords, u, v, distance_U, distance_V)
plot_in_degree_single_target_orig(target_gid, 'GC', 'MEC', 'MPPtoGC', forest_file, width_U)


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
