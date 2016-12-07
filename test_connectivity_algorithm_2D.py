from plot_results import *

"""
Distribute GC-BC and BC-GC synapses within a 60 um 1D section of the 6000 um S-T axis of the dentate.
Total probability across axonal extent should be normalized to hit target divergence counts.

Algorithm:
1. For each presynaptic population,
    i. Convert the (u, v) coordinates of their soma locations into (st, ml) where the distances are equivalent in arc
        length within the surface of the layer.
    ii. Specify the septal-temporal and medial-lateral axon widths.
2. For each cell, for each synapse type, choose the all the presynaptic sources with one function call rather than
    iterating. Populate a list of candidate sources based on relative (st, ml).
3. Luckily the shortest ml axis is > 3 mm long, so will fit the axonal extent of even the longest cell type. The
    normalized probability, then, shouldn't vary with septal-temporal position.

"""
rec_filename = '112016 GC MEC connectivity'

dx = 1.  # um
dy = 1.  # um
x_width = 6000.
y_width = 3000.  # this will be variable in the real surface
x = np.arange(0., x_width, dx)
y = np.arange(0., y_width, dy)
X, Y = np.meshgrid(x, y)

pop_size = {'GC': 1000000, 'MEC': 38000}  # of cells in a 60 um focal region

pop_density = {pop: pop_size[pop] / x_width / y_width for pop in pop_size}  # cell density per um^2
# Sloviter, Lomo, Frontiers, 2012. says GCs might only project 200 um S-T..., MEC less than 1500.
axon_width = {'GC': (900., 900.), 'MEC': (1500., 1500.)}  # full width in um (x, y)
# can't find good estimate of MEC divergence. circular reasoning: convergence onto GC * # of GCs / # of MEC cells
divergence = {'MEC': {'GC': 40000}}  # {source: {target: # of target cells contacted by 1 source cell}}
convergence = {'GC': {'MEC': 1500}}  # {target: {source: # of source cells that contact 1 target cell}}

pop_locs_X, pop_locs_Y = {}, {}
for pop in pop_size:
    unit_area = x_width * y_width / pop_size[pop]
    unit_radius = np.sqrt(0.9063 * unit_area / np.pi)
    x1 = np.arange(unit_radius, x_width, 2.*unit_radius)
    y1 = np.arange(unit_radius, y_width, 2.*np.sqrt(3.)*unit_radius)
    x2 = np.arange(2.*unit_radius, x_width, 2.*unit_radius)
    y2 = np.arange((1.+np.sqrt(3.))*unit_radius, y_width, 2.*np.sqrt(3.)*unit_radius)
    X1, Y1 = np.meshgrid(x1, y1)
    X2, Y2 = np.meshgrid(x2, y2)
    pop_locs_X[pop] = np.append(X1, X2)
    pop_locs_Y[pop] = np.append(Y1, Y2)


class AxonProb(object):
    """

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
            width = {'x': axon_width[source][0], 'y': axon_width[source][1]}
            sigma = {axis: width[axis] / 3. / np.sqrt(2.) for axis in width}
            for target in divergence[source]:
                surface_area = np.pi * width['x'] * width['y'] / 4.
                total_p = divergence[source][target] / surface_area / pop_density[target]
                this_volume = np.pi * sigma['x'] * sigma['y']
                self.p_dist[source][target] = np.vectorize(lambda target_x, target_y, source_x, source_y:
                                                   total_p *
                                                   np.exp(-(((source_x - target_x) / sigma['x'])**2. +
                                                            ((source_y - target_y) / sigma['y'])**2.)) /
                                                   this_volume)

    def choose_sources(self, target, source, target_x, target_y, num_syns, pop_locs_X, pop_locs_Y, local_random=None):
        """

        :param target: str
        :param source: str
        :param target_x: float
        :param target_y: float
        :param num_syns: int
        :param pop_locs_X: dict of float
        :param pop_locs_Y: dict of float
        :param local_random: :class:'np.random.RandomState'
        :return: list of int
        """
        if local_random is None:
            local_random = np.random.RandomState()
        source_indexes = np.where((np.abs(pop_locs_X[source] - target_x) <= axon_width[source][0]/2.) &
                                  (np.abs(pop_locs_Y[source] - target_y) <= axon_width[source][1]/2.))[0]
        source_p = self.p_dist[source][target](target_x, target_y, pop_locs_X[source][source_indexes],
                                               pop_locs_Y[source][source_indexes])
        source_p /= np.sum(source_p)
        return local_random.choice(source_indexes, num_syns, p=source_p)

p_connect = AxonProb(divergence, axon_width, pop_density)

local_np_random = np.random.RandomState()
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