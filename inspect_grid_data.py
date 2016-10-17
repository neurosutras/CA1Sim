import numpy as np
import matplotlib.pyplot as plt
import h5py

data_dir = 'data/'

"""
f:
    groups:
        key = str(grid_id)
        attrs:
            'xoff': r_off * cos(phi_off)
            'yoff': r_off * sin(phi_off)
            'lambda': lambda of module
            'theta': theta of module
"""

# rec_filename = '101416 generated grid data'
# rec_filename = '100716 grid data'
rec_filename = '101716 modified grid data'
grid_data = {}

with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
    for key in f:
        grid_data[int(key)] = {}
        for param_name, value in f[key].attrs.iteritems():
            grid_data[int(key)][param_name] = value

hist, edges = {}, {}
for param in grid_data[0].keys():
    hist[param], edges[param] = np.histogram([this_grid_data[param] for this_grid_data in grid_data.itervalues()], 200)

for param in hist:
    plt.plot(edges[param][1:], hist[param])
    plt.title(param)
    plt.show()
    plt.close()

plt.scatter(*zip(*[(this_grid_data['xoff'], this_grid_data['yoff']) for this_grid_data in grid_data.itervalues()]))
plt.show()
plt.close()