import numpy as np
import matplotlib.pyplot as plt
import h5py
import random
import time


data_dir = 'data/'

"""
f:
    groups:
        key = str(grid_id)
        attrs:
            'x_off': r_off * cos(phi_off)
            'y_off': r_off * sin(phi_off)
            'module': int(index of module)
            'lambda': arange[40., 400., 40.]
            'theta': [0., np/pi/3.]

"""

rec_filename = '101716 generated grid data'
local_random = random.Random()
local_random.seed(1e12)
num_modules = 10
num_grid_cells = 38000
lambda_array = np.arange(40., (num_modules+1)*40., 40.)
# every 60 degrees repeats in a hexagonal array
theta_array = [local_random.uniform(0., np.pi/3.) for i in range(num_modules)]

start_time = time.time()
with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:

    i_time = 0
    for grid_id in range(num_grid_cells):
        f.create_group(str(grid_id))
        this_module = local_random.choice(range(num_modules))
        this_lambda = lambda_array[this_module]
        f[str(grid_id)].attrs['lambda'] = this_lambda
        f[str(grid_id)].attrs['theta'] = theta_array[this_module]
        r_off = this_lambda * np.sqrt(local_random.random())
        phi_off = local_random.uniform(-np.pi, np.pi)
        x_off = r_off * np.cos(phi_off)
        y_off = r_off * np.sin(phi_off)
        f[str(grid_id)].attrs['xoff'] = x_off
        f[str(grid_id)].attrs['yoff'] = y_off
        i_time += 1
        if (i_time % 1000 == 0):
            print 'parameters generated for', i_time, 'grid cells'
print 'generating grid data for', i_time, 'grid cells took', time.time() - start_time, 's'

