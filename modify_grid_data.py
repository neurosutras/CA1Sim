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
            'xoff': r_off * cos(phi_off)
            'yoff': r_off * sin(phi_off)
            'lambda': arange[40., 400., 40.]
            'theta': [0., np/pi/3.]

"""

rec_filename = '101716 modified grid data'
local_random = random.Random()
local_random.seed(1e12)

start_time = time.time()
with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
    i_time = 0
    for grid_cell in f.values():
        this_lambda = grid_cell.attrs['lambda']
        r_off = this_lambda * np.sqrt(local_random.random())
        phi_off = local_random.uniform(-np.pi, np.pi)
        x_off = r_off * np.cos(phi_off)
        y_off = r_off * np.sin(phi_off)
        grid_cell.attrs['xoff'] = x_off
        grid_cell.attrs['yoff'] = y_off
        i_time += 1
        if (i_time % 1000 == 0):
            print 'parameters modified for', i_time, 'grid cells'
print 'modifying grid data for', i_time, 'grid cells took', time.time() - start_time, 's'
