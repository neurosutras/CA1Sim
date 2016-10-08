__author__ = 'Aaron D. Milstein'
import numpy as np
import matplotlib.pyplot as plt
import h5py
from ipyparallel import Client
from IPython.display import clear_output
import sys
import analyze_dentate_connectivity_engine
import time


data_dir = 'data/'

rec_filename = '100716 grid cell waveforms'
grid_data_filename = '100716 grid data'

start_time = time.time()

if len(sys.argv) > 1:
    cluster_id = sys.argv[1]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

grid_data = {}

with h5py.File(data_dir+grid_data_filename+'.hdf5', 'r') as f:
    for key in f:
        grid_data[int(key)] = {}
        for param_name, value in f[key].attrs.iteritems():
            grid_data[int(key)][param_name] = value

dv = c[:]
dv.clear()
dv.block = True
dv.execute('from analyze_dentate_connectivity_engine import *')
v = c.load_balanced_view()
result = v.map_async(analyze_dentate_connectivity_engine.compute_single_rate_map, grid_data.iteritems())
"""
result = v.map_async(analyze_dentate_connectivity_engine.compute_single_rate_map, [(i, grid_data[i]) for
                                                                                   i in range(16)])
"""
while not result.ready():
    time.sleep(30)
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
    clear_output()

with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
    for id, waveforms in result.get():
        f.create_group(str(id))
        for title, waveform in waveforms.iteritems():
            f[str(id)].create_dataset(title, compression='gzip', compression_opts=9, data=waveform)

print len(c), 'processes computed the set of grid cell rate maps in parallel in', time.time() - start_time, 's'