__author__ = 'milsteina'
from testMPI_engine import *
from ipyparallel import Client
import h5py
"""
This file is executed on the controller. The method 'testMPI' returns a dictionary from each engine, which is combined
and recorded to an output .hdf5 file.
"""
rec_filename = '052015 testMPI output'
data_dir = 'data/'

#c = Client(profile='mpi')
c = Client()
dv = c[:]
dv.block = True
dv.execute('from testMPI_engine import *')
result = dv.map(testMPI, range(len(c)))
with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
    for i in range(len(result)):
        g = f.create_group(str(i))
        for key in result[i]:
            g.attrs[key] = result[i][key]