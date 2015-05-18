__author__ = 'milsteina'
import os
from mpi4py import MPI
"""
This file is imported as a module on all ipengines. The method 'testMPI' returns a dictionary, which will be combined
and recorded to an output .hdf5 file by the controller.
"""


def testMPI(index):
    return dict(index=index, os_pid=os.getpid(), MPI_rank=MPI.COMM_WORLD.Get_rank())