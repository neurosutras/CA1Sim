__author__ = 'milsteina'
from function_lib import *
from ipyparallel import Client
import click
from mpi4py import MPI

"""
Tests that complete range of engines, potentially across multiple nodes on a cluster, are accessible.
Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

script_filename = 'test_ipyparallel.py'


@click.command()
@click.option("--cluster-id", type=str, required=True)
def main(cluster_id):
    """

    :param start: int
    :param cluster_id: str
    """
    c = Client(cluster_id=cluster_id)

    print c.ids

    dv = c[:]
    result = dv.map_sync(get_engine_ids, range(len(c)))
    print result


def get_engine_ids(index):
    """

    :return:
    """
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    pid = os.getpid()

    return {'index': index, 'rank': rank, 'size': size, 'pid': pid}


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1,sys.argv)+1):])