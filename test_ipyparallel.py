__author__ = 'milsteina'
from function_lib import *
from ipyparallel import Client
from ipyparallel import interactive
import click

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

    :param cluster_id: str
    """
    c = Client(cluster_id=cluster_id)
    dv = c[:]
    print 'Controller sees %i engines' % len(dv.targets)
    dv.execute('from test_ipyparallel import *', block=True)
    result = dv.map_sync(get_engine_ids, dv.targets)
    print result


@interactive
def get_engine_ids(target_id):
    """

    :return:
    """
    pid = os.getpid()
    return {'target_id': target_id, 'pid': pid}


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1,sys.argv)+1):])
