__author__ = 'Grace Ng'
import click
from function_lib import *
from ipyparallel import Client
import importlib


script_filename = 'to_console.py'

module_names = ['to_console2', 'to_console3']
module_config_func_names = ['config2', 'config3']
module_compute_func_names = ['compute2', 'compute3']


@click.command()
@click.option('--cluster-id', required=True)
def main(cluster_id):
    global c
    c = Client(cluster_id=cluster_id)
    module_config_func_refs, module_compute_func_refs = \
        config_controller(module_names, module_config_func_names, module_compute_func_names)
    c[:].execute('from to_console import *', block=True)
    result = c[:].apply_async(init_engines, module_names, module_config_func_refs)
    while not result.ready():
        time.sleep(1)
        sys.stdout.flush()
    for line in result.stdout:
        print line
    client_ranges = [xrange(i, i+2, 1) for i in xrange(0, len(c), 2)]
    for compute_func in module_compute_func_refs:
        result = map(compute_func, [c] * len(client_ranges), client_ranges)
        print result


def config_controller(module_names, module_config_func_names, module_compute_func_names):
    module_config_func_refs = []
    module_compute_func_refs = []
    for module_name in set(module_names):
        importlib.import_module(module_name)
    for module_name, config_func_name, compute_func_name in \
            zip(module_names, module_config_func_names, module_compute_func_names):
        this_config_func_ref = getattr(sys.modules[module_name], config_func_name)
        if not callable(this_config_func_ref):
            raise Exception('config_controller: %s.%s is not callable' % (module_name, config_func_name))
        else:
            module_config_func_refs.append(this_config_func_ref)
        this_compute_func_ref = getattr(sys.modules[module_name], compute_func_name)
        if not callable(this_compute_func_ref):
            raise Exception('compute_controller: %s.%s is not callable' % (module_name, compute_func_name))
        else:
            module_compute_func_refs.append(this_compute_func_ref)
    return module_config_func_refs, module_compute_func_refs


def init_engines(module_names, module_config_func_refs):
    for module_name in set(module_names):
        importlib.import_module(module_name)
    for module_name in set(module_names):
        this_init_func = getattr(sys.modules[module_name], 'init_engine')
        if not callable(this_init_func):
            raise Exception('init_engines: %s.%s is not callable' % (module_name, 'init_engine'))
        else:
            this_init_func(module_config_func_refs)


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1, sys.argv) + 1):])
