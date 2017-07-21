__author__ = 'Grace Ng'
from function_lib import *


module_name = 'to_console2'

context = Context()
context.update(locals())


def init_engine(module_config_func_refs):
    context.update(locals())
    print '%s.init_engine executed' % (context.module_name)


def config2(id):
    context.update(locals())
    print '%s.config2: id: %i stored' % (context.module_name, id)


def compute2(c, client_ranges):
    result = c[client_ranges].map_async(set_and_report, client_ranges)
    while not result.ready():
        time.sleep(1)
        sys.stdout.flush()
    for line in result.stdout:
        print line
    return result.get()


def set_and_report(id):
    for config_func in context.module_config_func_refs:
        config_func(id)
    print '%s.set_and_report: id: %i; context.id: %i' % (context.module_name, id, context.id)
    return id, context.id, context.module_name
