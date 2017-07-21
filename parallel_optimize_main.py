__author__ = 'Grace Ng'
import click
from ipyparallel import interactive
from ipyparallel import Client
from moopgen import *
from plot_results import *
import importlib

"""
Dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass

script_filename = 'parallel_optimize_main.py'

main_ctxt = Context()

@click.command()
@click.option("--cluster-id", type=str, default=None)
@click.option("--profile", type=str, default='default')
@click.option("--spines", is_flag=True)
@click.option("--mech-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None)
@click.option("--neurotree-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None)
@click.option("--neurotree-index", type=int, default=0)
@click.option("--param-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False),
              default='data/optimize_spiking_defaults.yaml')
@click.option("--param-gen", type=str, default=None)
@click.option("--update-params", "-up", multiple=True, type=str, default=None)
@click.option("--update-modules", "-um", multiple=True, type=str, default=None)
@click.option("--get-features", "-gf", multiple=True, type=str, default=None)
@click.option("--features-modules", "-fm", multiple=True, type=str, default=None)
@click.option("--objectives-modules", "-om", multiple=True, type=str, default=None)
@click.option("--group-sizes", "-gs", multiple=True, type=int, default=None)
@click.option("--pop-size", type=int, default=100)
@click.option("--wrap-bounds", is_flag=True)
@click.option("--seed", type=int, default=None)
@click.option("--max-iter", type=int, default=None)
@click.option("--path-length", type=int, default=1)
@click.option("--initial-step-size", type=float, default=0.5)
@click.option("--adaptive-step-factor", type=float, default=0.9)
@click.option("--survival-rate", type=float, default=0.2)
@click.option("--sleep", is_flag=True)
@click.option("--analyze", is_flag=True)
@click.option("--hot-start", is_flag=True)
@click.option("--storage-file-path", type=str, default=None)
@click.option("--export", is_flag=True)
@click.option("--output-dir", type=str, default='data')
@click.option("--export-file-path", type=str, default=None)
@click.option("--disp", is_flag=True)
def main(cluster_id, profile, spines, mech_file_path, neurotree_file_path, neurotree_index, param_file_path,
         param_gen, update_params, update_modules, get_features, features_modules, objectives_modules, group_sizes,
         pop_size, wrap_bounds, seed, max_iter, path_length, initial_step_size, adaptive_step_factor, survival_rate,
         sleep, analyze, hot_start, storage_file_path, export, output_dir, export_file_path, disp):
    """

    :param cluster_id: str
    :param profile: str
    :param spines: bool
    :param mech_file_path: str (path)
    :param neurotree_file_path: str (path)
    :param neurotree_index: int
    :param param_file_path: str (path)
    :param param_gen: str (must refer to callable in globals())
    :param update_params: tuple of str (must refer to callables in globals())
    :param update_modules: tuple of str
    :param get_features: tuple of str (must refer to callables in globals())
    :param features_modules: tuple of str
    :param objectives_modules: tuple of str
    :param group_sizes: tuple of int
    :param pop_size: int
    :param wrap_bounds: bool
    :param seed: int
    :param max_iter: int
    :param path_length: int
    :param initial_step_size: float in [0., 1.]
    :param adaptive_step_factor: float in [0., 1.]
    :param survival_rate: float
    :param sleep: bool
    :param analyze: bool
    :param hot_start: bool
    :param storage_file_path: str
    :param export: bool
    :param output_dir: str
    :param export_file_path: str
    :param disp: bool
    """
    init_main_controller(sleep, cluster_id, profile, path_length, spines, param_file_path, mech_file_path,
                         neurotree_file_path, neurotree_index, param_gen, get_features, features_modules,
                         objectives_modules, group_sizes, storage_file_path, output_dir, export_file_path, pop_size)
    main_ctxt.c[:].execute('from parallel_optimize_main import *', block=True)
    if sleep:
        time.sleep(120.)
    main_ctxt.c[:].apply_sync(setup_modules, main_ctxt.get_features, main_ctxt.features_modules,
                              main_ctxt.objectives_modules, main_ctxt.group_sizes, True, main_ctxt.param_names,
                              main_ctxt.mech_file_path, main_ctxt.neurotree_dict, spines, output_dir)
    global storage
    global x
    if not analyze:
        global this_param_gen
        if hot_start:
            this_param_gen = main_ctxt.param_gen_func(pop_size=pop_size, x0=param_dict_to_array(main_ctxt.x0,
                                                                                               main_ctxt.param_names),
                                                              bounds=main_ctxt.bounds, wrap_bounds=wrap_bounds,
                                                              seed=seed, max_iter=max_iter,
                                                              adaptive_step_factor=adaptive_step_factor,
                                                              survival_rate=survival_rate, disp=disp, hot_start=storage_file_path)
        else:
            this_param_gen = main_ctxt.param_gen_func(param_names=main_ctxt.param_names,
                                                              feature_names=main_ctxt.feature_names,
                                                              objective_names=main_ctxt.objective_names, pop_size=pop_size,
                                                              x0=param_dict_to_array(main_ctxt.x0, main_ctxt.param_names),
                                                              bounds=main_ctxt.bounds, wrap_bounds=wrap_bounds,
                                                              seed=seed, max_iter=max_iter, path_length=path_length,
                                                              initial_step_size=initial_step_size,
                                                              adaptive_step_factor=adaptive_step_factor,
                                                              survival_rate=survival_rate, disp=disp)
        run_optimization(disp)
        storage = this_param_gen.storage
        best_ind = storage.get_best(1, 'last')[0]
        x = param_array_to_dict(best_ind.x, storage.param_names)
        if disp:
            print 'Multi-Objective Optimization: Best params:'
            print x
    elif os.path.isfile(storage_file_path):
        storage = PopulationStorage(file_path=storage_file_path)
        print 'Analysis mode: history loaded from path: %s' % storage_file_path
        best_ind = storage.get_best(1, 'last')[0]
        x = param_array_to_dict(best_ind.x, storage.param_names)
        if disp:
            print 'Multi-Objective Optimization: Best params:'
            print x
    else:
        print 'Analysis mode: history not loaded'
        x = main_ctxt.x0
        if disp:
            print 'Multi-Objective Optimization: Loaded params:'
            print x
    sys.stdout.flush()
    if export:
        return export_traces(x, main_ctxt.group_sizes, export_file_path=main_ctxt.export_file_path, disp=disp)

def init_main_controller(sleep, cluster_id, profile, path_length, spines, param_file_path, mech_file_path,
                         neurotree_file_path, neurotree_index, param_gen, get_features, features_modules,
                         objectives_modules, group_sizes, storage_file_path, output_dir, export_file_path, pop_size):
    """

    :param sleep:
    :param cluster_id:
    :param profile:
    :param path_length:
    :param spines:
    :param param_file_path:
    :param mech_file_path:
    :param neurotree_file_path:
    :param neurotree_index:
    :param param_gen:
    :param get_features:
    :param features_modules:
    :param objectives_modules:
    :param group_sizes:
    :param storage_file_path:
    :param output_dir:
    :param export_file_path:
    :param pop_size:
    :return:
    """
    if sleep:
        time.sleep(300.)
    if cluster_id is not None:
        c = Client(cluster_id=cluster_id, profile=profile)
    else:
        c = Client(profile=profile)
    num_procs = len(c)
    main_ctxt.c = c
    main_ctxt.num_procs = num_procs

    if param_file_path is not None:
        params_dict = read_from_yaml(param_file_path)
        param_names = params_dict['param_names']
        x0 = params_dict['x0']
        bounds = [params_dict['bounds'][key] for key in param_names]
        feature_names = params_dict['feature_names']
        objective_names = params_dict['objective_names']
        target_val = params_dict['target_val']
        target_range = params_dict['target_range']
        optimization_title = params_dict['optimization_title']

        if param_gen is None:
            param_gen = params_dict['param_gen']
        if mech_file_path is None:
            mech_file_path = params_dict['mech_file_path']
        if neurotree_file_path is None:
            neurotree_file_path = params_dict['neurotree_file_path']
        if not get_features:
            get_features = params_dict['get_features']
        if not features_modules:
            features_modules = params_dict['features_modules']
        if not objectives_modules:
            objectives_modules = params_dict['objectives_modules']
        if not group_sizes:
            group_sizes = params_dict['group_sizes']
        if storage_file_path is None:
            storage_file_path = '%s/%s_%s_%s_optimization_history.hdf5' % \
                                (output_dir, datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)
        if export_file_path is None:
            export_file_path = '%s/%s_%s_%s_optimization_exported_traces.hdf5' % \
                               (output_dir, datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)
    else:
        raise NameError('A param_file_path containing default paramters must be provided.')
    main_ctxt.update(locals())
    main_ctxt.update({'x0': x0, 'param_names': param_names, 'bounds': bounds, 'feature_names': feature_names,
                    'objective_names': objective_names, 'target_val': target_val, 'target_range': target_range,
                    'path_length': path_length, 'spines': spines, 'pop_size': pop_size})

    if param_gen not in globals():
        raise NameError('Multi-Objective Optimization: %s has not been imported, or is not a valid class of parameter '
                        'generator.' % param_gen)
    param_gen_func = globals()[param_gen] # The variable 'param_gen_func' points to the actual generator function, while
                                          # param_gen points to the string name of the generator
    main_ctxt.param_gen_func = param_gen_func
    param_gen_func_name = param_gen

    neurotree_dict = read_from_pkl(neurotree_file_path)[neurotree_index]
    main_ctxt.neurotree_dict = neurotree_dict

    new_group_sizes = [group_sizes[i] for i in range(len(group_sizes))]
    blocks = range(len(group_sizes))
    un_utilized = range(len(group_sizes))
    for ind in range(len(group_sizes)):
        if new_group_sizes[ind] > num_procs:
            new_group_sizes[ind] = num_procs
            print 'Multi-Objective Optimization: group_sizes index %i adjusted to not exceed num_processes: %i' % (ind, num_procs)
        un_utilized[ind] = num_procs % new_group_sizes[ind]
        if un_utilized[ind] > 0:
            print 'Multi-Objective Optimization: index %i has %i unutilized processes' % (ind, un_utilized[ind])
        blocks[ind] = pop_size / (num_procs / new_group_sizes[ind]) * (len(get_features) / new_group_sizes[ind])
        if blocks[ind] * (num_procs / new_group_sizes[ind]) < pop_size:
            blocks[ind] += 1
    group_sizes = new_group_sizes
    main_ctxt.group_sizes = new_group_sizes

    print 'Multi-Objective Optimization %s: Generator: %s; Total processes: %i; Population size: %i; Group sizes: %s; ' \
          'Feature calculator: %s; Objective calculator: %s; Blocks / generation: %s' % \
          (optimization_title, param_gen_func_name, num_procs, pop_size, (','.join(str(x) for x in group_sizes)),
           (','.join(func for func in get_features)), (','.join(obj for obj in objectives_modules)),
           (','.join(str(x) for x in blocks)))
    sys.stdout.flush()
    setup_modules(get_features, features_modules, objectives_modules, group_sizes)

def setup_modules(get_features, features_modules, objectives_modules, group_sizes, config_engine=False,
                  param_names=None, mech_file_path=None, neurotree_dict=None, spines=None, output_dir=None):
    feat_module_refs = []
    feat_module_dict = {}
    for module in set(features_modules):
        m = importlib.import_module(module)
        feat_module_dict[module] = m
        if config_engine:
            global rec_filepath
            rec_filepath = output_dir+'/sim_output' + datetime.datetime.today().strftime('%m%d%Y%H%M') + \
                           '_pid' + str(os.getpid()) + '.hdf5'
            m.config_engine(param_names, mech_file_path, neurotree_dict, spines, rec_filepath)
    for module in features_modules:
        feat_module_refs.append(feat_module_dict[module])
    if len(get_features) != len(group_sizes):
        raise NameError('Number of arguments in get_features does not match number of arguments in group_sizes.')
    if len(get_features) != len(feat_module_refs):
        raise NameError('Number of arguments in get_features does not match number of imported modules.')

    get_features_func = [getattr(feat_module_refs[i], get_features[i]) for i in range(len(get_features))]
    main_ctxt.get_features_func = get_features_func
    for i in range(len(get_features_func)):
        if not callable(get_features_func[i]):
            raise NameError('Multi-Objective Optimization: get_features: %s for module %s is not a callable function.'
                            % (get_features[i], features_modules[i]))
    obj_module_refs = []
    for module in objectives_modules:
        m = sys.modules[module]
        obj_module_refs.append(m)
    get_objectives_func = [getattr(obj_module_refs[i], 'get_objectives') for i in range(len(obj_module_refs))]
    main_ctxt.get_objectives_func = get_objectives_func
    for i in range(len(get_objectives_func)):
        if not callable(get_objectives_func[i]):
            raise NameError('Multi-Objective Optimization: get_objectives for module %s is not a callable function.'
                            % objectives_modules[i])
    sys.stdout.flush()

def run_optimization(disp):
    """

    :param group_sizes: tuple of int
    :param path_length: int
    :param disp: bool
    :return:
    """
    for ind, generation in enumerate(this_param_gen()):
        if (ind > 0) and (ind % main_ctxt.path_length == 0):
            this_param_gen.storage.save(main_ctxt.storage_file_path, n=main_ctxt.path_length)
        features, objectives = get_all_features(generation, group_sizes=main_ctxt.group_sizes, disp=disp)
        this_param_gen.update_population(features, objectives)
    this_param_gen.storage.save(main_ctxt.storage_file_path, n=main_ctxt.path_length)

def get_all_features(generation, group_sizes=(1, 10), disp=False, export=False):
    """

    :param generation: list of arr
    :param group_sizes: list of int
    :param disp: bool
    :param export: False (for exporting voltage traces)
    :return: tuple of list of dict
    """
    pop_ids = range(len(generation))
    results = []
    curr_generation = {pop_id: generation[pop_id] for pop_id in pop_ids}
    final_features = {pop_id: {} for pop_id in pop_ids}
    for ind in range(len(main_ctxt.get_features_func)):
        next_generation = {}
        this_group_size = min(len(main_ctxt.c), group_sizes[ind])
        usable_procs = main_ctxt.num_procs - (main_ctxt.num_procs % this_group_size)
        client_ranges = [range(start, start + this_group_size) for start in range(0, usable_procs, this_group_size)]
        feature_function = main_ctxt.get_features_func[ind]
        indivs = [{'pop_id': pop_id, 'x': curr_generation[pop_id],
                   'features': final_features[pop_id]} for pop_id in curr_generation.keys()]
        while len(indivs) > 0 or len(results) > 0:
            num_groups = min(len(client_ranges), len(indivs))
            if num_groups > 0:
                results.extend(map(feature_function, [indivs.pop(0) for i in range(num_groups)], [main_ctxt.c] * num_groups,
                                   [client_ranges.pop(0) for i in range(num_groups)], [export] * num_groups))
            if np.any([this_result['async_result'].ready() for this_result in results]):
                for this_result in results:
                    if this_result['async_result'].ready():
                        client_ranges.append(this_result['client_range'])
                        if disp:
                            flush_engine_buffer(this_result['async_result'])
                        get_result = this_result['async_result'].get()
                        if get_result is None:
                            final_features[this_result['pop_id']] = None
                        else:
                            next_generation[this_result['pop_id']] = generation[this_result['pop_id']]
                            if 'filter_features' in this_result.keys():
                                filter_features_func = this_result['filter_features']
                                if not callable(filter_features_func):
                                    raise NameError(
                                        'Multi-Objective Optimization: filter_features function %s not callable' %
                                        filter_features_func)
                                new_features = filter_features_func(get_result, final_features[this_result['pop_id']])
                            else:
                                new_features = {key: value for result_dict in get_result for key, value in
                                                result_dict.iteritems()}
                            final_features[this_result['pop_id']].update(new_features)
                            if disp:
                                print 'Individual: %i, computing %s took %.2f s' % (this_result['pop_id'],
                                                                                    main_ctxt.get_features[ind],
                                                                                    this_result['async_result'].wall_time)
                        results.remove(this_result)
                    sys.stdout.flush()
            else:
                time.sleep(1.)
        curr_generation = next_generation
    features = [final_features[pop_id] for pop_id in pop_ids]
    objectives = []
    for i, this_features in enumerate(features):
        this_objectives = {}
        for j, objective_function in enumerate(main_ctxt.get_objectives_func):
            new_features, new_objectives = objective_function(this_features, main_ctxt.objective_names,
                                                              main_ctxt.target_val, main_ctxt.target_range)
            features[i] = new_features
            this_objectives.update(new_objectives)
        objectives.append(this_objectives)
    return features, objectives

def export_traces(x, group_sizes, export_file_path=None, discard=True, disp=False):
    """
    Run simulations on the engines with the given parameter values, have the engines export their results to .hdf5,
    and then read in and plot the results.

    :param x: dict
    :param group_sizes: tuple of int
    :param export_file_path: str
    :param discard: bool
    :param disp: bool
    """
    x_arr = param_dict_to_array(x, main_ctxt.param_names)
    exported_features, exported_objectives = get_all_features([x_arr], group_sizes=group_sizes, disp=disp, export=True)
    rec_file_path_list = [filepath for filepath in main_ctxt.c[:]['rec_filepath']
                          if os.path.isfile(filepath)]
    combine_hdf5_file_paths(rec_file_path_list, export_file_path)
    if discard:
        for rec_file_path in rec_file_path_list:
            os.remove(rec_file_path)
    print 'Multi-Objective Optimization: Exported traces to %s' % export_file_path
    sys.stdout.flush()
    return exported_features, exported_objectives

def plot_traces(export_file_path):
    with h5py.File(export_file_path, 'r') as f:
        for trial in f.itervalues():
            amplitude = trial.attrs['amp']
            fig, axes = plt.subplots(1)
            for rec in trial['rec'].itervalues():
                axes.plot(trial['time'], rec, label=rec.attrs['description'])
            axes.legend(loc='best', frameon=False, framealpha=0.5)
            axes.set_xlabel('Time (ms)')
            axes.set_ylabel('Vm (mV)')
            # axes.set_title('Optimize %s: I_inj amplitude %.2f' % (trial.attrs['description'], amplitude))
            clean_axes(axes)
            fig.tight_layout()
        plt.show()
        plt.close()


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1,sys.argv)+1):])

