__author__ = 'Grace Ng'
import click
from ipyparallel import interactive
from ipyparallel import Client
# from IPython.display import clear_output
from specify_cells3 import *
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
@click.option("--export-file-path", type=str, default=None)
@click.option("--disp", is_flag=True)
def main(cluster_id, profile, spines, mech_file_path, neurotree_file_path, neurotree_index, param_file_path,
         param_gen, get_features, features_modules, objectives_modules, group_sizes, pop_size, wrap_bounds, seed,
         max_iter, path_length, initial_step_size, adaptive_step_factor, survival_rate, sleep, analyze, hot_start,
         storage_file_path, export, export_file_path, disp):
    """

    :param cluster_id: str
    :param profile: str
    :param spines: bool
    :param mech_file_path: str (path)
    :param neurotree_file_path: str (path)
    :param neurotree_index: int
    :param param_file_path: str (path)
    :param param_gen: str (must refer to callable in globals())
    :param get_features: tuple of str (must refer to callable in globals())
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
    :param export_file_path: str
    :param disp: bool
    """
    global c
    if sleep:
        time.sleep(300.)
    if cluster_id is not None:
        c = Client(cluster_id=cluster_id, profile=profile)
    else:
        c = Client(profile=profile)

    global num_procs
    num_procs = len(c)

    global x0
    global param_names
    global bounds
    global feature_names
    global objective_names
    global target_val
    global target_range
    globals()['path_length'] = path_length
    globals()['spines'] = spines
    global neurotree_dict
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

        given_param_gen = params_dict['param_gen']
        given_mech_file_path = params_dict['mech_file_path']
        given_neurotree_file_path = params_dict['neurotree_file_path']
        given_get_features = params_dict['get_features']
        given_features_modules = params_dict['features_modules']
        given_objectives_modules = params_dict['objectives_modules']
        given_group_sizes = params_dict['group_sizes']

    if mech_file_path is None:
        mech_file_path = given_mech_file_path
    globals()['mech_file_path'] = mech_file_path
    if neurotree_file_path is None:
        neurotree_file_path = given_neurotree_file_path
    neurotree_dict = read_from_pkl(neurotree_file_path)[neurotree_index]
    if param_gen is None:
        param_gen = given_param_gen
    if not get_features:
        get_features = given_get_features
    if not features_modules:
        features_modules = given_features_modules
    if not objectives_modules:
        objectives_modules = given_objectives_modules
    if not group_sizes:
        group_sizes = given_group_sizes
    if storage_file_path is None:
        storage_file_path = 'data/%s_%s_%s_optimization_history.hdf5' % \
                       (datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)
    globals()['storage_file_path'] = storage_file_path
    if export_file_path is None:
        export_file_path = 'data/%s_%s_%s_optimization_exported_traces.hdf5' % \
                           (datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)
    globals()['export_file_path'] = export_file_path

    if param_gen not in globals():
        raise NameError('Multi-Objective Optimization: %s has not been imported, or is not a valid class of parameter '
                        'generator.' % param_gen)
    param_gen_func = globals()[param_gen]
    globals()['param_gen_func'] = param_gen_func #The global variable 'param_gen_func' points to the actual function
    global param_gen_func_name
    param_gen_func_name = param_gen #The global 'param_gen_func_name' points to the string with the name of the param_gen

    globals()['pop_size'] = pop_size
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
        blocks[ind] = pop_size / (num_procs / new_group_sizes[ind])
        if blocks[ind] * (num_procs / new_group_sizes[ind]) < pop_size:
            blocks[ind] += 1
    globals()['group_sizes'] = new_group_sizes
    group_sizes = new_group_sizes

    setup_module_functions(get_features, features_modules, objectives_modules, group_sizes)
    c[:].execute('from parallel_optimize_main import *', block=True)
    if sleep:
        time.sleep(120.)
    c[:].map_sync(setup_module_functions, [get_features] * num_procs, [features_modules] * num_procs,
                  [objectives_modules] * num_procs, [group_sizes] * num_procs)
    print 'Multi-Objective Optimization %s: Generator: %s; Total processes: %i; Population size: %i; Group sizes: %s; ' \
          'Feature calculator: %s; Objective calculator: %s; Blocks / generation: %s' % \
          (optimization_title, param_gen_func_name, num_procs, pop_size, (','.join(str(x) for x in group_sizes)),
           (','.join(func for func in get_features_names)), (','.join(obj for obj in objectives_modules)),
           (','.join(str(x) for x in blocks)))
    sys.stdout.flush()

    global storage
    global x
    if not analyze:
        if hot_start:
            globals()['param_gen'] = param_gen_func(pop_size, x0=param_dict_to_array(x0, param_names), bounds=bounds,
                                                    wrap_bounds=wrap_bounds, seed=seed, max_iter=max_iter,
                                                    adaptive_step_factor=adaptive_step_factor,
                                                    survival_rate=survival_rate, disp=disp, hot_start=storage_file_path)
        else:
            globals()['param_gen'] = param_gen_func(param_names, feature_names, objective_names, pop_size,
                                                    x0=param_dict_to_array(x0, param_names), bounds=bounds,
                                                    wrap_bounds=wrap_bounds, seed=seed, max_iter=max_iter,
                                                    path_length=path_length, initial_step_size=initial_step_size,
                                                    adaptive_step_factor=adaptive_step_factor,
                                                    survival_rate=survival_rate, disp=disp)
        this_param_gen = globals()['param_gen']
        run_optimization(group_sizes, path_length, disp)
        storage = this_param_gen.storage
        best_ind = storage.get_best(1, 'last')[0]
        x = param_array_to_dict(best_ind.x, param_names)
        if disp:
            print 'Multi-Objective Optimization: Best params:'
            print x
    elif os.path.isfile(storage_file_path):
        storage = PopulationStorage(file_path=storage_file_path)
        print 'Analysis mode: history loaded from path: %s' % storage_file_path
        best_ind = storage.get_best(1, 'last')[0]
        x = param_array_to_dict(best_ind.x, param_names)
        if disp:
            print 'Multi-Objective Optimization: Best params:'
            print x
    else:
        print 'Analysis mode: history not loaded'
        x = x0
        if disp:
            print 'Multi-Objective Optimization: Loaded params:'
            print x
    sys.stdout.flush()
    if export:
        export_traces(x, group_sizes, export_file_path=export_file_path, disp=disp)

@interactive
def setup_module_functions(get_features, features_modules, objectives_modules, group_sizes):
    global feat_module_refs
    feat_module_refs = []
    for module in features_modules:
        m = importlib.import_module(module)
        feat_module_refs.append(m)
    if len(get_features) != len(group_sizes):
        raise NameError('Number of arguments in get_features does not match number of arguments in group_sizes.')
    if len(get_features) != len(feat_module_refs):
        raise NameError('Number of arguments in get_features does not match number of imported modules.')
    get_features_func = [getattr(feat_module_refs[i], get_features[i]) for i in range(len(get_features))]
    global get_features_names
    get_features_names = [get_features[i] for i in range(len(get_features))]
    for i in range(len(get_features_func)):
        if not callable(get_features_func[i]):
            raise NameError('Multi-Objective Optimization: get_features: %s for module %s is not a callable function.'
                            % (get_features_names[i], features_modules[i]))
    global obj_module_refs
    obj_module_refs = []
    for module in objectives_modules:
        m = importlib.import_module(module)
        obj_module_refs.append(m)
    get_objectives_func = [getattr(obj_module_refs[i], 'get_objectives') for i in range(len(obj_module_refs))]
    for i in range(len(get_objectives_func)):
        if not callable(get_objectives_func[i]):
            raise NameError('Multi-Objective Optimization: get_objectives for module %s is not a callable function.'
                            % objectives_modules[i])
    globals()['get_features_func'] = get_features_func
    globals()['get_objectives_func'] = get_objectives_func
    global prev_job_type
    prev_job_type = None

@interactive
def run_optimization(group_sizes, path_length, disp):
    """

    :param group_sizes: tuple of int
    :param path_length: int
    :param disp: bool
    :return:
    """
    for ind, generation in enumerate(param_gen()):
        if (ind > 0) and (ind % path_length == 0):
            param_gen.storage.save(storage_file_path, n=path_length)
        features, objectives = compute_features(generation, group_sizes=group_sizes, disp=disp)
        param_gen.update_population(features, objectives)
    param_gen.storage.save(storage_file_path, n=path_length)

@interactive
def compute_features(generation, group_sizes=(1, 10), disp=False, export=False):
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
    for ind in range(len(get_features_func)):
        next_generation = {}
        this_group_size = min(len(c), group_sizes[ind])
        usable_procs = num_procs - (num_procs % this_group_size)
        client_ranges = [range(start, start + this_group_size) for start in range(0, usable_procs, this_group_size)]
        feature_function = get_features_func[ind]
        indivs = [{'pop_id': pop_id, 'x': curr_generation[pop_id],
                   'features': final_features[pop_id]} for pop_id in curr_generation.keys()]
        while len(indivs) > 0 or len(results) > 0:
            num_groups = min(len(client_ranges), len(indivs))
            if num_groups > 0:
                results.extend(map(feature_function, [indivs.pop(0) for i in range(num_groups)], [c] * num_groups,
                                   [client_ranges.pop(0) for i in range(num_groups)], [param_names] * num_groups,
                                   [mech_file_path] * num_groups, [neurotree_dict] * num_groups, [spines] * num_groups,
                                   [ind] * num_groups,
                                   [feat_module_refs[ind] for i in range(num_groups)], [export] * num_groups))
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
                                new_features = filter_features_func(get_result, feat_module_refs[ind],
                                                                    final_features[this_result['pop_id']])
                            else:
                                new_features = {key: value for result_dict in get_result for key, value in
                                                result_dict.iteritems()}
                            final_features[this_result['pop_id']].update(new_features)
                            if disp:
                                print 'Individual: %i, computing %s took %.2f s' % (this_result['pop_id'],
                                                                                    get_features_names[ind],
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
        for j, objective_function in enumerate(get_objectives_func):
            new_features, new_objectives = objective_function(obj_module_refs[j], this_features, objective_names,
                                                              target_val, target_range)
            features[i] = new_features
            this_objectives.update(new_objectives)
        objectives.append(this_objectives)
    return features, objectives

@interactive
def param_array_to_dict(x, param_names):
    """

    :param x: arr
    :param param_names: list
    :return:
    """
    return {param_name: x[ind] for ind, param_name in enumerate(param_names)}


@interactive
def param_dict_to_array(x_dict, param_names):
    """

    :param x_dict: dict
    :param param_names: list
    :return:
    """
    return np.array([x_dict[param_name] for param_name in param_names])

@interactive
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
    x_arr = param_dict_to_array(x, param_names)
    global exported_features
    global exported_objectives
    exported_features, exported_objectives = compute_features([x_arr], group_sizes=group_sizes, disp=disp, export=True)
    rec_file_path_list = [data_dir + filename + '.hdf5' for filename in c[:]['rec_filename']
                          if os.path.isfile(data_dir + filename + '.hdf5')]
    combine_hdf5_file_paths(rec_file_path_list, export_file_path)
    if discard:
        for rec_file_path in rec_file_path_list:
            os.remove(rec_file_path)
    print 'Multi-Objective Optimization: Exported traces to %s' % export_file_path
    sys.stdout.flush()

@interactive
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
            axes.set_title('Optimize %s: I_inj amplitude %.2f' % (trial.attrs['description'], amplitude))
            clean_axes(axes)
            fig.tight_layout()
        plt.show()
        plt.close()


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1,sys.argv)+1):])

