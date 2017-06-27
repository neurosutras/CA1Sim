__author__ = 'Grace Ng'
import click
from ipyparallel import interactive
from ipyparallel import Client
# from IPython.display import clear_output
import sys

try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass


script_filename = 'parallel_optimize_leak.py'

equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
#dt = 0.02
dt = 0.002
amp = 0.3
th_dvdt = 10.
v_init = -77.
v_active = -77.
i_holding = {'soma': 0., 'dend': 0., 'distal_dend': 0.}
soma_ek = -77.

default_mech_file_path = 'data/042717 GC optimizing spike stability.pkl'
default_neurotree_file_path = 'morphologies/121516_DGC_trees.pkl'
default_param_gen = 'BGen'
default_get_features = 'get_Rinp_features'
default_get_objectives = 'get_pas_objectives'

default_x0_dict = {'soma.g_pas': 1.050E-10, 'dend.g_pas slope': 1.058E-08, 'dend.g_pas tau': 3.886E+01}  # Error: 4.187E-09
default_bounds_dict = {'soma.g_pas': (1.0E-18, 1.0E-6), 'dend.g_pas slope': (1.0E-12, 1.0E-4),
                  'dend.g_pas tau': (25., 400.)}
default_feature_names = ['soma R_inp', 'dend R_inp', 'distal_dend R_inp']
default_objective_names = ['soma R_inp', 'dend R_inp', 'distal_dend R_inp']
default_target_val = {'soma R_inp': 295., 'dend R_inp': 375.}
default_target_range = {'soma R_inp': 0.5, 'dend R_inp': 1.}
default_optimization_title = 'leak'

# Load defaults from a file containing many parameters. Use write_to_file function in function_lib to generate this pkl file
# Note: click only recognizes this as a path if the string is copied as is into the command line; cannot type a variable
# name that stores this string
default_param_file_path = None #'data/leak_default_param_file.pkl'

@click.command()
@click.option("--cluster-id", type=str, default=None)
@click.option("--spines", is_flag=True)
@click.option("--mech-file-path", type=click.Path(exists=False, file_okay=True, dir_okay=False), default=None)
@click.option("--neurotree-file-path", type=click.Path(exists=False, file_okay=True, dir_okay=False), default=None)
@click.option("--neurotree-index", type=int, default=0)
@click.option("--param-file-path", type=click.Path(exists=False, file_okay=True, dir_okay=False), default=None)
@click.option("--param-gen", type=str, default=None)
@click.option("--get-features", type=str, default=None)
@click.option("--get-objectives", type=str, default=None)
@click.option("--group-size", type=int, default=3)
@click.option("--pop-size", type=int, default=100)
@click.option("--seed", type=int, default=None)
@click.option("--max-iter", type=int, default=30)
@click.option("--max-gens", type=int, default=None)
@click.option("--path-length", type=int, default=1)
@click.option("--adaptive-step-factor", type=float, default=0.9)
@click.option("--niter-success", type=int, default=None)
@click.option("--survival-rate", type=float, default=0.2)
@click.option("--disp", is_flag=True)
def main(cluster_id, spines, mech_file_path, neurotree_file_path, neurotree_index, param_file_path, param_gen,
         get_features, get_objectives, group_size, pop_size, seed, max_iter, max_gens, path_length,
         adaptive_step_factor, niter_success, survival_rate, disp):
    """

    :param cluster_id: str
    :param spines: bool
    :param mech_file_path: str (path)
    :param neurotree_file_path: str (path)
    :param neurotree_index: int
    :param param_file_path: str (path)
    :param param_gen: str (must refer to callable in globals())
    :param get_features: str (must refer to callable in globals())
    :param get_objectives: str (must refer to callable in globals())
    :param group_size: int
    :param pop_size: int
    :param seed: int
    :param max_iter: int
    :param max_gens: int
    :param path_length: int
    :param adaptive_step_factor: float in [0., 1.]
    :param niter_success: int
    :param survival_rate: float
    :param disp: bool
    """
    global c

    if cluster_id is not None:
        c = Client(cluster_id=cluster_id)
    else:
        c = Client()

    global num_procs
    num_procs = len(c)

    if mech_file_path is None:
        mech_file_path = default_mech_file_path

    if neurotree_file_path is None:
        neurotree_file_path = default_neurotree_file_path

    if param_file_path is None:
        param_file_path = default_param_file_path

    global x0
    global param_names
    global bounds
    global feature_names
    global objective_names
    global target_val
    global target_range
    global optimization_title

    param_names = default_x0_dict.keys()
    x0 = [default_x0_dict[key] for key in param_names]
    bounds = [default_bounds_dict[key] for key in param_names]
    feature_names = default_feature_names
    objective_names = default_objective_names
    target_val = default_target_val
    target_range = default_target_range
    optimization_title = default_optimization_title

    globals()['path_length'] = path_length

    if param_gen is None:
        param_gen = default_param_gen
    if param_gen not in globals():
        raise NameError('Multi-Objective Optimization: %s has not been imported, or is not a valid class of parameter '
                        'generator.' % param_gen)
    """
    global history_filename
    history_filename = '%s %s %s optimization history' % \
                       (datetime.datetime.today().strftime('%m%d%Y%H%M'), optimization_title, param_gen)
    """
    if get_features is None:
        get_features = default_get_features
    if get_features not in globals() or not callable(globals()[get_features]):
        raise NameError('Multi-Objective Optimization: get_features: %s has not been imported, or is not a callable '
                        'function.' % get_features)

    if get_objectives is None:
        get_objectives = default_get_objectives
    if get_objectives not in globals() or not callable(globals()[get_objectives]):
        raise NameError('Multi-Objective Optimization: get_objectives: %s has not been imported, or is not a callable '
                        'function.' % get_objectives)

    globals()['group_size'] = group_size
    globals()['pop_size'] = pop_size

    if group_size > num_procs:
        group_size = num_procs
        print 'Multi-Objective Optimization: group_size adjusted to not exceed num_processes: %i' % num_procs
    un_utilized = num_procs % group_size
    blocks = pop_size / (num_procs / group_size)
    if blocks * (num_procs / group_size) < pop_size:
        blocks += 1

    print 'Multi-Objective Optimization: %s; Total processes: %i; Population size: %i; Group size: %i; ' \
          'Feature calculator: %s; Objective calculator: %s; Blocks / generation: %i' % \
          (param_gen, num_procs, pop_size, group_size, get_features, get_objectives, blocks)
    if un_utilized > 0:
        print 'Multi-Objective Optimization: %i processes are unutilized' % un_utilized
    sys.stdout.flush()

def list_find (f, lst):
    i = 0
    for x in lst:
        if f(x):
            return i
        else:
            i += 1
    return None

if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1,sys.argv)+1):])