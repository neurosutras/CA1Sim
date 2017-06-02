__author__ = 'Grace Ng'
from ipyparallel import interactive
from ipyparallel import Client
import click
from IPython.display import clear_output
from specify_cells3 import *
import time
import os
import sys
import pprint
from moopgen import *

"""
Aims for spike initiation at initial segment by increasing nax density and decreasing activation V1/2 relative to soma,
axon_hill, and axon compartments. Extend linear kap gradient into basals and obliques, aim for 60% spike attenuation
at bifurcation of trunk and tuft.

Hierarchical optimization:
I) optimize g_pas for target rinp at soma, trunk bifurcation, and tuft bifurcation [without h].
II) optimize ghbar_h for target rinp at soma, trunk bifurcation, and tuft bifurcation, while also optimizing for v_rest
offset between soma and tuft, and EPSP shape changes between proximal and distal synapses measured at the soma.
III) optimize gbar_nax/nas/sh/sha, gkabar_kap/d, gkdrbar for target na spike threshold, AHP amp, and vm stability

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass


class History(object):
    def __init__(self):
        """

        """
        self.xlabels = []
        self.x_values = []
        self.error_values = []
        self.Rinp_values = {}

    def report_best(self):
        """
        Report the input parameters and output values with the lowest error.
        """
        lowest_Err = min(self.error_values)
        index = self.error_values.index(lowest_Err)
        best_x = self.x_values[index]
        best_Rinp_values = {section: self.Rinp_values[section][index] for section in self.Rinp_values}
        formatted_x = '[' + ', '.join(['%.3E' % xi for xi in best_x]) + ']'
        print 'best x: %s' % formatted_x
        print 'lowest Err: %.3E' % lowest_Err
        print 'Rinp:', ['%s: %.1f' % (section, Rinp) for (section, Rinp) in best_Rinp_values.iteritems()]
        return best_x

    def export_to_pkl(self, hist_filename):
        """
        Save the history to .pkl
        :param hist_filename: str
        """
        saved_history = {'xlabels': self.xlabels, 'x_values': self.x_values, 'error_values': self.error_values,
                         'Rinp_values': self.Rinp_values}
        write_to_pkl(data_dir+hist_filename+'.pkl', saved_history)

    def import_from_pkl(self, hist_filename):
        """
        Update a history object with data from a .pkl file
        :param hist_filename: str
        """
        previous_history = read_from_pkl(data_dir+hist_filename +'.pkl')
        self.xlabels = previous_history['xlabels']
        self.x_values = previous_history['x_values']
        self.error_values = previous_history['error_values']
        self.Rinp_values = previous_history['Rinp_values']

    def plot(self):
        """
        Remember to : also plot each value in x against error, and against input resistance
        """
        num_x_param = len(self.xlabels)
        num_plot_rows = math.floor(math.sqrt(num_x_param))
        print(num_plot_rows)
        num_plot_cols = math.ceil(num_x_param/num_plot_rows)
        print(num_plot_cols)

        #plot x-values against error
        plt.figure(1)
        for i, x_param in enumerate(self.xlabels):
            plt.subplot(num_plot_rows, num_plot_cols, i+1)
            x_param_vals = [x_val[i] for x_val in self.x_values]
            range_param_vals = max(x_param_vals) - min(x_param_vals)
            plt.scatter(x_param_vals, self.error_values)
            plt.xlim((min(x_param_vals)-0.1*range_param_vals, max(x_param_vals)+0.1*range_param_vals))
            plt.xlabel(x_param)
            plt.ylabel("Error values")
        #plt.show()
        #plt.close()

        plt.figure(2)
        colors = ['r', 'g', 'b', 'gray', 'darkviolet', 'goldenrod']
        #plot x-values against input resistance
        for i, x_param in enumerate(self.xlabels):
            plt.subplot(num_plot_rows, num_plot_cols, i+1)
            x_param_vals = [x_val[i] for x_val in self.x_values]
            range_param_vals = max(x_param_vals) - min(x_param_vals)
            for j, section in enumerate(self.Rinp_values):
                plt.scatter([x_param_vals], self.Rinp_values[section], label=section, color = colors[j])
            plt.xlim((min(x_param_vals) - 0.1 * range_param_vals, max(x_param_vals) + 0.1 * range_param_vals))
            plt.xlabel(x_param)
            plt.ylabel("Rinp values")
            plt.legend(loc='upper right', scatterpoints = 1, frameon=False, framealpha=0.5)
        plt.show()
        plt.close()


script_filename = 'parallel_optimize_leak.py'
history_filename = '030517 leak optimization history'

equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.02
amp = 0.3
th_dvdt = 10.
v_init = -77.
v_active = -77.
i_holding = {'soma': 0., 'dend': 0., 'distal_dend': 0.}
soma_ek = -77.

default_mech_file_path = data_dir + '042717 GC optimizing spike stability.pkl'
default_neurotree_file_path = morph_dir + '121516_DGC_trees.pkl'
default_param_gen = 'BGen'
default_get_features = 'get_Rinp_features'
default_get_objectives = 'get_objectives'

default_x0_dict = {'soma.g_pas': 1.050E-10, 'dend.g_pas slope': 1.058E-08, 'dend.g_pas tau': 3.886E+01}  # Error: 4.187E-09
default_bounds_dict = {'soma.g_pas': (1.0E-18, 1.0E-6), 'dend.g_pas slope': (1.0E-12, 1.0E-4),
                  'dend.g_pas tau': (25., 400.)}
default_feature_names = ['soma R_inp', 'dend R_inp', 'distal_dend R_inp']
default_objective_names = ['soma R_inp', 'dend R_inp', 'distal_dend R_inp']
# we should load defaults from a file if we're going to be running optimizations with many more parameters
default_param_file_path = None


@click.command()
@click.option("--cluster-id", type=str, default=None)
@click.option("--spines", type=bool, default=False)
@click.option("--mech-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None)
@click.option("--neurotree-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None)
@click.option("--neurotree-index", type=int, default=0)
@click.option("--param-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None)
@click.option("--adaptive-interval", type=int, default=1)
@click.option("--group-size", type=int, default=1)
@click.option("--pop-size", type=int, default=50)
@click.option("--param-gen", type=str, default=None)
@click.option("--get_features", type=str, default=None)
@click.option("--get_objectives", type=str, default=None)
def main(cluster_id, spines, mech_file_path, neurotree_file_path, neurotree_index, param_file_path, adaptive_interval,
         group_size, pop_size, param_gen, get_features, get_objectives):
    """

    :param cluster_id: str
    :param spines: bool
    :param mech_file_path: str (path)
    :param neurotree_file_path: str (path)
    :param neurotree_index: int
    :param param_file_path: str (path)
    :param adaptive_interval: int
    :param group_size: int
    :param pop_size: int
    :param param_gen: str (name of callable)
    :param get_features: str (name of callable)
    :param get_objectives: str (name of callable)
    """
    global c

    if cluster_id is not None:
        c = Client(cluster_id=cluster_id)
    else:
        c = Client()

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

    if param_file_path is not None:
        params_dict = read_from_pkl(param_file_path)
        param_names = params_dict['x0'].keys()
        x0 = [params_dict['x0'][key] for key in param_names]
        bounds = [params_dict['bounds'][key] for key in param_names]
        feature_names = params_dict['feature_names']
        objective_names = params_dict['objective_names']
    else:
        param_names = default_x0_dict.keys()
        x0 = [default_x0_dict[key] for key in param_names]
        bounds = [default_bounds_dict[key] for key in param_names]
        feature_names = default_feature_names
        objective_names = default_objective_names

    globals()['adaptive_interval'] = adaptive_interval

    if param_gen is None:
        param_gen = default_param_gen
    if param_gen not in globals():
        raise NameError('Multi-Objective Optimization: %s has not been imported, or is not a valid class of parameter '
                        'generator.' % param_gen)
    else:
        globals()['param_gen'] = param_gen

    if get_features is None:
        get_features = default_get_features
    if get_features not in globals() or not callable(globals()[get_features]):
        raise NameError('Multi-Objective Optimization: get_features: %s has not been imported, or is not a callable '
                        'function.' % get_features)
    globals()['get_features'] = globals()[get_features]

    if get_objectives is None:
        get_objectives = default_get_objectives
    if get_objectives not in globals() or not callable(globals()[get_objectives]):
        raise NameError('Multi-Objective Optimization: get_objectives: %s has not been imported, or is not a callable '
                        'function.' % get_objectives)
    globals()['get_objectives'] = globals()[get_objectives]

    if group_size > num_procs:
        group_size = num_procs
        print 'Multi-Objective Optimization: group_size adjusted to not exceed num_processes: %i' % num_procs
    un_utilized = num_procs % group_size
    iter_per_gen = pop_size / num_procs * group_size
    if iter_per_gen / group_size * num_procs < pop_size:
        iter_per_gen += 1

    print 'Multi-Objective Optimization: %s; Total processes: %i; Population size: %i; Group size: %i; ' \
          'Feature calculator: %s; Objective calculator: %s; Iterations / generation: %i' % \
          (param_gen, num_procs, pop_size, group_size, get_features, get_objectives, iter_per_gen)
    if un_utilized > 0:
        print 'Multi-Objective Optimization: %i processes are unutilized' % un_utilized

    c[:].execute('from parallel_optimize_leak import *', block=True)
    c[:].map_sync(init_engine, [spines] * num_procs, [mech_file_path] * num_procs, [neurotree_file_path] * num_procs,
                  [neurotree_index] * num_procs)


@interactive
def init_engine(spines=False, mech_file_path=None, neurotree_file_path=None, neurotree_index=0):
    """

    :param spines: bool
    :param mech_file_path: str
    :param neurotree_file_path: str
    :param neurotree_index: int
    """
    if mech_file_path is None:
        mech_file_path = default_mech_file_path

    globals()['spines'] = spines

    if neurotree_file_path is None:
        neurotree_file_path = default_neurotree_file_path
    neurotree_dict = read_from_pkl(neurotree_file_path)[neurotree_index]

    global cell
    cell = DG_GC(neurotree_dict=neurotree_dict, mech_file_path=mech_file_path, full_spines=spines)
    # in order to use a single engine to compute many different objectives, we're going to need a way to reset the cell
    # to the original state by re-initializing the mech_dict from the mech_file_path

    # get the thickest apical dendrite ~200 um from the soma
    candidate_branches = []
    candidate_diams = []
    candidate_locs = []
    for branch in cell.apical:
        if ((cell.get_distance_to_node(cell.tree.root, branch, 0.) >= 200.) &
                (cell.get_distance_to_node(cell.tree.root, branch, 1.) > 300.) & (not cell.is_terminal(branch))):
            candidate_branches.append(branch)
            for seg in branch.sec:
                loc = seg.x
                if cell.get_distance_to_node(cell.tree.root, branch, loc) > 250.:
                    candidate_diams.append(branch.sec(loc).diam)
                    candidate_locs.append(loc)
                    break
    index = candidate_diams.index(max(candidate_diams))
    dend = candidate_branches[index]
    dend_loc = candidate_locs[index]

    # get the most distal terminal branch > 300 um from the soma
    candidate_branches = []
    candidate_end_distances = []
    for branch in (branch for branch in cell.apical if cell.is_terminal(branch)):
        if cell.get_distance_to_node(cell.tree.root, branch, 0.) >= 300.:
            candidate_branches.append(branch)
            candidate_end_distances.append(cell.get_distance_to_node(cell.tree.root, branch, 1.))
    index = candidate_end_distances.index(max(candidate_end_distances))
    distal_dend = candidate_branches[index]
    distal_dend_loc = 1.

    global rec_locs
    rec_locs = {'soma': 0., 'dend': dend_loc, 'distal_dend': distal_dend_loc}
    global rec_nodes
    rec_nodes = {'soma': cell.tree.root, 'dend': dend, 'distal_dend': distal_dend}
    global rec_filename
    rec_filename = 'sim_output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'_pid'+str(os.getpid())

    global sim
    sim = QuickSim(duration, verbose=False)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)
    sim.parameters['spines'] = spines


@interactive
def run_optimization(group_size, group_cores, pop_size, pop_cores, generator, test_function, error_function,
                                       x0, bounds):
    global check_bounds
    check_bounds = CheckBounds(bounds[0], bounds[1])
    global generation
    generation = generator.__call__()
    results = run_generation(group_size, group_cores, pop_size, pop_cores, test_function, generation)

    possibles = globals().copy()
    possibles.update(locals())
    method = possibles.get(error_function)
    if not method:
        raise NotImplementedError("Method %s not implemented" % error_function)
    global errors
    errors = method(results)


@interactive
def run_generation(group_size, group_cores, pop_size, pop_cores, test_function, generation):
    c.load_balanced_view()
    start_time = time.time()
    global results
    results = []
    final_results = []
    possibles = globals().copy()
    possibles.update(locals())
    method = possibles.get(test_function)
    if not method:
        raise NotImplementedError("Method %s not implemented" % test_function)

    num_operating_groups = math.floor(pop_cores/group_cores)
    group_indexes = range(0, pop_size)
    for operating_group_ind in range(0, num_operating_groups):
        dv = c[(operating_group_ind*group_cores):(operating_group_ind *group_cores+group_cores-1)]
        group_ind = group_indexes.pop(0)
        formatted_result = method(dv, generation[group_ind].x, group_ind)
        results.append(formatted_result)
    while len(group_indexes) > 0:
        while not np.any(formatted_result[1].ready() for formatted_result in results):
            for formatted_result in results:
                result = formatted_result[1]
                """
                last = []
                time.sleep(1.)
                clear_output()
                for i, stdout in enumerate([stdout for stdout in result.stdout if stdout][-len(sec_list):]):
                    line = stdout.splitlines()[-1]
                    if line not in last:
                        print line
                        last.append(line)
                if len(last) > len(sec_list):
                    last = last[-len(sec_list):]
                sys.stdout.flush()
                """
                last_buffer_len = 0
                time.sleep(1.)
                clear_output()
                stdout = result.stdout
                if stdout:
                    # lines = stdout.splitlines()
                    # above fails because stdout is a list of strings, not a string. Instead, try:
                    lines = []
                    for line_ind in enumerate(stdout):
                        split_lines = stdout[line_ind].splitlines()
                        for line in split_lines:
                            lines.append(line)
                    if len(lines) > last_buffer_len:
                        for line in lines[last_buffer_len:]:
                            print line
                        last_buffer_len = len(lines)
                sys.stdout.flush()
        for operating_group_ind in np.where(formatted_result[1].ready() for formatted_result in results):
            old_result = results[operating_group_ind][1].get()
            old_result_group_index = results[operating_group_ind][0]
            final_results.append([old_result_group_index, result])
            dv = c[(operating_group_ind*group_cores):(operating_group_ind*group_cores+group_cores-1)]
            group_ind = group_indexes.pop(0)
            new_result = method(dv, generation[group_ind].x, group_ind)
            results[operating_group_ind] = new_result
    print('Simulation took %.3f s' % (time.time() - start_time))
    return final_results


@interactive
def test_pas(dv, x, group_ind):
    """
    Distribute simulations across available engines for optimization of leak conductance density gradient.
    :param x: array (soma.g_pas, dend.g_pas slope, dend.g_pas tau, dend.g_pas xhalf)
    :return: float
    """
    if not check_bounds.within_bounds(x, 'pas'):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    sec_list = ['soma', 'dend', 'distal_dend']
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    print 'Process %i using current x: %s: %s' % (os.getpid(), str(param_names), formatted_x)
    result = dv.map_async(get_Rinp_for_section, sec_list, [x for entry in sec_list])
    formatted_result = [group_ind, result]
    return formatted_result

@interactive
def pas_error(results):
    # the target values and acceptable ranges
    target_val = {}
    target_range = {}
    target_val['pas'] = {'soma': 295., 'dend': 375.}
    target_range['pas'] = {'soma': 0.5, 'dend': 1.}
    target_val['v_rest'] = {'soma': v_init, 'tuft_offset': 0.}
    target_range['v_rest'] = {'soma': 0.25, 'tuft_offset': 0.1}

    for result in results:
        group_index = result[0]
        Err = {}
        final_result = {}
        for dict in result[1]:
            final_result.update(dict)
        for section in target_val['pas']:
            Err[section] = ((target_val['pas'][section] - final_result[section]) / target_range['pas'][section]) ** 2.
        section = 'distal_dend'
        # add catch for decreasing terminal end input resistance too much
        if final_result['distal_dend'] < final_result['dend']:
            Err += ((final_result['dend'] - final_result['distal_dend']) / target_range['pas']['dend']) ** 2.
            #fix this!!

        for section in objective_names:
            generation[group_index].objectives[section] = Err[section]

        print 'Process %i: %s: %s; soma R_inp: %.1f, dend R_inp: %.1f, distal_dend R_inp: %.1f; Err: %.3E' % (os.getpid(),
                                                                                    str(xlabels['pas']), formatted_x,
                                                                                    final_result['soma'],
                                                                                    final_result['dend'],
                                                                                    final_result['distal_dend'], Err)
        errors.append(Err)
    return errors


@interactive
def offset_vm(description, vm_target=None):
    """

    :param description: str
    :param vm_target: float
    """
    if vm_target is None:
        vm_target = v_init
    sim.modify_stim(0, amp=0.)
    node = rec_nodes[description]
    loc = rec_locs[description]
    rec_dict = sim.get_rec(description)
    sim.modify_stim(1, node=node, loc=loc, amp=0.)
    rec = rec_dict['vec']
    offset = True
    sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    sim.modify_stim(1, amp=i_holding[description])
    sim.run(vm_target)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < vm_target - 0.5:
        i_holding[description] += 0.01
        while offset:
            if sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(1, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest < vm_target - 0.5:
                i_holding[description] += 0.01
            else:
                offset = False
    elif v_rest > vm_target + 0.5:
        i_holding[description] -= 0.01
        while offset:
            if sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(1, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                i_holding[description] -= 0.01
            else:
                offset = False
    sim.tstop = duration
    return v_rest


@interactive
def update_mech_dict(x):
    update_pas_exp(x)
    cell.export_mech_dict(cell.mech_file_path)


@interactive
def update_pas_exp(x):
    """

    x0 = [2.28e-05, 1.58e-06, 58.4]
    :param x: array [soma.g_pas, dend.g_pas slope, dend.g_pas tau]
    """
    if spines is False:
        cell.reinit_mechanisms(reset_cable=True)
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.modify_mech_param('apical', 'pas', 'g', origin='soma', slope=x[1], tau=x[2])
    for sec_type in ['axon_hill', 'axon', 'ais', 'apical', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')
    if spines is False:
        cell.correct_for_spines()


def print_gpas_cm_values():
    sec_types = ['apical']
    gpas_values = {s: [] for s in sec_types}
    cm_values = {s: [] for s in sec_types}
    for i in [0, 10, 20]:
        node = cell.get_nodes_of_subtype('apical')[i]
        for i, segment in enumerate(node.sec):
            node.sec.push()
            h.pop_section()
            gpas_values['apical'].append(node.sec(segment.x).g_pas)
            cm_values['apical'].append(node.sec(segment.x).cm)
    print 'g_pas: '
    pprint.pprint(gpas_values)
    print 'cm '
    pprint.pprint(cm_values)


@interactive
def get_objectives(*args):
    """

    :param args:
    :return:
    """
    pass


@interactive
def get_Rinp_features(*args):
    """

    :param args:
    :return:
    """
    pass


@interactive
def get_Rinp_for_section(section, x):
    """
    Inject a hyperpolarizing step current into the specified section, and return the steady-state input resistance.
    :param section: str
    :return: dict: {str: float}
    """
    start_time = time.time()
    sim.tstop = duration
    sim.parameters['section'] = section
    sim.parameters['target'] = 'Rinp'
    sim.parameters['optimization'] = 'pas'
    amp = -0.05
    update_pas_exp(x)
    cell.zero_na()
    offset_vm(section)
    loc = rec_locs[section]
    node = rec_nodes[section]
    rec = sim.get_rec(section)
    sim.modify_stim(0, node=node, loc=loc, amp=amp, dur=stim_dur)
    sim.run(v_init)
    Rinp = get_Rinp(np.array(sim.tvec), np.array(rec['vec']), equilibrate, duration, amp)[2]
    result = {section: Rinp}
    print 'Process:', os.getpid(), 'calculated Rinp for %s in %.1f s, Rinp: %.1f' % (section, time.time() - start_time,
                                                                                    Rinp)
    return result


@interactive
def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
        sim.export_to_file(f)

def plot_best(dv, x=None, discard=True):
    """
    Run simulations on the engines with the last best set of parameters, have the engines export their results to .hdf5,
    and then read in and plot the results.
    :param x: array
    """
    if x is None:
        if hist.x_values:
            pas_error(hist.report_best())
        else:
            raise Exception('Please specify input parameters (history might be empty).')
    else:
        pas_error(x)
    dv.execute('export_sim_results()')
    rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir + filename + '.hdf5')]
    for rec_filename in rec_file_list:
        with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
            for trial in f.itervalues():
                target = trial.attrs['target']
                section = trial.attrs['section']
                optimization = trial.attrs['optimization']
                fig, axes = plt.subplots(1)
                for rec in trial['rec'].itervalues():
                    axes.plot(trial['time'], rec, label=rec.attrs['description'])
                axes.legend(loc='best', frameon=False, framealpha=0.5)
                axes.set_xlabel('Time (ms)')
                axes.set_ylabel('Vm (mV)')
                axes.set_title('Optimize %s: %s (%s)' % (optimization, target, section))
                clean_axes(axes)
                fig.tight_layout()
                plt.show()
                plt.close()
    if discard:
        for rec_filename in rec_file_list:
            os.remove(data_dir + rec_filename + '.hdf5')


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_filename) != -1,sys.argv)+1):])