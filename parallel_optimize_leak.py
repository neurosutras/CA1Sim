__author__ = 'Grace Ng'
from function_lib import *
from ipyparallel import interactive
import click
from IPython.display import clear_output
from plot_results import *
from specify_cells3 import *
import scipy.optimize as optimize
import time
import os
import sys
import pprint
import mkl

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


try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass


neurotree_filename = '121516_DGC_trees.pkl'
neurotree_dict = read_from_pkl(morph_dir+neurotree_filename)
history_filename = '030517 leak optimization history'


equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.01
amp = 0.3
th_dvdt = 10.
v_init = -77.
v_active = -77.
i_holding = {'soma': 0., 'dend': 0., 'distal_dend': 0.}
soma_ek = -77.

# the target values and acceptable ranges
target_val = {}
target_range = {}
target_val['pas'] = {'soma': 295., 'dend': 375.}
target_range['pas'] = {'soma': 0.5, 'dend': 1.}
target_val['v_rest'] = {'soma': v_init, 'tuft_offset': 0.}
target_range['v_rest'] = {'soma': 0.25, 'tuft_offset': 0.1}

x0 = {}
xlabels = {}
xmin = {}
xmax = {}

check_bounds = CheckBounds(xmin, xmax)
xlabels['pas'] = ['soma.g_pas', 'dend.g_pas slope', 'dend.g_pas tau']
hist = History()
hist.xlabels = xlabels['pas']

"""
explore_niter = 2100  # max number of iterations to run
polish_niter = 400
take_step = Normalized_Step(x0['pas'], xmin['pas'], xmax['pas'])
minimizer_kwargs = dict(method=null_minimizer)
"""


default_mech_filename = '042717 GC optimizing spike stability'

@click.command()
@click.option("--spines", type=bool, default=False)
@click.option("--mech_filename", type=str, default=None)
@click.option("--cluster-id", type=str, default=None)
@click.option("--group-size", type=int, default=1)
def main(spines, mech_filename, cluster_id, group_size):
    """

    :param spines: bool
    :param mech_filename: str
    :param cluster_id: str
    :param group_size: int
    :return:
    """
    from ipyparallel import Client
    global c
    global all_engines
    global sub_clients

    if cluster_id is not None:
        c = Client(cluster_id=cluster_id)
    else:
        c = Client()

    if mech_filename is None:
        mech_filename = default_mech_filename

    all_engines = c[:]
    all_engines.block = True

    sub_clients = c[::group_size+1]
    sub_clients.block = True
    sub_clients.execute('from parallel_optimize_leak import *', block=True)
    result = sub_clients.map_sync(init_sub_client, [spines for id in sub_clients.targets],
                                  [mech_filename for id in sub_clients.targets], [cluster_id for id in sub_clients.targets],
                                  sub_clients.targets, [group_size for id in sub_clients.targets])
    if spines:
        x0['pas'] = [3.80E-08, 8.08E-07, 6.78E+01]  # Err: 2.527E-09
        xmin['pas'] = [1.0E-18, 1.0E-12, 25.]
        xmax['pas'] = [1.0E-7, 1.0E-4, 400.]
    else:
        x0['pas'] = [1.050E-10, 1.058E-08, 3.886E+01]  # Error: 4.187E-09
        xmin['pas'] = [1.0E-18, 1.0E-12, 25.]
        xmax['pas'] = [1.0E-6, 1.0E-4, 400.]
    #Need to distribute x values to subclients, then evaluate result to determine next sequence
    """
    result = optimize.basinhopping(pas_error, x0['pas'], niter=explore_niter, niter_success=explore_niter,
                                   disp=True, interval=40, minimizer_kwargs=minimizer_kwargs, take_step=take_step)
    """

    # best_x = hist.report_best()
    # hist.export_to_pkl(history_filename)


    # dv['x'] = hist.report_best()
    # dv['x'] = x0['pas']
    # c[0].apply(parallel_optimize_leak_engine.update_mech_dict)
    sys.stdout.flush()
    # plot_best(x0['pas'])
    # print result


@interactive
def init_sub_client(spines=False, mech_filename=None, cluster_id=None, sub_client_id=0, group_size=1):
    """

    :param mech_filename: str
    :param cluster_id: str
    :param sub_client_id: int
    :param group_size: int
    :return:
    """
    from ipyparallel import Client

    global c
    global dv

    if cluster_id is not None:
        c = Client(cluster_id=cluster_id)
    else:
        c = Client()

    if mech_filename is None:
        mech_filename = default_mech_filename

    dv = c[sub_client_id + 1:sub_client_id + 1 + group_size]
    dv.execute('from parallel_optimize_leak import *', block=True)

    dv.map_sync(init_engine, [spines for int in range(group_size)],
                                  [mech_filename for int in range(group_size)])
    global_start_time = time.time()
    # time.sleep(60)
    c.load_balanced_view()

    if 'pas' not in xmin.keys():
        if spines:
            xmin['pas'] = [1.0E-18, 1.0E-12, 25.]
            xmax['pas'] = [1.0E-7, 1.0E-4, 400.]
        else:
            xmin['pas'] = [1.0E-18, 1.0E-12, 25.]
            xmax['pas'] = [1.0E-6, 1.0E-4, 400.]
    return os.getpid(), dv.targets


@interactive
def init_engine(spines=False, mech_filename=None):
    """

    :param spines: bool
    :param mech_filename: str
    :return:
    """
    if mech_filename is None:
        mech_filename = default_mech_filename

    globals()['spines'] = spines
    global sim_description
    if spines:
        sim_description = 'with_spines'
    else:
        sim_description = 'no_spines'

    global cell
    cell = DG_GC(neurotree_dict=neurotree_dict[0], mech_filename=mech_filename, full_spines=spines)
    cell.zero_na()

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
    rec_filename = str(time.strftime('%m%d%Y', time.gmtime())) + '_' + str(time.strftime('%H%M%S', time.gmtime())) + \
                   '_pid' + str(os.getpid()) + '_sim_output'

    global sim
    sim = QuickSim(duration, verbose=False)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
    sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

    for description, node in rec_nodes.iteritems():
        sim.append_rec(cell, node, loc=rec_locs[description], description=description)

@interactive
def pas_error(x):
    """
    Distribute simulations across available engines for optimization of leak conductance density gradient.
    :param x: array (soma.g_pas, dend.g_pas slope, dend.g_pas tau, dend.g_pas xhalf)
    :return: float
    """
    if not check_bounds.within_bounds(x, 'pas'):
        print 'Aborting: Invalid parameter values.'
        return 1e9
    start_time = time.time()

    hist.x_values.append(x)

    sec_list = ['soma', 'dend', 'distal_dend']
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    print 'Process %i using current x: %s: %s' % (os.getpid(), str(xlabels['pas']), formatted_x)
    result = dv.map_async(get_Rinp_for_section, sec_list, [x for entry in sec_list])
    last = []
    while not result.ready():
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
    while not result.ready():
        time.sleep(1.)
        clear_output()
        stdout = result.stdout
        if stdout:
            lines = stdout.splitlines()
            #fails because stdout is a list of strings, not a string
            if len(lines) > last_buffer_len:
                for line in lines[last_buffer_len:]:
                    print line
                last_buffer_len = len(lines)
        sys.stdout.flush()
    """
    result = result.get()

    Err = 0.

    final_result = {}
    for dict in result:
        final_result.update(dict)
    for section in final_result:
        if section not in hist.Rinp_values:
            hist.Rinp_values[section] = []
    for section in target_val['pas']:
        Err += ((target_val['pas'][section] - final_result[section]) / target_range['pas'][section]) ** 2.
        hist.Rinp_values[section].append(final_result[section])
    section = 'distal_dend'
    hist.Rinp_values[section].append(final_result[section])
    # add catch for decreasing terminal end input resistance too much
    if final_result['distal_dend'] < final_result['dend']:
        Err += ((final_result['dend'] - final_result['distal_dend']) / target_range['pas']['dend']) ** 2.
    hist.error_values.append(Err)

    print('Simulation took %.3f s' % (time.time() - start_time))
    print 'Process %i: %s: %s; soma R_inp: %.1f, dend R_inp: %.1f, distal_dend R_inp: %.1f; Err: %.3E' % (os.getpid(),
                                                                                str(xlabels['pas']), formatted_x,
                                                                                final_result['soma'],
                                                                                final_result['dend'],
                                                                                final_result['distal_dend'], Err)
    return Err


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
    cell.export_mech_dict(cell.mech_filename)


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
    sim.parameters['description'] = sim_description
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

def plot_best(x=None, discard=True):
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