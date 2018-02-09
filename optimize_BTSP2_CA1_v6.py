"""
These methods aim to optimize a single parameterization of a model of bidirectional, state-dependent behavioral time 
scale synaptic plasticity to account for the width and amplitude of all place fields in an experimental data set from 
the Magee lab that includes: 
1) Silent cells converted into place cells by spontaneous plateaus
2) Silent cells converted into place cells by experimentally induced plateaus
3) Existing place cells that shift their place field locations after an experimentally induced plateau

Features/assumptions of the phenomenological model:
1) Synaptic weights in a silent cell are all = 1 prior to field induction 1. w(t0) = 1
2) Activity at each synapse generates a long duration 'local signal', or an 'eligibility trace' for plasticity.
3) Dendritic plateaus generate a long duration 'global signal', or an 'availability trace' for plasticity.
4) Synaptic weights w1 at time t1 after a plateau are a function of the initial weights w0 at time t0, and the two
plasticity signals.

Features/assumptions of the mechanistic model:
1) Synaptic strength is equivalent to the number of AMPA-Rs at a synapse (quantal size). 
2) Dendritic plateaus generate a global signal that increases the size of a pool of mobile AMPA-Rs, available for stable
incorporation into synapses.
3) Activity at each synapse generates a local signal that increases the number of eligible slots to capture mobile 
AMPA-Rs.
4) Both signals interact at already potentiated synapses to destabilize captured AMPA-Rs, returning them to the mobile 
pool, and reducing the number of eligible slots. 
5) AMPAR-s can be in 2 states (Markov-style kinetic scheme):

       rMC0 * f1(global_signal * local_signal)
M (mobile) <----------------------> C (captured by a synapse)
       rCM0 * f2(global_signal * local_signal)

6) At rest 100% of non-synaptic receptors are in state M, mobile and available for synaptic capture.
7) The forward transition rate from state M (mobile AMPA-Rs) to state C (captured AMPA-Rs) is proportional to the
product of local signal and global signal, and inversely proportional to the current weight (occupancy of state C).
8) The reverse transition rate from state C to state M depends on a non-monotonic function of the product of local
signal and global signal, and proportional to the current weight (occupancy of state C).
"""
__author__ = 'milsteina'
from BTSP_utils import *
from nested.optimize_utils import *
from scipy.optimize import minimize
import click


context = Context()


def config_interactive(config_file_path=None, output_dir=None, temp_output_path=None, export=False,
                       export_file_path=None, label=None, disp=True, verbose=2, **kwargs):
    """

    :param config_file_path: str (.yaml file path)
    :param output_dir: str (dir path)
    :param temp_output_path: str (.hdf5 file path)
    :param export: bool
    :param export_file_path: str (.hdf5 file path)
    :param label: str
    :param disp: bool
    :param verbose: int
    """

    if config_file_path is not None:
        context.config_file_path = config_file_path
    if 'config_file_path' not in context() or context.config_file_path is None or \
            not os.path.isfile(context.config_file_path):
        raise Exception('config_file_path specifying required parameters is missing or invalid.')
    config_dict = read_from_yaml(context.config_file_path)
    context.param_names = config_dict['param_names']
    if 'default_params' not in config_dict or config_dict['default_params'] is None:
        context.default_params = {}
    else:
        context.default_params = config_dict['default_params']
    for param in context.default_params:
        config_dict['bounds'][param] = (context.default_params[param], context.default_params[param])
    context.bounds = [config_dict['bounds'][key] for key in context.param_names]
    if 'rel_bounds' not in config_dict or config_dict['rel_bounds'] is None:
        context.rel_bounds = None
    else:
        context.rel_bounds = config_dict['rel_bounds']
    if 'x0' not in config_dict or config_dict['x0'] is None:
        context.x0 = None
    else:
        context.x0 = config_dict['x0']
        context.x0_dict = context.x0
        context.x0_array = param_dict_to_array(context.x0_dict, context.param_names)
    context.feature_names = config_dict['feature_names']
    context.objective_names = config_dict['objective_names']
    context.target_val = config_dict['target_val']
    context.target_range = config_dict['target_range']
    context.optimization_title = config_dict['optimization_title']
    context.kwargs = config_dict['kwargs']  # Extra arguments to be passed to imported sources
    context.kwargs['verbose'] = verbose
    context.update(context.kwargs)

    missing_config = []
    if 'update_context' not in config_dict or config_dict['update_context'] is None:
        missing_config.append('update_context')
    else:
        context.update_context_dict = config_dict['update_context']
    if 'get_features_stages' not in config_dict or config_dict['get_features_stages'] is None:
        missing_config.append('get_features_stages')
    else:
        context.stages = config_dict['get_features_stages']
    if 'get_objectives' not in config_dict or config_dict['get_objectives'] is None:
        missing_config.append('get_objectives')
    else:
        context.get_objectives_dict = config_dict['get_objectives']
    if missing_config:
        raise Exception('config_file at path: %s is missing the following required fields: %s' %
                        (context.config_file_path, ', '.join(str(field) for field in missing_config)))

    if label is not None:
        context.label = label
    if 'label' not in context() or context.label is None:
        label = ''
    else:
        label = '_' + context.label

    if output_dir is not None:
        context.output_dir = output_dir
    if 'output_dir' not in context():
        context.output_dir = None
    if context.output_dir is None:
        output_dir_str = ''
    else:
        output_dir_str = context.output_dir + '/'

    if temp_output_path is not None:
        context.temp_output_path = temp_output_path
    if 'temp_output_path' not in context() or context.temp_output_path is None:
        context.temp_output_path = '%s%s_pid%i_%s%s_temp_output.hdf5' % \
                                   (output_dir_str, datetime.datetime.today().strftime('%Y%m%d%H%M'), os.getpid(),
                                    context.optimization_title, label)

    context.export = export
    if export_file_path is not None:
        context.export_file_path = export_file_path
    if 'export_file_path' not in context() or context.export_file_path is None:
        context.export_file_path = '%s%s_%s%s_interactive_exported_output.hdf5' % \
                                   (output_dir_str, datetime.datetime.today().strftime('%Y%m%d%H%M'),
                                    context.optimization_title, label)

    context.update_context_funcs = []
    for source, func_name in context.update_context_dict.iteritems():
        if source == os.path.basename(__file__).split('.')[0]:
            try:
                func = globals()[func_name]
                if not isinstance(func, collections.Callable):
                    raise Exception('update_context function: %s not callable' % func_name)
                context.update_context_funcs.append(func)
            except:
                raise Exception('update_context function: %s not found' % func_name)
    if not context.update_context_funcs:
        raise Exception('update_context function not found')

    context.disp=disp
    config_worker(context.update_context_funcs, context.param_names, context.default_params, context.target_val,
                  context.target_range, context.temp_output_path, context.export_file_path, context.output_dir,
                  context.disp, **context.kwargs)
    update_source_contexts(context.x0_array)


def config_controller(export_file_path, output_dir, **kwargs):
    """

    :param export_file_path: str
    :param output_dir: str
    """
    processed_export_file_path = export_file_path.replace('.hdf5', '_processed.hdf5')
    context.update(locals())
    context.update(kwargs)
    context.data_path = context.output_dir + '/' + context.data_file_name
    init_context()


def config_worker(update_context_funcs, param_names, default_params, target_val, target_range, temp_output_path,
                  export_file_path, output_dir, disp, data_file_name, **kwargs):
    """
    :param update_context_funcs: list of function references
    :param param_names: list of str
    :param default_params: dict
    :param target_val: dict
    :param target_range: dict
    :param temp_output_path: str
    :param export_file_path: str
    :param output_dir: str (dir path)
    :param disp: bool
    :param data_file_name: str (path)
    """
    context.update(kwargs)
    param_indexes = {param_name: i for i, param_name in enumerate(param_names)}
    processed_export_file_path = export_file_path.replace('.hdf5', '_processed.hdf5')
    context.update(locals())
    context.data_path = context.output_dir + '/' + context.data_file_name
    init_context()


def init_context():
    """

    """

    if context.data_path is None or not os.path.isfile(context.data_path):
        raise IOError('init_context: invalid data_path: %s' % context.data_path)
    with h5py.File(context.data_path, 'r') as f:
        dt = f['defaults'].attrs['dt']  # ms
        input_field_width = f['defaults'].attrs['input_field_width']  # cm
        input_field_peak_rate = f['defaults'].attrs['input_field_peak_rate']  # Hz
        num_inputs = f['defaults'].attrs['num_inputs']
        track_length = f['defaults'].attrs['track_length']  # cm
        binned_dx = f['defaults'].attrs['binned_dx']  # cm
        generic_dx = f['defaults'].attrs['generic_dx']  # cm
        default_run_vel = f['defaults'].attrs['default_run_vel']  # cm/s
        generic_position_dt = f['defaults'].attrs['generic_position_dt']  # ms
        default_interp_dx = f['defaults'].attrs['default_interp_dx']  # cm
        ramp_scaling_factor = f['defaults'].attrs['ramp_scaling_factor']
        binned_x = f['defaults']['binned_x'][:]
        generic_x = f['defaults']['generic_x'][:]
        generic_t = f['defaults']['generic_t'][:]
        default_interp_t = f['defaults']['default_interp_t'][:]
        default_interp_x = f['defaults']['default_interp_x'][:]
        extended_x = f['defaults']['extended_x'][:]
        input_rate_maps = f['defaults']['input_rate_maps'][:]
        peak_locs = f['defaults']['peak_locs'][:]
        if 'data_keys' not in context() or context.data_keys is None:
            data_keys = [(int(cell_id), int(induction)) for cell_id in f['data'] for induction in f['data'][cell_id]]
        spont_cell_id_list = [int(cell_id) for cell_id in f['data'] if f['data'][cell_id].attrs['spont']]
    self_consistent_cell_ids = [cell_id for (cell_id, induction) in context.data_keys if induction == 1 and
                                 (cell_id, 2) in context.data_keys]
    down_dt = 10.  # ms, to speed up optimization
    context.update(locals())
    # context.input_EPSPs = compute_EPSP_matrix(input_rate_maps, input_x=binned_x, output_x=generic_x)
    context.sm = StateMachine(dt=down_dt)
    context.cell_id = None
    context.induction = None


def import_data(cell_id, induction):
    """

    :param cell_id: int
    :param induction: int
    """
    cell_id = int(cell_id)
    induction = int(induction)
    if cell_id == context.cell_id and induction == context.induction:
        return
    cell_key = str(cell_id)
    induction_key = str(induction)
    with h5py.File(context.data_path, 'r') as f:
        if cell_key not in f['data'] or induction_key not in f['data'][cell_key]:
            raise KeyError('optimize_BTSP2_CA1: no data found for cell_id: %s, induction: %s' %
                           (cell_key, induction_key))
        else:
            context.cell_id = cell_id
            context.induction = induction
        this_group = f['data'][cell_key][induction_key]
        context.induction_locs = this_group.attrs['induction_locs']
        context.induction_durs = this_group.attrs['induction_durs']
        context.exp_ramp_raw = {'after': this_group['raw']['exp_ramp']['after'][:]}
        if 'before' in this_group['raw']['exp_ramp']:
            context.exp_ramp_raw['before'] = this_group['raw']['exp_ramp']['before'][:]
        context.position = {}
        context.t = {}
        context.current = []
        for category in this_group['processed']['position']:
            context.position[category] = []
            context.t[category] = []
            for i in xrange(len(this_group['processed']['position'][category])):
                lap_key = str(i)
                context.position[category].append(this_group['processed']['position'][category][lap_key][:])
                context.t[category].append(this_group['processed']['t'][category][lap_key][:])
        for i in xrange(len(this_group['processed']['current'])):
            lap_key = str(i)
            context.current.append(this_group['processed']['current'][lap_key][:])
        context.mean_position = this_group['processed']['mean_position'][:]
        context.mean_t = this_group['processed']['mean_t'][:]
        context.exp_ramp = {'after': this_group['processed']['exp_ramp']['after'][:]}
        context.exp_ramp_vs_t = {'after': this_group['processed']['exp_ramp_vs_t']['after'][:]}
        if 'before' in this_group['processed']['exp_ramp']:
            context.exp_ramp['before'] = this_group['processed']['exp_ramp']['before'][:]
            context.exp_ramp_vs_t['before'] = this_group['processed']['exp_ramp_vs_t']['before'][:]
        context.LSA_ramp = {}
        context.LSA_ramp['after'] = this_group['processed']['LSA_ramp']['after'][:]
        if 'before' in this_group['processed']['LSA_ramp']:
            context.LSA_ramp['before'] = this_group['processed']['LSA_ramp']['before'][:]
        context.LSA_weights = {}
        context.LSA_weights['before'] = this_group['processed']['LSA_weights']['before'][:]
        context.LSA_weights['after'] = this_group['processed']['LSA_weights']['after'][:]
        context.complete_run_vel = this_group['complete']['run_vel'][:]
        context.complete_run_vel_gate = this_group['complete']['run_vel_gate'][:]
        context.complete_position = this_group['complete']['position'][:]
        context.complete_t = this_group['complete']['t'][:]
        context.induction_gate = this_group['complete']['induction_gate'][:]
    context.mean_induction_start_loc = np.mean(context.induction_locs)
    context.mean_induction_dur = np.mean(context.induction_durs)
    mean_induction_start_index = np.where(context.mean_position >= context.mean_induction_start_loc)[0][0]
    mean_induction_stop_index = np.where(context.mean_t >= context.mean_t[mean_induction_start_index] +
                                         context.mean_induction_dur)[0][0]
    context.mean_induction_stop_loc = context.mean_position[mean_induction_stop_index]
    induction_start_times = []
    induction_stop_times = []
    track_start_times = []
    track_stop_times = []
    running_position = 0.
    running_t = 0.
    for i, (this_induction_loc, this_induction_dur) in enumerate(zip(context.induction_locs, context.induction_durs)):
        this_induction_start_index = np.where(context.complete_position >= this_induction_loc + running_position)[0][0]
        this_induction_start_time = context.complete_t[this_induction_start_index]
        this_induction_stop_time = this_induction_start_time + this_induction_dur
        track_start_times.append(running_t)
        running_t += len(context.t['induction'][i]) * context.dt
        track_stop_times.append(running_t)
        induction_start_times.append(this_induction_start_time)
        induction_stop_times.append(this_induction_stop_time)
        running_position += context.track_length
    context.induction_start_times = np.array(induction_start_times)
    context.induction_stop_times = np.array(induction_stop_times)
    context.track_start_times = np.array(track_start_times)
    context.track_stop_times = np.array(track_stop_times)
    context.complete_rate_maps = get_complete_rate_maps(context.input_rate_maps, context.binned_x)
    context.down_t = np.arange(context.complete_t[0], context.complete_t[-1] + context.down_dt / 2., context.down_dt)
    context.down_induction_gate = np.interp(context.down_t, context.complete_t, context.induction_gate)
    if context.disp:
        print 'optimize_BTSP2_CA1: process: %i loaded data for cell: %i, induction: %i' % \
              (os.getpid(), cell_id, induction)


def update_source_contexts(x, local_context=None):
    """

    :param x: array
    :param local_context: :class:'Context'
    """
    if local_context is None:
        local_context = context
    for update_func in local_context.update_context_funcs:
        update_func(x, local_context)


def update_model_params(x, local_context):
    """

    :param x: array
    :param local_context: :class:'Context'
    """
    if local_context is None:
        local_context = context
    local_context.update(param_array_to_dict(x, local_context.param_names))


def plot_data():
    """

    """
    fig, axes = plt.subplots(1)
    for group in context.position:
        for i, this_position in enumerate(context.position[group]):
            this_t = context.t[group][i]
            axes.plot(this_t / 1000., this_position, label=group + str(i))
    axes.set_xlabel('Time (s)')
    axes.set_ylabel('Position (cm)')
    axes.set_title('Interpolated position')
    axes.legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
    fig.tight_layout()
    clean_axes(axes)

    fig, axes = plt.subplots(1)
    axes2 = axes.twinx()
    axes.plot(context.complete_t / 1000., context.complete_run_vel)
    axes2.plot(context.complete_t / 1000., context.complete_run_vel_gate, c='k')
    axes.set_xlabel('Time (s)')
    axes.set_ylabel('Running speed (cm/s)')
    clean_axes(axes)
    axes2.tick_params(direction='out')
    axes2.spines['top'].set_visible(False)
    axes2.spines['left'].set_visible(False)
    axes2.get_xaxis().tick_bottom()
    axes2.get_yaxis().tick_right()
    fig.tight_layout()

    fig, axes = plt.subplots(2, 2)
    axes[1][0].plot(context.binned_x, context.exp_ramp['after'])
    if 'before' in context.exp_ramp:
        axes[1][0].plot(context.binned_x, context.exp_ramp['before'])
        axes[1][0].plot(context.binned_x, context.exp_ramp_raw['before'])
    axes[1][0].plot(context.binned_x, context.exp_ramp_raw['after'])
    axes[1][0].set_xlabel('Position (cm)')
    axes[0][0].set_xlabel('Position (cm)')
    axes[1][0].set_ylabel('Ramp amplitude (mV)')
    axes[1][1].set_ylabel('Ramp amplitude (mV)')
    axes[1][1].set_xlabel('Time (s)')
    axes[0][1].set_xlabel('Time (s)')
    axes[0][0].set_ylabel('Induction current (nA)')
    axes[0][1].set_ylabel('Induction gate (a.u.)')
    for i, this_position in enumerate(context.position['induction']):
        this_t = context.t['induction'][i]
        this_current = context.current[i]
        this_induction_gate = np.zeros_like(this_current)
        indexes = np.where(this_current >= 0.5 * np.max(this_current))[0]
        this_induction_gate[indexes] = 1.
        start_index = indexes[0]
        this_induction_loc = context.induction_locs[i]
        this_induction_dur = context.induction_durs[i]
        axes[0][0].plot(this_position, this_current, label='Lap %i: Loc: %i cm, Dur: %i ms' %
                                                           (i, this_induction_loc, this_induction_dur))
        axes[0][1].plot(np.subtract(this_t, this_t[start_index]) / 1000., this_induction_gate)
    mean_induction_index = np.where(context.mean_position >= context.mean_induction_start_loc)[0][0]
    mean_induction_onset = context.mean_t[mean_induction_index]
    peak_val, ramp_width, peak_shift, ratio, start_loc, peak_loc, end_loc = \
        calculate_ramp_features(context.exp_ramp['after'], context.mean_induction_start_loc)
    start_index, peak_index, end_index = get_indexes_from_ramp_bounds_with_wrap(context.binned_x, start_loc, peak_loc,
                                                                                end_loc)
    axes[1][0].scatter(context.binned_x[[start_index, peak_index, end_index]],
                       context.exp_ramp['after'][[start_index, peak_index, end_index]])
    start_index, peak_index, end_index = get_indexes_from_ramp_bounds_with_wrap(context.mean_position, start_loc,
                                                                                peak_loc, end_loc)
    this_shifted_t = np.subtract(context.mean_t, mean_induction_onset) / 1000.
    axes[1][1].plot(this_shifted_t, context.exp_ramp_vs_t['after'])
    axes[1][1].scatter(this_shifted_t[[start_index, peak_index, end_index]],
                       context.exp_ramp_vs_t['after'][[start_index, peak_index, end_index]])
    if 'before' in context.exp_ramp_vs_t:
        axes[1][1].plot(this_shifted_t, context.exp_ramp_vs_t['before'])
    axes[0][0].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
    clean_axes(axes)
    fig.tight_layout()

    fig, axes = plt.subplots(1, 2)
    x_start = context.mean_induction_start_loc
    x_end = context.mean_induction_stop_loc
    max_ramp = max(np.max(context.LSA_ramp['after']), np.max(context.exp_ramp['after']),
                   np.max(context.exp_ramp_raw['after']))
    max_weights = np.max(context.LSA_weights['after'])
    axes[0].plot(context.binned_x, context.LSA_ramp['after'])
    axes[1].plot(context.peak_locs, context.LSA_weights['after'] + 1., label='After induction')
    if 'before' in context.exp_ramp:
        axes[0].plot(context.binned_x, context.LSA_ramp['before'])
        axes[1].plot(context.peak_locs, context.LSA_weights['before'] + 1., label='Before induction')
        max_weights = max(max_weights, np.max(context.LSA_weights['before']))
    max_weights += 1
    axes[0].plot(context.binned_x, context.exp_ramp['after'])
    axes[0].plot(context.binned_x, context.exp_ramp_raw['after'])
    if 'before' in context.exp_ramp:
        axes[0].plot(context.binned_x, context.exp_ramp['before'])
        axes[0].plot(context.binned_x, context.exp_ramp_raw['before'])
        max_ramp = max(max_ramp, np.max(context.LSA_ramp['before']), np.max(context.exp_ramp['before']),
                       np.max(context.exp_ramp_raw['before']))
    axes[0].hlines(max_ramp * 1.1, xmin=x_start, xmax=x_end, linewidth=2, colors='k')
    axes[1].hlines(max_weights * 1.1, xmin=x_start, xmax=x_end, linewidth=2, colors='k')
    axes[0].set_ylim(-1., max_ramp * 1.2)
    axes[1].set_ylim(0.5, max_weights * 1.2)
    axes[0].set_xlabel('Position (cm)')
    axes[1].set_xlabel('Position (cm)')
    axes[0].set_ylabel('Ramp amplitude (mV)')
    axes[1].set_ylabel('Candidate synaptic weights (a.u.)')
    axes[1].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
    fig.suptitle('Cell: %i, Induction: %i' % (context.cell_id, context.induction))
    clean_axes(axes)
    fig.tight_layout()
    plt.show()
    plt.close()


def get_indexes_from_ramp_bounds_with_wrap(x, start, peak, end, min):
    """

    :param x: array
    :param start: float
    :param peak: float
    :param end: float
    :param min: float
    :return: tuple of float: (start_index, peak_index, end_index, min_index)
    """
    peak_index = np.where(x >= peak)[0]
    if np.any(peak_index):
        peak_index = peak_index[0]
    else:
        peak_index = len(x) - 1
    min_index = np.where(x <= min)[0]
    if np.any(min_index):
        min_index = min_index[0]
    else:
        min_index = len(x) - 1
    if start < peak:
        start_index = np.where(x[:peak_index] <= start)[0][-1]
    else:
        start_index = peak_index + np.where(x[peak_index:] <= start)[0]
        if np.any(start_index):
            start_index = start_index[-1]
        else:
            start_index = len(x) - 1
    if end < peak:
        end_index = np.where(x > end)[0][0]
    else:
        end_index = peak_index + np.where(x[peak_index:] > end)[0]
        if np.any(end_index):
            end_index = end_index[0]
        else:
            end_index = len(x) - 1
    return start_index, peak_index, end_index, min_index


def calculate_ramp_features(ramp, induction_loc, offset=False, smooth=False):
    """

    :param ramp: array
    :param induction_loc: float
    :param offset: bool
    :param smooth: bool
    :return tuple of float
    """
    binned_x = context.binned_x
    track_length = context.track_length
    default_interp_x = context.default_interp_x
    extended_binned_x = np.concatenate([binned_x - track_length, binned_x, binned_x + track_length])
    if smooth:
        local_ramp = signal.savgol_filter(ramp, 21, 3, mode='wrap')
    else:
        local_ramp = np.array(ramp)
    extended_binned_ramp = np.concatenate([local_ramp] * 3)
    extended_interp_x = np.concatenate([default_interp_x - track_length, default_interp_x,
                                        default_interp_x + track_length])
    extended_ramp = np.interp(extended_interp_x, extended_binned_x, extended_binned_ramp)
    interp_ramp = extended_ramp[len(default_interp_x):2 * len(default_interp_x)]
    baseline_indexes = np.where(interp_ramp <= np.percentile(interp_ramp, 10.))[0]
    baseline = np.mean(interp_ramp[baseline_indexes])
    if offset:
        interp_ramp -= baseline
        extended_ramp -= baseline
    peak_index = np.where(interp_ramp == np.max(interp_ramp))[0][0] + len(interp_ramp)
    peak_val = extended_ramp[peak_index]
    peak_x = extended_interp_x[peak_index]
    start_index = np.where(extended_ramp[:peak_index] <=
                           0.15 * (peak_val - baseline) + baseline)[0][-1]
    end_index = peak_index + np.where(extended_ramp[peak_index:] <= 0.15 *
                                      (peak_val - baseline) + baseline)[0][0]
    start_loc = float(start_index % len(default_interp_x)) / float(len(default_interp_x)) * track_length
    end_loc = float(end_index % len(default_interp_x)) / float(len(default_interp_x)) * track_length
    peak_loc = float(peak_index % len(default_interp_x)) / float(len(default_interp_x)) * track_length
    min_index = np.where(interp_ramp == np.min(interp_ramp))[0][0] + len(interp_ramp)
    min_val = extended_ramp[min_index]
    min_loc = float(min_index % len(default_interp_x)) / float(len(default_interp_x)) * track_length
    peak_shift = peak_x - induction_loc
    if peak_shift > track_length / 2.:
        peak_shift = -(track_length - peak_shift)
    elif peak_shift < -track_length / 2.:
        peak_shift += track_length
    ramp_width = extended_interp_x[end_index] - extended_interp_x[start_index]
    before_width = induction_loc - start_loc
    if induction_loc < start_loc:
        before_width += track_length
    after_width = end_loc - induction_loc
    if induction_loc > end_loc:
        after_width += track_length
    ratio = before_width / after_width
    return peak_val, ramp_width, peak_shift, ratio, start_loc, peak_loc, end_loc, min_val, min_loc


def wrap_around_and_compress(waveform, interp_x):
    """

    :param waveform: array of len(3 * interp_x)
    :param interp_x: array
    :return: array of len(interp_x)
    """
    before = np.array(waveform[:len(interp_x)])
    after = np.array(waveform[2 * len(interp_x):])
    within = np.array(waveform[len(interp_x):2 * len(interp_x)])
    waveform = within[:len(interp_x)] + before[:len(interp_x)] + after[:len(interp_x)]
    return waveform


def subtract_baseline(waveform, baseline=None):
    """

    :param waveform: array
    :param baseline: float
    :return: array
    """
    new_waveform = np.array(waveform)
    if baseline is None:
        baseline = np.mean(new_waveform[np.where(new_waveform <= np.percentile(new_waveform, 10.))[0]])
    new_waveform -= baseline
    return new_waveform, baseline


def get_complete_rate_maps(input_rate_maps, input_x):
    """
    :param input_rate_maps: array
    :param input_x: array (x resolution of input)
    :return: list of array
    """
    complete_rate_maps = []
    for j in xrange(len(input_rate_maps)):
        this_complete_rate_map = np.array([])
        for group in ['pre', 'induction', 'post']:
            for i, this_position in enumerate(context.position[group]):
                this_rate_map = np.interp(this_position, input_x, input_rate_maps[j])
                this_complete_rate_map = np.append(this_complete_rate_map, this_rate_map)
        this_complete_rate_map = np.multiply(this_complete_rate_map, context.complete_run_vel_gate)
        if len(this_complete_rate_map) != len(context.complete_run_vel_gate):
            print 'get_complete_rate_maps: mismatched array length'
        complete_rate_maps.append(this_complete_rate_map)
    return complete_rate_maps


def get_signal_filters(local_signal_rise, local_signal_decay, global_signal_rise, global_signal_decay, dt=None,
                       plot=False):
    """
    :param local_signal_rise: float
    :param local_signal_decay: float
    :param global_signal_rise: float
    :param global_signal_decay: float
    :param dt: float
    :param plot: bool
    :return: array, array
    """
    max_time_scale = max(local_signal_rise + local_signal_decay, global_signal_rise + global_signal_decay)
    if dt is None:
        dt = context.dt
    filter_t = np.arange(0., 6. * max_time_scale, dt)
    local_filter = np.exp(-filter_t / local_signal_decay) - np.exp(-filter_t / local_signal_rise)
    peak_index = np.where(local_filter == np.max(local_filter))[0][0]
    decay_indexes = np.where(local_filter[peak_index:] < 0.005 * np.max(local_filter))[0]
    if np.any(decay_indexes):
        local_filter = local_filter[:peak_index + decay_indexes[0]]
    local_filter /= np.sum(local_filter)
    local_filter_t = filter_t[:len(local_filter)]
    global_filter = np.exp(-filter_t / global_signal_decay) - np.exp(-filter_t / global_signal_rise)
    peak_index = np.where(global_filter == np.max(global_filter))[0][0]
    decay_indexes = np.where(global_filter[peak_index:] < 0.005 * np.max(global_filter))[0]
    if np.any(decay_indexes):
        global_filter = global_filter[:peak_index + decay_indexes[0]]
    global_filter /= np.sum(global_filter)
    global_filter_t = filter_t[:len(global_filter)]
    if plot:
        fig, axes = plt.subplots(1)
        axes.plot(local_filter_t / 1000., local_filter / np.max(local_filter), color='k',
                  label='Local signal filter')
        axes.plot(global_filter_t / 1000., global_filter / np.max(global_filter), color='r',
                  label='Global signal filter')
        axes.set_xlabel('Time (s)')
        axes.set_ylabel('Normalized filter amplitude')
        axes.set_title('Plasticity signal filters')
        axes.legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
        axes.set_xlim(-0.5, max(5000., local_filter_t[-1], global_filter_t[-1]) / 1000.)
        clean_axes(axes)
        fig.tight_layout()
        plt.show()
        plt.close()
    return local_filter_t, local_filter, global_filter_t, global_filter


def get_local_signal(rate_map, local_filter, dt):
    """

    :param rate_map: array
    :param local_filter: array
    :param dt: float
    :return: array
    """
    return np.convolve(0.001 * dt * rate_map, local_filter)[:len(rate_map)]


def get_global_signal(induction_gate, global_filter):
    """

    :param induction_gate: array
    :param global_filter: array
    :return: array
    """
    return np.convolve(induction_gate, global_filter)[:len(induction_gate)]


def get_local_signal_population(local_filter):
    """

    :param local_filter:
    :return:
    """
    local_signals = []
    for i in xrange(len(context.peak_locs)):
        rate_map = np.interp(context.down_t, context.complete_t, context.complete_rate_maps[i])
        local_signals.append(get_local_signal(rate_map, local_filter, context.down_dt))
    return local_signals


def get_args_static_signal_amplitudes():
    """
    A nested map operation is required to compute model signal amplitudes. The arguments to be mapped are the same
    (static) for each set of parameters.
    :return: list of list
    """
    with h5py.File(context.data_path, 'r') as f:
        data_keys = []
        for cell_key in f['data']:
            for induction_key in f['data'][cell_key]:
                data_keys.append((int(cell_key), int(induction_key)))
    return zip(*data_keys)


def compute_features_signal_amplitudes(x, cell_id=None, induction=None, export=False, plot=False, full_output=False):
    """

    :param x: array
    :param cell_id: int
    :param induction: int
    :param export: bool
    :param plot: bool
    :param full_output: bool
    :return: dict
    """
    import_data(cell_id, induction)
    update_source_contexts(x, context)
    print 'Process: %i: computing signal_amplitude features for cell_id: %i, induction: %i with x: %s' % \
          (os.getpid(), context.cell_id, context.induction, ', '.join('%.3E' % i for i in x))
    start_time = time.time()
    local_filter_t, local_filter, global_filter_t, global_filter = \
        get_signal_filters(context.local_signal_rise, context.local_signal_decay, context.global_signal_rise,
                           context.global_signal_decay, context.down_dt, plot)
    global_signal = get_global_signal(context.down_induction_gate, global_filter)
    local_signals = get_local_signal_population(local_filter)
    local_signal_peaks = [np.max(local_signal) for local_signal in local_signals]
    if plot:
        fig, axes = plt.subplots(1)
        hist, edges = np.histogram(local_signal_peaks, density=True)
        bin_width = edges[1] - edges[0]
        axes.plot(edges[:-1]+bin_width/2., hist * bin_width)
        axes.set_xlabel('Plasticity peak local signal amplitude (a.u.)')
        axes.set_ylabel('Probability')
        axes.set_title('Local signal amplitude distribution')
        clean_axes(axes)
        fig.tight_layout()
        plt.show()
        plt.close()
    result = {'local_signal_peaks': local_signal_peaks,
              'global_signal_peak': np.max(global_signal)
              }
    print 'Process: %i: computing signal_amplitude features for cell_id: %i, induction: %i took %.1f s' % \
          (os.getpid(), context.cell_id, context.induction, time.time() - start_time)
    return {cell_id: {induction: result}}


def filter_features_signal_amplitudes(primitives, features, export=False, plot=False):
    """

    :param primitives: list of dict (each dict contains results from a single simulation)
    :param current_features: dict
    :param export: bool
    :param plot: bool
    :return: dict
    """
    all_inputs_local_signal_peaks = np.array([], dtype='float32')
    each_cell_local_signal_peak = []
    each_cell_global_signal_peak = []
    for this_dict in primitives:
        for cell_id in this_dict:
            for induction_id in this_dict[cell_id]:
                each_cell_global_signal_peak.append(this_dict[cell_id][induction_id]['global_signal_peak'])
                each_cell_local_signal_peak.append(np.max(this_dict[cell_id][induction_id]['local_signal_peaks']))
                all_inputs_local_signal_peaks = np.append(all_inputs_local_signal_peaks,
                                                          this_dict[cell_id][induction_id]['local_signal_peaks'])
    if plot:
        hist, edges = np.histogram(all_inputs_local_signal_peaks, bins=10, density=True)
        bin_width = edges[1] - edges[0]
        fig, axes = plt.subplots(1)
        axes.plot(edges[:-1] + bin_width / 2., hist * bin_width)
        axes.set_xlabel('Plasticity peak local signal amplitude (a.u.)')
        axes.set_ylabel('Probability')
        axes.set_title('Local signal amplitude distribution (all inputs, all cells)')
        clean_axes(axes)
        fig.tight_layout()

        hist, edges = np.histogram(each_cell_local_signal_peak, bins=min(10, len(primitives)), density=True)
        bin_width = edges[1] - edges[0]
        fig, axes = plt.subplots(1)
        axes.plot(edges[:-1] + bin_width / 2., hist * bin_width)
        axes.set_xlabel('Plasticity peak local signal amplitude (a.u.)')
        axes.set_ylabel('Probability')
        axes.set_title('Local signal amplitude distribution (each cell)')
        clean_axes(axes)
        fig.tight_layout()

        hist, edges = np.histogram(each_cell_global_signal_peak, bins=min(10, len(primitives)), density=True)
        bin_width = edges[1] - edges[0]
        fig, axes = plt.subplots(1)
        axes.plot(edges[:-1] + bin_width / 2., hist * bin_width)
        axes.set_xlabel('Plasticity peak global signal amplitude (a.u.)')
        axes.set_ylabel('Probability')
        axes.set_title('Global signal amplitude distribution (each cell)')
        clean_axes(axes)
        fig.tight_layout()
        plt.show()
        plt.close()
    return {'local_signal_max': np.max(each_cell_local_signal_peak),
            'global_signal_max': np.max(each_cell_global_signal_peak)
            }


def get_delta_weights_LSA_error(delta_weights, target_ramp, input_matrix, ramp_x, input_x, bounds=None):
    """

    :param delta_weights: array
    :param target_ramp: array
    :param input_matrix: array
    :param ramp_x: array
    :param input_x: array
    :param bounds: array
    :return: float
    """
    if bounds is not None:
        min_weight, max_weight = bounds
        if np.min(delta_weights) < min_weight or np.max(delta_weights) > max_weight:
            return 1e9
    model_ramp = delta_weights.dot(input_matrix)
    if len(model_ramp) != len(ramp_x):
        model_ramp = np.interp(ramp_x, input_x, model_ramp)
    Err = 0.
    for i in xrange(len(target_ramp)):
        Err += ((target_ramp[i] - model_ramp[i]) / 0.1) ** 2.
    for delta in np.diff(np.insert(delta_weights, 0, delta_weights[-1])):
        # Err += (delta / 0.01) ** 2.
        Err += (delta / 0.0075) ** 2.
    return Err


def get_delta_weights_LSA(target_ramp, input_rate_maps, initial_delta_weights=None, bounds=None, beta=2., ramp_x=None,
                          input_x=None, plot=False, verbose=1):
    """
    Uses least square approximation to estimate a set of weights to match any arbitrary place field ramp, agnostic
    about underlying kernel, induction velocity, etc.
    :param target_ramp: dict of array
    :param input_rate_maps: array; x=default_interp_x
    :param initial_delta_weights: array
    :param bounds: tuple of float
    :param beta: float; regularization parameter
    :param ramp_x: array (spatial resolution of ramp)
    :param input_x: array (spatial resolution of input_rate_maps)
    :param plot: bool
    :param verbose: int
    :return: tuple of array
    """
    if ramp_x is None:
        ramp_x = context.binned_x
    if input_x is None:
        input_x = context.binned_x
    if len(target_ramp) != len(input_x):
        exp_ramp = np.interp(input_x, ramp_x, target_ramp)
    else:
        exp_ramp = np.array(target_ramp)

    ramp_amp, ramp_width, peak_shift, ratio, start_loc, peak_loc, end_loc, min_val, min_loc = {}, {}, {}, {}, {}, {}, \
                                                                                              {}, {}, {}
    ramp_amp['target'], ramp_width['target'], peak_shift['target'], ratio['target'], start_loc['target'], \
        peak_loc['target'], end_loc['target'], min_val['target'], min_loc['target'] = \
        calculate_ramp_features(target_ramp, context.mean_induction_start_loc)

    input_matrix = np.multiply(input_rate_maps, context.ramp_scaling_factor)
    if initial_delta_weights is None:
        [U, s, Vh] = np.linalg.svd(input_matrix)
        V = Vh.T
        D = np.zeros_like(input_matrix)
        D[np.where(np.eye(*D.shape))] = s / (s ** 2. + beta ** 2.)
        input_matrix_inv = V.dot(D.conj().T).dot(U.conj().T)
        initial_delta_weights = exp_ramp.dot(input_matrix_inv)
    initial_ramp = initial_delta_weights.dot(input_matrix)
    if bounds is None:
        bounds = (0., 3.)
    result = minimize(get_delta_weights_LSA_error, initial_delta_weights,
                      args=(target_ramp, input_matrix, ramp_x, input_x, bounds), method='L-BFGS-B',
                      bounds=[bounds] * len(initial_delta_weights), options={'disp': verbose > 1, 'maxiter': 100})
    delta_weights = result.x
    model_ramp = delta_weights.dot(input_matrix)
    if len(model_ramp) != len(ramp_x):
        model_ramp = np.interp(ramp_x, input_x, model_ramp)
    ramp_amp['model'], ramp_width['model'], peak_shift['model'], ratio['model'], start_loc['model'], \
        peak_loc['model'], end_loc['model'], min_val['model'], min_loc['model'] = \
        calculate_ramp_features(model_ramp, context.mean_induction_start_loc)

    if verbose > 1:
        print 'Process: %i: get_delta_weights_LSA:' % os.getpid()
        print 'exp: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, peak_loc: %.1f, ' \
              'end_loc: %.1f' % (ramp_amp['target'], ramp_width['target'], peak_shift['target'], ratio['target'],
                                 start_loc['target'], peak_loc['target'], end_loc['target'])
        print 'model: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, peak_loc: %.1f' \
              ', end_loc: %.1f' % (ramp_amp['model'], ramp_width['model'], peak_shift['model'], ratio['model'],
                                   start_loc['model'], peak_loc['model'], end_loc['model'])
    sys.stdout.flush()

    if plot:
        x_start = context.mean_induction_start_loc
        x_end = context.mean_induction_stop_loc
        ylim = max(np.max(target_ramp), np.max(model_ramp))
        ymin = min(np.min(target_ramp), np.min(model_ramp))
        fig, axes = plt.subplots(1)
        axes.plot(ramp_x, target_ramp, label='Experiment', color='k')
        axes.plot(ramp_x, initial_ramp, label='Model (Initial)', color='c')
        axes.plot(ramp_x, model_ramp, label='Model (LSA)', color='c')
        axes.hlines(ylim + 0.2, xmin=x_start, xmax=x_end, linewidth=2, colors='k')
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Ramp amplitude (mV)')
        axes.set_xlim([0., context.track_length])
        axes.set_ylim([math.floor(ymin), max(math.ceil(ylim), ylim + 0.4)])
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        axes.set_title('Vm ramp')
        clean_axes(axes)
        fig.tight_layout()
        ylim = np.max(delta_weights) + 1.
        ymin = np.min(delta_weights) + 1.
        fig1, axes1 = plt.subplots(1)
        axes1.plot(context.peak_locs, initial_delta_weights + 1., c='r', label='Model (Initial)')
        axes1.plot(context.peak_locs, delta_weights + 1., c='c', label='Model (LSA)')
        axes1.hlines(ylim + 0.2, xmin=x_start, xmax=x_end, linewidth=2, colors='k')
        axes1.set_xlabel('Location (cm)')
        axes1.set_ylabel('Candidate synaptic weights (a.u.)')
        axes1.set_xlim([0., context.track_length])
        axes1.set_ylim([math.floor(ymin), max(math.ceil(ylim), ylim + 0.4)])
        axes1.legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes1)
        fig1.tight_layout()
        plt.show()
        plt.close()

    return model_ramp, delta_weights


def calculate_model_ramp(initial_weights=None, local_signal_peak=None, global_signal_peak=None, export=False,
                         plot=False):
    """

    :param initial_weights: array
    :param local_signal_peak: float
    :param global_signal_peak: float
    :param export: bool
    :param plot: bool
    :return: dict
    """
    local_filter_t, local_filter, global_filter_t, global_filter = \
        get_signal_filters(context.local_signal_rise, context.local_signal_decay, context.global_signal_rise,
                           context.global_signal_decay, context.down_dt, plot)
    global_signal = get_global_signal(context.down_induction_gate, global_filter) / global_signal_peak
    local_signals = np.divide(get_local_signal_population(local_filter), local_signal_peak)
    signal_xrange = np.linspace(0., 1., 10000)
    pot_rate = sigmoid_segment(context.rMC_slope, context.rMC_th)
    raw_depot_rate = lambda signal: (np.exp(-signal / context.rCM_decay) -
                                     np.exp(-signal / context.rCM_rise)) ** context.rCM_k
    depot_rate_peak_loc = np.log(context.rCM_decay / context.rCM_rise) * context.rCM_decay * context.rCM_rise / \
                        (context.rCM_decay - context.rCM_rise)
    depot_rate_norm_factor = raw_depot_rate(depot_rate_peak_loc)
    depot_rate = lambda signal: (1. / depot_rate_norm_factor) * (np.exp(-signal / context.rCM_decay) -
                                                                 np.exp(-signal / context.rCM_rise)) ** context.rCM_k
    if plot:
        fig, axes = plt.subplots(1)
        axes.plot(signal_xrange, pot_rate(signal_xrange), label='Potentiation rate')
        axes.plot(signal_xrange, depot_rate(signal_xrange), label='Depotentiation rate')
        axes.set_xlabel('Normalized plasticity signal amplitude (a.u.)')
        axes.set_ylabel('Normalized rate')
        axes.set_title('Plasticity signal transformations')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes)
        fig.tight_layout()
        #plt.show()
        #plt.close()
    weights = []
    peak_weight = context.peak_delta_weight + 1.
    if initial_weights is None:
        description = 'model_ramp_features'
        initial_delta_weights = context.LSA_weights['before']
        # re-compute initial weights if they are out of the current weight bounds
        if not np.all((context.min_delta_weight < initial_delta_weights) &
                      (initial_delta_weights < context.peak_delta_weight)):
            initial_ramp, initial_delta_weights = get_delta_weights_LSA(context.exp_ramp['before'],
                                                                        context.input_rate_maps,
                                                                        initial_delta_weights = initial_delta_weights,
                                                                        bounds=(context.min_delta_weight,
                                                                                context.peak_delta_weight),
                                                                        verbose=context.verbose)
        initial_weights = np.divide(np.add(initial_delta_weights, 1.), peak_weight)
    else:
        description = 'model_ramp_self_consistent_features'
        initial_weights = np.divide(initial_weights, peak_weight)
    for i in xrange(len(context.peak_locs)):
        # normalize total number of receptors
        initial_weight = initial_weights[i]
        available = 1. - initial_weight
        context.sm.update_states({'M': available, 'C': initial_weight})
        local_signal = local_signals[i]
        dual_signal_product = np.multiply(global_signal, local_signal)
        context.sm.update_rates(
            {'M': {'C': context.rMC0 * pot_rate(dual_signal_product)},
             'C': {'M': context.rCM0 * depot_rate(dual_signal_product)}})
        context.sm.reset()
        context.sm.run()
        if i == 100:
            example_weight_dynamics = np.array(context.sm.states_history['C'][:-1]) * peak_weight
            example_local_signal = np.array(local_signal)
            example_dual_signal_product = np.array(dual_signal_product)
            if plot:
                fig, axes = plt.subplots(3, sharex=True)
                ymax0 = max(np.max(local_signal), np.max(global_signal))
                bar_loc0 = ymax0 * 1.05
                ymax1 = np.max(dual_signal_product)
                bar_loc1 = ymax1 * 1.05
                axes[0].plot(context.down_t / 1000., local_signal, c='k', label='Local signal')
                axes[0].plot(context.down_t / 1000., global_signal, c='r', label='Global signal')
                # axes[0].set_xlim([-1., context.track_stop_times[0] / 1000. + 1.])
                axes[0].set_ylim([-0.1 * ymax0, 1.1 * ymax0])
                axes[0].hlines([bar_loc0] * len(context.induction_start_times),
                               xmin=context.induction_start_times / 1000.,
                               xmax=context.induction_stop_times / 1000., linewidth=2)
                axes[0].set_xlabel('Time (s)')
                axes[0].set_ylabel('Dual plasticity\nsignal amplitudes')
                axes[0].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
                axes[1].plot(context.down_t / 1000., example_dual_signal_product)
                axes[1].set_ylim([-0.1 * ymax1, 1.1 * ymax1])
                axes[1].hlines([bar_loc1] * len(context.induction_start_times),
                               xmin=context.induction_start_times / 1000.,
                               xmax=context.induction_stop_times / 1000., linewidth=2)
                axes[1].set_xlabel('Time (s)')
                axes[1].set_ylabel('Plasticity\nsignal product')
                axes[2].plot(context.down_t / 1000., example_weight_dynamics)
                axes[2].set_ylim([0., peak_weight * 1.1])
                axes[2].hlines([peak_weight * 1.05] * len(context.induction_start_times),
                               xmin=context.induction_start_times / 1000.,
                               xmax=context.induction_stop_times / 1000., linewidth=2)
                axes[2].set_ylabel('Synaptic weight\n(example\nsingle input)')
                axes[2].set_xlabel('Time (s)')
                clean_axes(axes)
                fig.tight_layout(h_pad=2.)
                #plt.show()
                #plt.close()
        weights.append(context.sm.states['C'] * peak_weight)
    initial_weights = np.multiply(initial_weights, peak_weight)
    weights = np.array(weights)
    delta_weights = np.subtract(weights, initial_weights)
    ramp_amp, ramp_width, peak_shift, ratio, start_loc, peak_loc, end_loc, min_val, min_loc = {}, {}, {}, {}, {}, {}, \
                                                                                              {}, {}, {}
    target_ramp = context.exp_ramp['after']
    ramp_amp['target'], ramp_width['target'], peak_shift['target'], ratio['target'], start_loc['target'], \
    peak_loc['target'], end_loc['target'], min_val['target'], min_loc['target'] = \
        calculate_ramp_features(target_ramp, context.mean_induction_start_loc)
    model_ramp = get_model_ramp(np.subtract(weights, 1.))
    ramp_amp['model'], ramp_width['model'], peak_shift['model'], ratio['model'], start_loc['model'], \
    peak_loc['model'], end_loc['model'], min_val['model'], min_loc['model'] = \
        calculate_ramp_features(model_ramp, context.mean_induction_start_loc)
    if context.disp:
        print 'Process: %i; cell: %i; induction: %i:' % (os.getpid(), context.cell_id, context.induction)
        print 'exp: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, peak_loc: %.1f, ' \
              'end_loc: %.1f, min_val: %.1f, min_loc: %.1f' % \
              (ramp_amp['target'], ramp_width['target'], peak_shift['target'], ratio['target'], start_loc['target'],
               peak_loc['target'], end_loc['target'], min_val['target'], min_loc['target'])
        print 'model: amp: %.1f, ramp_width: %.1f, peak_shift: %.1f, asymmetry: %.1f, start_loc: %.1f, peak_loc: %.1f' \
              ', end_loc: %.1f, min_val: %.1f, min_loc: %.1f' % \
              (ramp_amp['model'], ramp_width['model'], peak_shift['model'], ratio['model'], start_loc['model'],
               peak_loc['model'], end_loc['model'], min_val['model'], min_loc['model'])
        sys.stdout.flush()

    result = {'delta_amp': ramp_amp['model'] - ramp_amp['target'],
              'delta_width': ramp_width['model'] - ramp_width['target'],
              'delta_peak_shift': peak_shift['model'] - peak_shift['target'],
              'delta_asymmetry': ratio['model'] - ratio['target'],
              'delta_min_val': min_val['model'] - min_val['target']}
    abs_delta_min_loc = abs(min_loc['model'] - min_loc['target'])
    if min_loc['model'] <= min_loc['target']:
        if abs_delta_min_loc > context.track_length / 2.:
            delta_min_loc = context.track_length - abs_delta_min_loc
        else:
            delta_min_loc = -abs_delta_min_loc
    else:
        if abs_delta_min_loc > context.track_length / 2.:
            delta_min_loc = -(context.track_length - abs_delta_min_loc)
        else:
            delta_min_loc = abs_delta_min_loc
    result['delta_min_loc'] = delta_min_loc
    if plot:
        bar_loc = max(10., np.max(model_ramp) + 1., np.max(target_ramp) + 1.) * 0.95
        fig, axes = plt.subplots(2)
        axes[1].plot(context.peak_locs, delta_weights)
        axes[1].hlines(peak_weight * 1.05, xmin=context.mean_induction_start_loc, xmax=context.mean_induction_stop_loc)
        axes[0].plot(context.binned_x, target_ramp, label='Experiment')
        axes[0].plot(context.binned_x, model_ramp, label='Model')
        axes[0].hlines(bar_loc, xmin=context.mean_induction_start_loc, xmax=context.mean_induction_stop_loc)
        axes[1].set_ylabel('Change in\nsynaptic weight')
        axes[1].set_xlabel('Location (cm)')
        axes[0].set_ylabel('Subthreshold\ndepolarization (mV)')
        axes[0].set_xlabel('Location (cm)')
        axes[0].legend(loc='best', frameon=False, framealpha=0.5, handlelength=1)
        axes[0].set_ylim([min(-1., np.min(model_ramp) - 1., np.min(target_ramp) - 1.),
                          max(10., np.max(model_ramp) + 1., np.max(target_ramp) + 1.)])
        axes[1].set_ylim([-peak_weight, peak_weight * 1.1])
        clean_axes(axes)
        fig.suptitle('Cell_id: %i, Induction: %i' % (context.cell_id, context.induction))
        fig.tight_layout()
        plt.show()
        plt.close()
    result['model_weights'] = weights
    result['initial_weights'] = initial_weights
    result['model_ramp'] = np.array(model_ramp)
    result['target_ramp'] = np.array(target_ramp)
    if export:
        with h5py.File(context.processed_export_file_path, 'a') as f:
            this_description = 'model_ramp_context'
            if this_description not in f:
                f.create_group(this_description)
                group = f[this_description]
                group.create_dataset('peak_locs', compression='gzip', compression_opts=9, data=context.peak_locs)
                group.create_dataset('binned_x', compression='gzip', compression_opts=9, data=context.binned_x)
                group.create_dataset('signal_xrange', compression='gzip', compression_opts=9,
                                     data=signal_xrange)
                group.create_dataset('pot_rate', compression='gzip', compression_opts=9,
                                     data=pot_rate(signal_xrange))
                group.create_dataset('depot_rate', compression='gzip', compression_opts=9,
                                     data=depot_rate(signal_xrange))
                group.attrs['peak_weight'] = peak_weight
            if description not in f:
                f.create_group(description)
            cell_key = str(context.cell_id)
            induction_key = str(context.induction)
            if cell_key not in f[description]:
                f[description].create_group(cell_key)
            if induction_key not in f[description][cell_key]:
                f[description][cell_key].create_group(induction_key)
            group = f[description][cell_key][induction_key]
            group.create_dataset('target_ramp', compression='gzip', compression_opts=9, data=target_ramp)
            group.create_dataset('model_ramp', compression='gzip', compression_opts=9, data=model_ramp)
            group.create_dataset('model_weights', compression='gzip', compression_opts=9, data=weights)
            group.create_dataset('initial_weights', compression='gzip', compression_opts=9, data=initial_weights)
            group.create_dataset('example_local_signal', compression='gzip', compression_opts=9,
                                 data=example_local_signal)
            group.create_dataset('example_dual_signal_product', compression='gzip', compression_opts=9,
                                 data=example_dual_signal_product)
            group.create_dataset('global_signal', compression='gzip', compression_opts=9, data=global_signal)
            group.create_dataset('down_t', compression='gzip', compression_opts=9, data=context.down_t)
            group.create_dataset('example_weight_dynamics', compression='gzip', compression_opts=9,
                                 data=example_weight_dynamics)
            group.attrs['mean_induction_start_loc'] = context.mean_induction_start_loc
            group.attrs['mean_induction_stop_loc'] = context.mean_induction_stop_loc
            group.attrs['induction_start_times'] = context.induction_start_times
            group.attrs['induction_stop_times'] = context.induction_stop_times
            group.attrs['track_start_times'] = context.track_start_times
            group.attrs['track_stop_times'] = context.track_stop_times
    return {context.cell_id: {context.induction: result}}


def get_args_dynamic_model_ramp(x, features):
    """
    A nested map operation is required to compute model_ramp features. The arguments to be mapped depend on each set of
    parameters and prior features (dynamic).
    :param x: array
    :param features: dict
    :return: list of list
    """
    group_size = len(context.data_keys)
    return [list(item) for item in zip(*context.data_keys)] + [[features['local_signal_max']] * group_size] + \
           [[features['global_signal_max']] * group_size]


def compute_features_model_ramp(x, cell_id=None, induction=None, local_signal_peak=None, global_signal_peak=None,
                                export=False, plot=False):
    """

    :param x: array
    :param cell_id: int
    :param induction: int
    :param local_signal_peak: float
    :param global_signal_peak: float
    :param export: bool
    :param plot: bool
    :return: dict
    """
    import_data(cell_id, induction)
    update_source_contexts(x, context)
    start_time = time.time()
    print 'Process: %i: computing model_ramp_features for cell_id: %i, induction: %i with x: %s' % \
          (os.getpid(), context.cell_id, context.induction, ', '.join('%.3E' % i for i in x))
    result = calculate_model_ramp(local_signal_peak=local_signal_peak, global_signal_peak=global_signal_peak,
                                  export=export, plot=plot)
    print 'Process: %i: computing model_ramp_features for cell_id: %i, induction: %i took %.1f s' % \
          (os.getpid(), context.cell_id, context.induction, time.time() - start_time)
    return result


def filter_features_model_ramp(primitives, current_features, export=False):
    """

    :param primitives: list of dict (each dict contains results from a single simulation)
    :param current_features: dict
    :param export: bool
    :return: dict
    """
    features = {}
    groups = ['spont', 'exp1', 'exp2']
    feature_names = ['delta_amp', 'delta_width', 'delta_peak_shift', 'delta_asymmetry', 'delta_min_loc',
                      'delta_min_val']
    features['raw'] = {}
    for this_result_dict in primitives:
        for cell_id in this_result_dict:
            cell_id = int(cell_id)
            if cell_id not in features['raw']:
                features['raw'][cell_id] = {}
            for induction in this_result_dict[cell_id]:
                induction = int(induction)
                if induction not in features['raw'][cell_id]:
                    features['raw'][cell_id][induction] = {}
                if cell_id in context.spont_cell_id_list:
                    group = 'spont'
                else:
                    group = 'exp' + str(induction)
                for feature_name in feature_names:
                    key = group + '_' + feature_name
                    if key not in features:
                        features[key] = []
                    features[key].append(this_result_dict[cell_id][induction][feature_name])
                    features['raw'][cell_id][induction][feature_name] = \
                        this_result_dict[cell_id][induction][feature_name]
                residuals = np.subtract(this_result_dict[cell_id][induction]['model_ramp'],
                                        this_result_dict[cell_id][induction]['target_ramp'])
                features['raw'][cell_id][induction]['residuals'] = residuals
                features['raw'][cell_id][induction]['model_weights'] = \
                    this_result_dict[cell_id][induction]['model_weights']
                delta_weights = np.subtract(this_result_dict[cell_id][induction]['model_weights'],
                                        this_result_dict[cell_id][induction]['initial_weights'])
                features['raw'][cell_id][induction]['delta_weights'] = delta_weights
    for feature_name in feature_names:
        for group in groups:
            key = group + '_' + feature_name
            if key in features and len(features[key]) > 0:
                features[key] = np.mean(features[key])
            else:
                features[key] = 0.
    return features


def get_args_dynamic_self_consistent_model_ramp(x, features):
    """
    A nested map operation is required to compute model_ramp features. The arguments to be mapped depend on each set of
    parameters and prior features (dynamic).
    :param x: array
    :param features: dict
    :return: list of list
    """
    self_consistent_initial_weights = []
    self_consistent_data_keys = []
    for cell_id in context.self_consistent_cell_ids:
        try:
            this_model_weights = features['raw'][cell_id][1]['model_weights']
        except Exception:
            raise KeyError('get_args_dynamic_self_consistent_model_ramp: model_weights for induction 1 not found in '
                           'features')
        self_consistent_initial_weights.append(this_model_weights)
        self_consistent_data_keys.append((cell_id, 2))
    group_size = len(self_consistent_data_keys)
    return [list(item) for item in zip(*self_consistent_data_keys)] + \
           [[features['local_signal_max']] * group_size] + \
           [[features['global_signal_max']] * group_size] + [self_consistent_initial_weights]


def compute_features_self_consistent_model_ramp(x, cell_id=None, induction=None, local_signal_peak=None,
                                                global_signal_peak=None, self_consistent_initial_weights=None,
                                                export=False, plot=False):
    """

    :param x: array
    :param cell_id: int
    :param induction: int
    :param local_signal_peak: float
    :param global_signal_peak: float
    :param self_consistent_initial_weights: array
    :param export: bool
    :param plot: bool
    :return: dict
    """
    import_data(cell_id, induction)
    update_source_contexts(x, context)
    start_time = time.time()
    print 'Process: %i: computing self_consistent_model_ramp_features for cell_id: %i, induction: %i with x: %s' % \
          (os.getpid(), context.cell_id, context.induction, ', '.join('%.3E' % i for i in x))
    result = calculate_model_ramp(initial_weights=self_consistent_initial_weights, local_signal_peak=local_signal_peak,
                                  global_signal_peak=global_signal_peak, export=export, plot=plot)
    print 'Process: %i: computing self_consistent_model_ramp_features for cell_id: %i, induction: %i took %.1f s' % \
          (os.getpid(), context.cell_id, context.induction, time.time() - start_time)
    return result


def filter_features_self_consistent_model_ramp(primitives, current_features, export=False):
    """

    :param primitives: list of dict (each dict contains results from a single simulation)
    :param current_features: dict
    :param export: bool
    :return: dict
    """
    features = {}
    features['self_consistent'] = {}
    for this_result_dict in primitives:
        for cell_id in this_result_dict:
            for induction in this_result_dict[cell_id]:
                cell_id = int(cell_id)
                induction = int(induction)
                if cell_id not in context.self_consistent_cell_ids:
                    raise KeyError('filter_features_self_consistent_model_ramp: received invalid cell_id: %i' % cell_id)
                if cell_id not in features['self_consistent']:
                    features['self_consistent'][cell_id] = {}
                if induction not in features['self_consistent'][cell_id]:
                    features['self_consistent'][cell_id][induction] = {}
                try:
                    residuals = np.subtract(this_result_dict[cell_id][induction]['model_ramp'],
                                            this_result_dict[cell_id][induction]['target_ramp'])
                except Exception:
                    raise KeyError('filter_features_self_consistent_model_ramp: could not compute residuals for '
                                   'cell_id: %i' % cell_id)
                features['self_consistent'][cell_id][induction]['residuals'] = residuals
    return features


def get_objectives(features):
    """

    :param features: dict
    :return: tuple of dict
    """
    objectives = {}
    feature_names = ['delta_amp', 'delta_width', 'delta_peak_shift', 'delta_asymmetry', 'delta_min_loc',
                     'delta_min_val']
    for cell_id in features['raw']:
        cell_id = int(cell_id)
        for induction in features['raw'][cell_id]:
            induction = int(induction)
            if cell_id in context.spont_cell_id_list:
                group = 'spont'
            else:
                group = 'exp' + str(induction)
            for feature_name in feature_names:
                objective_name = group + '_' + feature_name
                if objective_name not in objectives:
                    objectives[objective_name] = []
                objectives[objective_name].append((features['raw'][cell_id][induction][feature_name] /
                                                  context.target_range[objective_name]) ** 2.)
            feature_name = 'residuals'
            objective_name = group + '_' + feature_name
            if objective_name not in objectives:
                objectives[objective_name] = []
            this_residual_sum = np.sum([(val / context.target_range[objective_name]) ** 2. for
                                        val in features['raw'][cell_id][induction][feature_name]])
            objectives[objective_name].append(this_residual_sum)
            feature_name = 'delta_weights'
            objective_name = 'smoothness'
            if objective_name not in objectives:
                objectives[objective_name] = []
            this_weights = features['raw'][cell_id][induction][feature_name]
            this_delta = np.diff(np.insert(this_weights, 0, this_weights[-1]))
            this_smoothness = np.sum([(val / context.target_range[objective_name]) ** 2. for val in this_delta])
            objectives[objective_name].append(this_smoothness)
    for cell_id in features['self_consistent']:
        cell_id = int(cell_id)
        for induction in features['self_consistent'][cell_id]:
            induction = int(induction)
            feature_name = 'residuals'
            objective_name = 'self_consistent' + '_' + feature_name
            if objective_name not in objectives:
                objectives[objective_name] = []
            this_residual_sum = np.sum([(val / context.target_range[objective_name]) ** 2. for
                                        val in features['self_consistent'][cell_id][induction][feature_name]])
            objectives[objective_name].append(this_residual_sum)
    objective_names = context.target_val.keys()
    for objective_name in objective_names:
        if objective_name not in objectives:
            objectives[objective_name] = 0.
        else:
            objectives[objective_name] = np.mean(objectives[objective_name])
    return features, objectives


def get_model_ramp(delta_weights, input_x=None, ramp_x=None, plot=False):
    """

    :param delta_weights: array
    :param input_x: array (x resolution of inputs)
    :param ramp_x: array (x resolution of ramp)
    :param plot: bool
    :return: array
    """
    if input_x is None:
        input_x = context.binned_x
    if ramp_x is None:
        ramp_x = context.binned_x
    model_ramp = np.multiply(delta_weights.dot(context.input_rate_maps), context.ramp_scaling_factor)
    if len(model_ramp) != len(ramp_x):
        model_ramp = np.interp(ramp_x, input_x, model_ramp)
    if plot:
        fig, axes = plt.subplots(1)
        max_ramp = max(np.max(model_ramp), np.max(context.exp_ramp['after']))
        axes.hlines(max_ramp * 1.2,
                       xmin=context.mean_induction_start_loc,
                       xmax=context.mean_induction_stop_loc, linewidth=2)
        axes.plot(context.binned_x, context.exp_ramp['after'], label='Experiment')
        axes.plot(context.binned_x, model_ramp, label='Model')
        axes.set_ylim(-0.5, max_ramp * 1.4)
        axes.set_xlabel('Location (cm)')
        axes.set_ylabel('Ramp amplitude (mV)')
        axes.set_title('Cell_id: %i, Induction: %i' % (context.cell_id, context.induction))
        clean_axes(axes)
        fig.tight_layout()
        plt.show()
        plt.close()
    return model_ramp


def get_features_interactive(x, plot=False):
    """

    :param x:
    :param plot:
    :return: dict
    """
    features = {}
    args = get_args_static_signal_amplitudes()
    group_size = len(args[0])
    sequences = [[x] * group_size] + args + [[context.export] * group_size]  # + [[plot] * group_size]
    primitives = map(compute_features_signal_amplitudes, *sequences)
    new_features = filter_features_signal_amplitudes(primitives, features, context.export, plot)
    features.update(new_features)

    args = get_args_dynamic_model_ramp(x, features)
    group_size = len(args[0])
    sequences = [[x] * group_size] + args + [[context.export] * group_size] + [[plot] * group_size]
    primitives = map(compute_features_model_ramp, *sequences)
    new_features = filter_features_model_ramp(primitives, features, context.export)
    features.update(new_features)

    args = get_args_dynamic_self_consistent_model_ramp(x, features)
    group_size = len(args[0])
    sequences = [[x] * group_size] + args + [[context.export] * group_size] + [[plot] * group_size]
    primitives = map(compute_features_self_consistent_model_ramp, *sequences)
    new_features = filter_features_self_consistent_model_ramp(primitives, features, context.export)
    features.update(new_features)

    return features


def get_model_ramp_error(x, check_bounds=None, plot=False, full_output=False):
    """

    :param x:
    :param check_bounds: callable
    :param plot: bool
    :param full_output: bool
    :return: float
    """
    if check_bounds is not None:
        if not check_bounds(x):
            return 1e9
    features = get_features_interactive(x, plot=plot)
    features, objectives = get_objectives(features)
    Err = np.sum(objectives.values())
    print 'Process: %i: absolute sum of model_ramp objectives: %.4E' % (os.getpid(), Err)
    if full_output:
        return Err, features, objectives
    else:
        return Err


@click.command()
@click.option("--config-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False),
              default='config/optimize_BTSP2_CA1_v6_cell1_config.yaml')
@click.option("--output-dir", type=click.Path(exists=True, file_okay=False, dir_okay=True), default='data')
@click.option("--export", is_flag=True)
@click.option("--export-file-path", type=str, default=None)
@click.option("--label", type=str, default=None)
@click.option("--disp", is_flag=True)
@click.option("--verbose", type=int, default=1)
@click.option("--serial-optimize", is_flag=True)
@click.option("--plot", is_flag=True)
@click.option("--debug", is_flag=True)
def main(config_file_path, output_dir, export, export_file_path, label, disp, verbose, serial_optimize, plot, debug):
    """

    :param config_file_path: str (path)
    :param output_dir: str (path)
    :param export: bool
    :param export_file_path: str
    :param label: str
    :param disp: bool
    :param verbose: int
    :param serial_optimize: bool
    :param plot: bool
    :param debug: bool
    """
    # requires a global variable context: :class:'Context'

    context.update(locals())
    config_interactive(config_file_path=config_file_path, output_dir=output_dir, export=export,
                       export_file_path=export_file_path, label=label, disp=disp, verbose=verbose)
    rel_bounds_handler = RelativeBoundedStep(context.x0_array, context.param_names, context.bounds, context.rel_bounds)
    x1_array = context.x0_array

    if debug:
        features = get_features_interactive(x1_array, plot=plot)
        features, objectives = get_objectives(features)
        print 'features:'
        pprint.pprint({key: val for (key, val) in features.iteritems() if key in context.feature_names})
        print 'objectives'
        pprint.pprint({key: val for (key, val) in objectives.iteritems() if key in context.objective_names})

    if serial_optimize:
        result = minimize(get_model_ramp_error, x1_array, method='Nelder-Mead',
                                   args=(rel_bounds_handler.check_bounds,),
                                   options={'disp': verbose > 1, 'maxiter': 20})
        x1_array = result.x
    context.update(locals())


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(os.path.basename(__file__)) != -1, sys.argv) + 1):],
         standalone_mode=False)
