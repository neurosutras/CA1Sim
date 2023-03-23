__author__ = 'Aaron D. Milstein'
"""
In this version of the simulation, phase precession of CA3 inputs is implemented using the method from Chadwick et al.,
Elife, 2015, which uses a circular gaussian with a phase sensitivity factor that effectively compresses the range of
phases within each theta cycle that each input is active, which will reduce jitter across within-cycle input sequences.

"""
from specify_cells import *
import click


@click.command()
@click.option("--mech-filename", type=str, default='20220808_default_biophysics.yaml')
@click.option("--synapses-seed", type=int, default=0)
@click.option("--trial-seed", type=int, default=0)
@click.option("--label", type=str, default=None)
@click.option("--mod_inh", type=int, default=0)
@click.option("--sim_duration", type=float, default=None)
@click.option("--field_center", type=float, default=0.6)
@click.option("--spines", type=bool, default=True)
@click.option("--export", is_flag=True)
@click.option("--plot", is_flag=True)
@click.option("--interactive", is_flag=True)
@click.option("--debug", is_flag=True)
def main(mech_filename, synapses_seed, trial_seed, label, mod_inh, sim_duration, field_center, spines, export, plot,
         interactive, debug):
    """

    :param mech_filename:
        .yaml file must be located in the data subdirectory
    :param synapses_seed:
        a unique random seed can be used to shuffle the number and locations of synapses, and the place field locs of
        the presynaptic CA3 inputs (like simulating a different cell with the same morphology)
    :param trial_seed:
        a unique random seed shuffles the input spike times and synaptic release probabilities to allows simulation of
        multiple independent trials
    :param label: append a label when exporting data to .hdf5
    :param mod_inh:
        whether to decrease the firing rate of inhibitory inputs, mimicking the optogenetic silencing in
        Grienberger, Milstein et al., Nat. Neurosci., 2017. (0 = no, 1 = out of field at track start, 2 = in field,
        3 = entire length of track)
    :param sim_duration: float (ms) - sim duration can be truncated during testing
    :param field_center: float, value between 0 and 1, where along the track the CA1 place field peaks
    :param spines: bool, whether to include explicit spine neck and head compartments for every excitatory synapse
    :param export: bool, whether to export to .hdf5
    :param plot: bool, whether to plot
    :param interactive: bool, whether to enable live object inspection after simulation
    :param debug: bool, does not run simulation in debug mode
    """
    morph_filename = 'EB2-late-bifurcation.swc'

    num_exc_syns = 3200
    num_inh_syns = 600

    # the synaptic AMPAR conductances at in-field inputs are multiplied by a factor with this value at the peak of the
    # field, and decays with cosine spatial modulation away from the field
    mod_weights = 2.5

    if label is None:
        rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())+'-seed'+\
                       str(synapses_seed)+'-e'+str(num_exc_syns)+'-i'+str(num_inh_syns)+'-mod_inh'+str(mod_inh)+\
                       '-trial'+str(trial_seed) + '.hdf5'
    else:
        rec_filename = 'output' + datetime.datetime.today().strftime('%m%d%Y%H%M') + '-pid' + str(os.getpid()) + \
                       '-' + label + '-seed' + str(synapses_seed) + '-e' + str(num_exc_syns) + '-i' + str(num_inh_syns) + \
                       '-mod_inh' + str(mod_inh) + '-trial' + str(trial_seed) + '.hdf5'


    def get_dynamic_theta_phase_force(phase_ranges, peak_loc, input_field_duration, stim_t, dt):
        """
        Expects a list of tuples containing times and phases relative to peak_loc and the non-modulated phase preference
        (zero degrees). Returns a waveform of phase vs time.
        :param phase_ranges: list of tuple (ms, degrees)
        :param peak_loc:
        :param input_field_duration:
        :param stim_t:
        :param dt:
        :return: :class: 'np.array'
        """
        start_phase_val = phase_ranges[0][1] * 2. * np.pi / 360.  # convert degrees to radians
        end_phase_val = phase_ranges[-1][1] * 2. * np.pi / 360.  # convert degrees to radians
        phase_force = np.ones_like(stim_t) * start_phase_val
        phase_gradient = np.array([])
        for i in range(len(phase_ranges)-1):
            t0 = phase_ranges[i][0]
            t1 = phase_ranges[i+1][0]
            phase0 = phase_ranges[i][1] * 2. * np.pi / 360.  # convert degrees to radians
            phase1 = phase_ranges[i+1][1] * 2. * np.pi / 360.
            del_t = t1 - t0
            del_phase = phase1 - phase0
            if abs(del_phase) > 0.:
                del_phase = del_phase / del_t * dt
                this_range_piece = np.arange(phase0, phase1, del_phase)
            else:
                this_range_piece = np.ones(int(del_t / dt)) * phase0
            phase_gradient = np.append(phase_gradient, this_range_piece)
        if stim_t[0] <= peak_loc-input_field_duration*0.5 <= stim_t[-1]:
            phase_start = np.where(peak_loc-input_field_duration*0.5 >= stim_t)[0]
            if np.any(phase_start):
                phase_start = phase_start[-1]
                phase_end = min(len(stim_t), phase_start+len(phase_gradient))
                phase_force[:phase_start] = start_phase_val
                phase_force[phase_start:phase_end] = phase_gradient[:phase_end-phase_start]
                phase_force[phase_end:] = end_phase_val
        elif stim_t[0] <= peak_loc+input_field_duration*0.5 <= stim_t[-1]:
            phase_end = np.where(peak_loc+input_field_duration*0.5 >= stim_t)[0]
            if np.any(phase_end):
                phase_end = phase_end[-1]
                phase_start = max(0, phase_end-len(phase_gradient))
                phase_force[:phase_start] = start_phase_val
                phase_force[phase_start:phase_end] = phase_gradient[-(phase_end-phase_start):]
                phase_force[phase_end:] = end_phase_val
        return phase_force


    def run_trial(simiter):
        """

        :param simiter: int
        """
        local_random.seed(simiter)
        global_phase_offset = local_random.uniform(-np.pi, np.pi)
        if not debug and export:
            with h5py.File(data_dir+rec_filename, 'a') as f:
                f.create_group(str(simiter))
                f[str(simiter)].create_group('train')
                f[str(simiter)].create_group('inh_train')
                f[str(simiter)].attrs['phase_offset'] = global_phase_offset / 2. / np.pi * global_theta_cycle_duration
        exc_rate_maps = {}
        if mod_inh > 0:
            if mod_inh == 1:
                mod_inh_start = int(track_equilibrate / dt)
                mod_inh_stop = mod_inh_start + int(inhibitory_manipulation_duration * input_field_duration / dt)
            elif mod_inh == 2:
                mod_inh_start = int((track_equilibrate + modulated_field_center - 0.3 * input_field_duration) / dt)
                mod_inh_stop = mod_inh_start + int(inhibitory_manipulation_duration * input_field_duration / dt)
            elif mod_inh == 3:
                mod_inh_start = 0
                mod_inh_stop = len(stim_t)
            sim.parameters['mod_inh_start'] = stim_t[mod_inh_start]
            sim.parameters['mod_inh_stop'] = stim_t[mod_inh_stop-1]
        index = 0
        for group in stim_exc_syns:
            exc_rate_maps[group] = []
            for i, syn in enumerate(stim_exc_syns[group]):
                # the stochastic sequence used for each synapse is unique for each trial,
                # up to 1000 input spikes per spine
                if excitatory_stochastic:
                    syn.randObj.seq(rand_exc_seq_locs[group][i]+int(simiter*1e3))
                gauss_force = excitatory_peak_rate[group] * np.exp(-((stim_t - peak_locs[group][i]) / gauss_sigma)**2.)
                if group in excitatory_precession_range:
                    phase_force = get_dynamic_theta_phase_force(excitatory_precession_range[group], peak_locs[group][i],
                                                                input_field_duration, stim_t, stim_dt)
                    theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] * np.cos(phase_force +
                                            excitatory_theta_phase_offset[group] - 2. * np.pi * stim_t /
                                            global_theta_cycle_duration + global_phase_offset))
                else:
                    theta_force = np.exp(excitatory_theta_phase_tuning_factor[group] *
                                     np.cos(excitatory_theta_phase_offset[group] - 2. * np.pi * stim_t /
                                            global_theta_cycle_duration + global_phase_offset))
                theta_force -= np.min(theta_force)
                theta_force /= np.max(theta_force)
                theta_force *= excitatory_theta_modulation_depth[group]
                theta_force += 1. - excitatory_theta_modulation_depth[group]
                stim_force = np.multiply(gauss_force, theta_force)
                if not debug:
                    train = get_inhom_poisson_spike_times_by_thinning(stim_force, stim_t, dt=stim_dt,
                                                                      generator=local_random)
                    syn.source.play(h.Vector(np.add(train, equilibrate + track_equilibrate)))
                    if export:
                        with h5py.File(data_dir+rec_filename, 'a') as f:
                            f[str(simiter)]['train'].create_dataset(str(index), compression='gzip', compression_opts=9, data=train)
                            f[str(simiter)]['train'][str(index)].attrs['group'] = group
                            f[str(simiter)]['train'][str(index)].attrs['index'] = syn.node.index
                            f[str(simiter)]['train'][str(index)].attrs['type'] = syn.node.parent.parent.type
                            f[str(simiter)]['train'][str(index)].attrs['peak_loc'] = peak_locs[group][i]
                else:
                    exc_rate_maps[group].append(stim_force)
                index += 1
        if not debug:
            index = 0
            for group in stim_inh_syns:
                inh_peak_rate = 2. * inhibitory_mean_rate[group] / (2. - inhibitory_theta_modulation_depth[group])
                inhibitory_theta_force = np.exp(inhibitory_theta_phase_tuning_factor[group] *
                                                np.cos(inhibitory_theta_phase_offset[group] - 2. * np.pi * stim_t /
                                                       global_theta_cycle_duration + global_phase_offset))
                inhibitory_theta_force -= np.min(inhibitory_theta_force)
                inhibitory_theta_force /= np.max(inhibitory_theta_force)
                inhibitory_theta_force *= inhibitory_theta_modulation_depth[group]
                inhibitory_theta_force += 1. - inhibitory_theta_modulation_depth[group]
                inhibitory_theta_force *= inh_peak_rate
                for syn in stim_inh_syns[group]:
                    stim_force = np.array(inhibitory_theta_force)
                    if mod_inh > 0 and group in inhibitory_manipulation_offset:
                        # inhibitory manipulation subtracts from the mean firing rate, but maintains the same theta modulation
                        # depth
                        mod_inh_multiplier = 1. - inhibitory_manipulation_offset[group] / inhibitory_mean_rate[group]
                        stim_force[mod_inh_start:mod_inh_stop] *= mod_inh_multiplier
                    train = get_inhom_poisson_spike_times_by_thinning(stim_force, stim_t, dt=stim_dt,
                                                                      generator=local_random)
                    syn.source.play(h.Vector(np.add(train, equilibrate + track_equilibrate)))
                    if export:
                        with h5py.File(data_dir+rec_filename, 'a') as f:
                            f[str(simiter)]['inh_train'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                                        data=train)
                            f[str(simiter)]['inh_train'][str(index)].attrs['group'] = group
                            f[str(simiter)]['inh_train'][str(index)].attrs['index'] = syn.node.index
                            f[str(simiter)]['inh_train'][str(index)].attrs['loc'] = syn.loc
                            f[str(simiter)]['inh_train'][str(index)].attrs['type'] = syn.node.type
                    index += 1
            sim.run(v_init)
            if export:
                sim.export_to_file(data_dir + rec_filename, simiter)
                with h5py.File(data_dir+rec_filename, 'a') as f:
                    if excitatory_stochastic:
                        f[str(simiter)].create_group('successes')
                        index = 0
                        for group in stim_exc_syns:
                            for syn in stim_exc_syns[group]:
                                f[str(simiter)]['successes'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                            data=np.subtract(syn.netcon('AMPA_KIN').get_recordvec().to_python(),
                                                             equilibrate + track_equilibrate))
                                index += 1
                    # save the spike output of the cell, removing the equilibration offset
                    f[str(simiter)].create_dataset('output', compression='gzip', compression_opts=9,
                                                data=np.subtract(cell.spike_detector.get_recordvec().to_python(),
                                                                 equilibrate + track_equilibrate))
        if debug:
            return exc_rate_maps

    NMDA_type = 'NMDA_KIN5'

    equilibrate = 250.  # time to steady-state

    global_theta_cycle_duration = 150.  # (ms)
    input_field_width = 20  # (theta cycles per 6 standard deviations)
    input_field_duration = input_field_width * global_theta_cycle_duration
    track_length = 2.5  # field widths
    track_duration = track_length * input_field_duration
    track_equilibrate = 2. * global_theta_cycle_duration
    if sim_duration is None:
        duration = equilibrate + track_equilibrate + track_duration
    else:
        duration = equilibrate + sim_duration

    excitatory_peak_rate = {'CA3': 40., 'ECIII': 40.}
    excitatory_theta_modulation_depth = {'CA3': 0.7, 'ECIII': 0.7}
    # From Chadwick et al., ELife 2015
    excitatory_theta_phase_tuning_factor = {'CA3': 0.8, 'ECIII': 0.8}
    excitatory_precession_range = {}  # (ms, degrees)
    excitatory_precession_range['CA3'] = [(-input_field_duration*0.5, 180.), (-input_field_duration*0.35, 180.),
                                          (input_field_duration*0.35, -180.), (input_field_duration*0.5, -180.)]
    excitatory_theta_phase_offset = {}
    excitatory_theta_phase_offset['CA3'] = 165. / 360. * 2. * np.pi  # radians
    excitatory_theta_phase_offset['ECIII'] = 0. / 360. * 2. * np.pi  # radians
    excitatory_stochastic = 1
    inhibitory_manipulation_offset = {'perisomatic': 9., 'axo-axonic': 9., 'apical dendritic': 9.,
                                        'distal apical dendritic': 9., 'tuft feedback': 9.}
    inhibitory_manipulation_duration = 0.6  # Ratio of input_field_duration
    inhibitory_mean_rate = {'perisomatic': 25., 'axo-axonic': 25., 'apical dendritic': 25., 'distal apical dendritic': 25.,
                            'tuft feedforward': 25., 'tuft feedback': 25.}
    inhibitory_theta_modulation_depth = {'perisomatic': 0.5, 'axo-axonic': 0.5, 'apical dendritic': 0.5,
                                         'distal apical dendritic': 0.5, 'tuft feedforward': 0.5, 'tuft feedback': 0.5}
    inhibitory_theta_phase_tuning_factor = {'perisomatic': 0.6, 'axo-axonic': 0.6, 'apical dendritic': 0.6,
                                         'distal apical dendritic': 0.6, 'tuft feedforward': 0.6, 'tuft feedback': 0.6}
    inhibitory_precession_range = {}
    inhibitory_theta_phase_offset = {}
    inhibitory_theta_phase_offset['perisomatic'] = 135. / 360. * 2. * np.pi  # Like PV+ Basket
    inhibitory_theta_phase_offset['axo-axonic'] = 45. / 360. * 2. * np.pi  # Vargas et al., ELife, 2014
    inhibitory_theta_phase_offset['apical dendritic'] = 200. / 360. * 2. * np.pi  # Like PYR-layer Bistratified
    inhibitory_theta_phase_offset['distal apical dendritic'] = 180. / 360. * 2. * np.pi  # Like SR/SLM Border Cells
    inhibitory_theta_phase_offset['tuft feedforward'] = 340. / 360. * 2. * np.pi  # Like Neurogliaform
    inhibitory_theta_phase_offset['tuft feedback'] = 200. / 360. * 2. * np.pi  # Like SST+ O-LM
    inhibitory_stochastic = 0

    stim_dt = 0.025
    dt = 0.025
    v_init = -67.

    exc_syn_types = ['AMPA_KIN', NMDA_type]
    inh_syn_types = ['GABA_A_KIN']

    local_random = random.Random()

    # choose a subset of synapses to stimulate with inhomogeneous poisson rates
    local_random.seed(synapses_seed)

    cell = CA1_Pyr(morph_filename, mech_filename, full_spines=spines)
    cell.set_terminal_branch_na_gradient()

    trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
    if trunk_bifurcation:
        trunk_branches = [branch for branch in trunk_bifurcation[0].children if branch.type == 'trunk']
        # get where the thickest trunk branch gives rise to the tuft
        trunk = max(trunk_branches, key=lambda node: node.sec(0.).diam)
        trunk = next((node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type
                                                                                        for child in node.children)))
    else:
        trunk_bifurcation = [node for node in cell.trunk if 'tuft' in (child.type for child in node.children)]
        trunk = trunk_bifurcation[0]

    exc_syns_sec_types = ['basal', 'trunk', 'apical', 'tuft']
    inh_syns_sec_types = ['soma', 'ais', 'basal', 'trunk', 'apical', 'tuft']
    stim_exc_syns = {'CA3': [], 'ECIII': []}
    stim_inh_syns = {'perisomatic': [], 'axo-axonic': [], 'apical dendritic': [], 'distal apical dendritic': [],
                     'tuft feedforward': [], 'tuft feedback': []}
    stim_successes = []
    peak_locs = {'CA3': [], 'ECIII': []}

    if spines:
        exc_syn_locs_by_sec_type = {}
        for sec_type in exc_syns_sec_types:
            exc_syn_locs_by_sec_type[sec_type] = []
            for node in cell.get_nodes_of_subtype(sec_type):
                for spine in node.spines:
                    exc_syn_locs_by_sec_type[sec_type].append((spine, 0.5))
    else:
        exc_syn_locs_by_sec_type = cell.get_excitatory_syn_locs(sec_type_list=exc_syns_sec_types)

    # get the fraction of excitatory synapses contained in each sec_type
    exc_syn_count = {sec_type: len(exc_syn_locs_by_sec_type[sec_type]) for sec_type in exc_syn_locs_by_sec_type}
    total_exc_syns = np.sum(list(exc_syn_count.values()))
    fraction_exc_syns = {sec_type: float(exc_syn_count[sec_type]) / float(total_exc_syns) for sec_type in exc_syn_count}

    for sec_type in exc_syn_locs_by_sec_type:
        exc_syn_locs = local_random.sample(exc_syn_locs_by_sec_type[sec_type],
                                           int(num_exc_syns*fraction_exc_syns[sec_type]))
        syn_list = cell.insert_synapses_at_syn_locs(exc_syn_locs, exc_syn_types, stochastic=excitatory_stochastic)
        if sec_type == 'tuft':
            stim_exc_syns['ECIII'].extend(syn_list)
        else:
            stim_exc_syns['CA3'].extend(syn_list)

    inh_syn_locs_by_sec_type = cell.get_inhibitory_syn_locs(sec_type_list=inh_syns_sec_types)

    # get the fraction of inhibitory synapses contained in each sec_type
    inh_syn_count = {sec_type: len(inh_syn_locs_by_sec_type[sec_type]) for sec_type in inh_syn_locs_by_sec_type}
    total_inh_syns = np.sum(list(inh_syn_count.values()))
    fraction_inh_syns = {sec_type: float(inh_syn_count[sec_type]) / float(total_inh_syns) for sec_type in inh_syn_count}

    for sec_type in inh_syn_locs_by_sec_type:
        inh_syn_locs = local_random.sample(inh_syn_locs_by_sec_type[sec_type],
                                           int(num_inh_syns * fraction_inh_syns[sec_type]))
        syn_list = cell.insert_synapses_at_syn_locs(inh_syn_locs, inh_syn_types, stochastic=inhibitory_stochastic)

        if sec_type == 'tuft':
            for syn in syn_list:
                if cell.is_terminal(syn.node):
                    # GABAergic synapses on terminal tuft branches are about 25% feedforward
                    group = local_random.choice(
                        ['tuft feedforward', 'tuft feedback', 'tuft feedback', 'tuft feedback'])
                else:
                    # GABAergic synapses on intermediate tuft branches are about 50% feedforward
                    group = local_random.choice(['tuft feedforward', 'tuft feedback'])
                stim_inh_syns[group].append(syn)
        elif sec_type == 'trunk':
            for syn in syn_list:
                distance = cell.get_distance_to_node(cell.tree.root, syn.node, syn.loc)
                if distance <= 50.:
                    group = 'perisomatic'
                elif distance <= 150.:
                    group = local_random.choice(['apical dendritic', 'apical dendritic', 'distal apical dendritic'])
                else:
                    group = local_random.choice(
                        ['apical dendritic', 'distal apical dendritic', 'distal apical dendritic'])
                stim_inh_syns[group].append(syn)
        elif sec_type == 'basal':
            for syn in syn_list:
                distance = cell.get_distance_to_node(cell.tree.root, syn.node, syn.loc)
                group = 'perisomatic' if distance <= 50. and not cell.is_terminal(syn.node) else 'apical dendritic'
                stim_inh_syns[group].append(syn)
        elif sec_type == 'soma':
            group = 'perisomatic'
            stim_inh_syns[group].extend(syn_list)
        elif sec_type == 'apical':
            for syn in syn_list:
                distance = cell.get_distance_to_node(cell.tree.root, cell.get_dendrite_origin(syn.node), loc=1.)
                if distance <= 150.:
                    group = local_random.choice(['apical dendritic', 'apical dendritic', 'distal apical dendritic'])
                else:
                    group = local_random.choice(
                        ['apical dendritic', 'distal apical dendritic', 'distal apical dendritic'])
                stim_inh_syns[group].append(syn)
        elif sec_type == 'ais':
            group = 'axo-axonic'
            stim_inh_syns[group].extend(syn_list)

    cell.init_synaptic_mechanisms()

    sim = QuickSim(duration, cvode=False, dt=dt)
    sim.parameters['equilibrate'] = equilibrate
    sim.parameters['track_equilibrate'] = track_equilibrate
    sim.parameters['global_theta_cycle_duration'] = global_theta_cycle_duration
    sim.parameters['input_field_duration'] = input_field_duration
    sim.parameters['track_length'] = track_length
    sim.parameters['duration'] = duration
    sim.parameters['stim_dt'] = stim_dt
    sim.append_rec(cell, cell.tree.root, description='soma', loc=0.)
    sim.append_rec(cell, trunk_bifurcation[0], description='proximal_trunk', loc=1.)
    sim.append_rec(cell, trunk, description='distal_trunk', loc=1.)
    spike_output_vec = h.Vector()
    cell.spike_detector.record(spike_output_vec)

    stim_t = np.arange(-track_equilibrate, track_duration, dt)

    gauss_sigma = global_theta_cycle_duration * input_field_width / 3. / np.sqrt(2.)  # contains 99.7% gaussian area

    rand_exc_seq_locs = {}
    for group in stim_exc_syns:
        rand_exc_seq_locs[group] = []
        if stim_exc_syns[group]:
            peak_locs[group] = np.arange(-0.75 * input_field_duration, (0.75 + track_length) * input_field_duration,
                              (1.5 + track_length) * input_field_duration / int(len(stim_exc_syns[group])))
            peak_locs[group] = peak_locs[group][:len(stim_exc_syns[group])]

    for group in stim_exc_syns:
        for syn in stim_exc_syns[group]:
            if excitatory_stochastic:
                success_vec = h.Vector()
                stim_successes.append(success_vec)
                syn.netcon('AMPA_KIN').record(success_vec)
                rand_exc_seq_locs[group].append(syn.randObj.seq())


    # modulate the weights of inputs with peak_locs along this stretch of the track
    modulated_field_center = track_duration * field_center
    cos_mod_weight = {}
    peak_mod_weight = mod_weights
    tuning_amp = (peak_mod_weight - 1.) / 2.
    tuning_offset = tuning_amp + 1.

    for group in stim_exc_syns:
        this_cos_mod_weight = tuning_amp * np.cos(2. * np.pi / (input_field_duration * 1.2) * (peak_locs[group] -
                                                                            modulated_field_center)) + tuning_offset
        left = np.where(peak_locs[group] >= modulated_field_center - input_field_duration * 1.2 / 2.)[0][0]
        right = np.where(peak_locs[group] > modulated_field_center + input_field_duration * 1.2 / 2.)[0][0]
        cos_mod_weight[group] = np.array(this_cos_mod_weight)
        cos_mod_weight[group][:left] = 1.
        cos_mod_weight[group][right:] = 1.
        peak_locs[group] = list(peak_locs[group])
        cos_mod_weight[group] = list(cos_mod_weight[group])
        indexes = list(range(len(peak_locs[group])))
        local_random.shuffle(indexes)
        peak_locs[group] = list(map(peak_locs[group].__getitem__, indexes))
        cos_mod_weight[group] = list(map(cos_mod_weight[group].__getitem__, indexes))
        for i, syn in enumerate(stim_exc_syns[group]):
            syn.netcon('AMPA_KIN').weight[0] = cos_mod_weight[group][i]

    if not debug:
        run_trial(trial_seed)
        if plot:
            sim.plot()
    else:
        exc_rate_maps = run_trial(trial_seed)
        if plot:
            for pop in ['CA3', 'ECIII']:
                fig = plt.figure()
                plt.imshow(exc_rate_maps[pop], aspect='auto', interpolation='none',
                           extent=(-track_equilibrate/1000., track_duration/1000., len(exc_rate_maps[pop])+0.5, -0.5))
                plt.xlabel('Time (sec)')
                plt.ylabel('Presynaptic unit ID')
                fig.suptitle('%s Firing Rates' % pop)
                plt.colorbar()
                fig.show()

    if interactive:
        globals().update(locals())


if __name__ == '__main__':
    main(standalone_mode=False)