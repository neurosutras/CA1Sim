__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
from scipy import signal
import random
import sys
"""

"""
morph_filename = 'EB2-late-bifurcation.swc'
#mech_filename = '103115 interim dendritic excitability ampa nmda_kin3'
mech_filename = '112915_less_excitable'

if len(sys.argv) > 1:
    synapses_seed = int(sys.argv[1])
else:
    synapses_seed = 1
if len(sys.argv) > 2:
    num_exc_syns = int(sys.argv[2])
else:
    num_exc_syns = 3600
if len(sys.argv) > 3:
    num_inh_syns = int(sys.argv[3])
else:
    num_inh_syns = 400
# whether to modulate the peak rate of all inhibitory inputs (0 = no, 1 = out of field at track start, 2 = in field)
# input_field_width)
if len(sys.argv) > 4:
    mod_inh = int(sys.argv[4])
else:
    mod_inh = 0
# allows parallel computation of multiple trials for the same spines with the same peak_locs, but with different
# input spike trains and stochastic synapses for each trial
if len(sys.argv) > 5:
    trial_seed = int(sys.argv[5])
else:
    trial_seed = None

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())+'-seed'+\
               str(synapses_seed)+'-e'+str(num_exc_syns)+'-i'+str(num_inh_syns)+'-mod_inh'+str(mod_inh)+\
               '-cal_fb'+str(trial_seed)


def get_instantaneous_spike_probability(rate, dt=0.1, generator=None):
    """

    :param rate: float (Hz)
    :param dt: float (ms)
    :param generator: :class:'random.Random'
    :return: bool
    """
    if generator is None:
        generator = random
    x = generator.uniform(0, 1)
    rate /= 1000.
    p = 1 - np.exp(-rate * dt)
    return bool(x < p)


def get_inhom_poisson_spike_times(rate, t, dt=0.1, refractory=3., generator=None):
    """

    :param rate: instantaneous rates in time (Hz)
    :param t: corresponding time values (ms)
    :param dt: temporal resolution for spike times (ms)
    :param refractory: absolute deadtime following a spike (ms)
    :param generator: :class:'random.Random()'
    :return: list of m spike times (ms)
    """
    if generator is None:
        generator = random
    interp_t = np.arange(t[0], t[-1]+dt, dt)
    interp_rate = np.interp(interp_t, t, rate)
    spike_times = []
    i = 0
    while i < len(interp_t):
        if get_instantaneous_spike_probability(interp_rate[i], dt, generator):
            spike_times.append(interp_t[i])
            i += int(refractory / dt)
        else:
            i += 1
    return spike_times


def run_n_trials(n):
    """

    :param n: int
    """
    global trials
    for simiter in range(trials, trials + n):
        stim_exc_trains = {}
        stim_inh_trains = {}
        local_random.seed(simiter)
        global_phase_offset = local_random.uniform(-np.pi, np.pi)
        if mod_inh > 0:
            if mod_inh == 1:
                mod_inh_start = int(track_equilibrate / dt)
            elif mod_inh == 2:
                mod_inh_start = int((track_equilibrate + modulated_field_center - 0.2 * input_field_duration) / dt)
            sim.parameters['mod_inh_start'] = stim_t[mod_inh_start]
            mod_inh_stop = mod_inh_start + int(inhibitory_manipulation_duration * input_field_duration / dt)
            sim.parameters['mod_inh_stop'] = stim_t[mod_inh_stop]
        for group in stim_exc_syns.keys():  # ['CA3']:  #
            stim_exc_trains[group] = []
            for i, syn in enumerate(stim_exc_syns[group]):
                # the stochastic sequence used for each synapse is unique for each trial,
                # up to 1000 input spikes per spine
                if excitatory_stochastic:
                    syn.randObj.seq(rand_exc_seq_locs[group][i]+int(simiter*1e3))
                gauss_force = excitatory_peak_rate * np.exp(-((stim_t - peak_locs[group][i]) / gauss_sigma)**2.)
                if group == 'ECIII':
                    theta_force = excitatory_theta_offset + excitatory_theta_amp * np.cos(2. * np.pi /
                                            global_theta_cycle_duration * stim_t - global_phase_offset -
                                            excitatory_theta_phase_offset['ECIII'])
                else:
                    unit_phase_offset = peak_locs[group][i] * theta_compression_factor
                    theta_force = excitatory_theta_offset + excitatory_theta_amp * np.cos(2. * np.pi /
                                            unit_theta_cycle_duration * (stim_t - unit_phase_offset) -
                                            global_phase_offset - excitatory_theta_phase_offset['CA3'])
                stim_force = np.multiply(gauss_force, theta_force)
                train = get_inhom_poisson_spike_times(stim_force, stim_t, dt=stim_dt, generator=local_random)
                stim_exc_trains[group].append(train)
                syn.source.play(h.Vector(np.add(train, equilibrate + track_equilibrate)))
        for group in stim_inh_syns.keys():
            stim_inh_trains[group] = []
            for i, syn in enumerate(stim_inh_syns[group]):
                inhibitory_theta_amp = inhibitory_peak_rate[group] * inhibitory_theta_modulation_depth[group] / 2.
                inhibitory_theta_offset = inhibitory_peak_rate[group] - inhibitory_theta_amp
                inhibitory_phase_offset = inhibitory_theta_phase_offset[group]
                inhibitory_theta_force = inhibitory_theta_offset + inhibitory_theta_amp * np.cos(2. * np.pi /
                                                        global_theta_cycle_duration * stim_t - global_phase_offset -
                                                        inhibitory_phase_offset)
                if mod_inh > 0 and group in inhibitory_manipulation_strength:
                    inhibitory_theta_amp = inhibitory_manipulation_strength[group] * inhibitory_peak_rate[group] * \
                                           inhibitory_theta_modulation_depth[group] / 2.
                    inhibitory_theta_offset = inhibitory_manipulation_strength[group] * inhibitory_peak_rate[group] - \
                                              inhibitory_theta_amp
                    inhibitory_theta_force[mod_inh_start:mod_inh_stop] = inhibitory_theta_offset + \
                                                            inhibitory_theta_amp * np.cos(2. * np.pi /
                                                            global_theta_cycle_duration *
                                                            stim_t[mod_inh_start:mod_inh_stop] - global_phase_offset -
                                                            inhibitory_phase_offset)
                train = get_inhom_poisson_spike_times(inhibitory_theta_force, stim_t, dt=stim_dt,
                                                      generator=local_random)
                stim_inh_trains[group].append(train)
                syn.source.play(h.Vector(np.add(train, equilibrate + track_equilibrate)))
        sim.run(v_init)
        with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
            sim.export_to_file(f, simiter)
            f[str(simiter)].create_group('train')
            f[str(simiter)].create_group('inh_train')
            f[str(simiter)].create_group('successes')
            f[str(simiter)].attrs['phase_offset'] = global_phase_offset / 2. / np.pi * global_theta_cycle_duration
            index = 0
            for group in stim_exc_trains.keys():
                for i, train in enumerate(stim_exc_trains[group]):
                    f[str(simiter)]['train'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                            data=train)
                    f[str(simiter)]['train'][str(index)].attrs['group'] = group
                    f[str(simiter)]['train'][str(index)].attrs['index'] = stim_exc_syns[group][i].node.index
                    f[str(simiter)]['train'][str(index)].attrs['type'] = stim_exc_syns[group][i].node.parent.parent.type
                    if excitatory_stochastic:
                        f[str(simiter)]['successes'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                data=np.subtract(stim_exc_syns[group][i].netcon('AMPA_KIN').get_recordvec().to_python(),
                                                 equilibrate + track_equilibrate))
                    f[str(simiter)]['train'][str(index)].attrs['peak_loc'] = peak_locs[group][i]
                    index += 1
            index = 0
            for group in stim_inh_trains.keys():
                for i, train in enumerate(stim_inh_trains[group]):
                    f[str(simiter)]['inh_train'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                                data=train)
                    f[str(simiter)]['inh_train'][str(index)].attrs['group'] = group
                    f[str(simiter)]['inh_train'][str(index)].attrs['index'] = stim_inh_syns[group][i].node.index
                    f[str(simiter)]['inh_train'][str(index)].attrs['loc'] = stim_inh_syns[group][i].loc
                    f[str(simiter)]['inh_train'][str(index)].attrs['type'] = stim_inh_syns[group][i].node.type
                    index += 1
            # save the spike output of the cell, removing the equilibration offset
            f[str(simiter)].create_dataset('output', compression='gzip', compression_opts=9,
                                        data=np.subtract(cell.spike_detector.get_recordvec().to_python(),
                                                         equilibrate + track_equilibrate))
    trials += n


def plot_waveform_phase_vs_time(t, x, time_offset=0.):
    """

    :param :
    """
    phase_array = []
    peak_array = []
    peak_locs = signal.argrelmax(x)[0]
    peak_times = t[peak_locs]
    peak_array.append(peak_times)
    peak_times = np.subtract(peak_times, time_offset)
    peak_phases = np.mod(peak_times, global_theta_cycle_duration)
    peak_phases /= global_theta_cycle_duration
    peak_phases *= 360.
    phase_array.append(peak_phases)
    return peak_times, peak_phases


NMDA_type = 'NMDA_KIN3'

equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # (ms)
input_field_width = 20  # (theta cycles per 6 standard deviations)
excitatory_phase_extent = 450.  # (degrees)
# Geissler...Buzsaki, PNAS 2010
unit_theta_cycle_duration = global_theta_cycle_duration * input_field_width / (input_field_width +
                                                                               (excitatory_phase_extent / 360.))
input_field_duration = input_field_width * global_theta_cycle_duration
track_length = 2.5  # field widths
track_duration = track_length * input_field_duration
track_equilibrate = 2. * global_theta_cycle_duration
duration = equilibrate + track_equilibrate + track_duration  # input_field_duration
excitatory_peak_rate = 40.
excitatory_theta_modulation_depth = 0.8
theta_compression_factor = 1. - unit_theta_cycle_duration / global_theta_cycle_duration
excitatory_theta_phase_offset = {}
excitatory_theta_phase_offset['CA3'] = 160. / 360. * 2. * np.pi  # radians
excitatory_theta_phase_offset['ECIII'] = 0. / 360. * 2. * np.pi  # radians
excitatory_stochastic = 0
inhibitory_peak_rate = {}
inhibitory_theta_modulation_depth = {}
inhibitory_theta_phase_offset = {}
inhibitory_manipulation_strength = {'perisomatic': 0.5, 'apical dendritic': 0.8, 'tuft feedback': 0.8}
                                    # 'tuft feedforward': 1.,
inhibitory_manipulation_duration = 0.6  # Ratio of input_field_duration
inhibitory_peak_rate['perisomatic'] = 40.
inhibitory_peak_rate['apical dendritic'] = 40.
inhibitory_peak_rate['tuft feedforward'] = 40.
inhibitory_peak_rate['tuft feedback'] = 40.
inhibitory_theta_modulation_depth['perisomatic'] = 0.5
inhibitory_theta_modulation_depth['apical dendritic'] = 0.5
inhibitory_theta_modulation_depth['tuft feedforward'] = 0.5
inhibitory_theta_modulation_depth['tuft feedback'] = 0.5
inhibitory_theta_phase_offset['perisomatic'] = 145. / 360. * 2. * np.pi  # Like PV+ Basket
inhibitory_theta_phase_offset['apical dendritic'] = 215. / 360. * 2. * np.pi
                                                                            # Like Bistratified (Mixed PV+, CCK+, NPY+)
inhibitory_theta_phase_offset['tuft feedforward'] = 345. / 360. * 2. * np.pi  # Like Neurogliaform
inhibitory_theta_phase_offset['tuft feedback'] = 215. / 360. * 2. * np.pi  # Like SST+ O-LM

stim_dt = 0.02
dt = 0.02
v_init = -67.

syn_types = ['AMPA_KIN', NMDA_type]

local_random = random.Random()

# choose a subset of synapses to stimulate with inhomogeneous poisson rates
local_random.seed(synapses_seed)

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.insert_inhibitory_synapses_in_subset()

trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
if trunk_bifurcation:
    trunk_branches = [branch for branch in trunk_bifurcation[0].children if branch.type == 'trunk']
    # get where the thickest trunk branch gives rise to the tuft
    trunk = max(trunk_branches, key=lambda node: node.sec(0.).diam)
    trunk = (node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type
                                                                                    for child in node.children)).next()
else:
    trunk_bifurcation = [node for node in cell.trunk if 'tuft' in (child.type for child in node.children)]
    trunk = trunk_bifurcation[0]

all_exc_syns = {sec_type: [] for sec_type in ['basal', 'trunk', 'apical', 'tuft']}
#all_exc_syns = {sec_type: [] for sec_type in ['basal', 'trunk', 'apical']}
#all_exc_syns = {sec_type: [] for sec_type in ['tuft']}
all_inh_syns = {sec_type: [] for sec_type in ['soma', 'basal', 'trunk', 'apical', 'tuft']}
#all_inh_syns = {sec_type: [] for sec_type in ['soma', 'basal', 'trunk', 'apical']}
#all_inh_syns = {sec_type: [] for sec_type in ['tuft']}
stim_exc_syns = {'CA3': [], 'ECIII': []}
stim_inh_syns = {'perisomatic': [], 'apical dendritic': [], 'tuft feedforward': [], 'tuft feedback': []}
stim_successes = []
peak_locs = {'CA3': [], 'ECIII': []}

# place synapses in trunk for inheritance of mechanisms (for testing)
if 'trunk' not in all_exc_syns:
    for node in cell.trunk:
        for spine in node.spines:
            syn = Synapse(cell, spine, syn_types, stochastic=excitatory_stochastic)

# place synapses in every spine
for sec_type in all_exc_syns:
    for node in cell.get_nodes_of_subtype(sec_type):
        for spine in node.spines:
            syn = Synapse(cell, spine, syn_types, stochastic=excitatory_stochastic)
            all_exc_syns[sec_type].append(syn)
cell.init_synaptic_mechanisms()

# collate inhibitory synapses
for sec_type in all_inh_syns:
    for node in cell.get_nodes_of_subtype(sec_type):
        for syn in node.synapses:
            if 'GABA_A_KIN' in syn._syn:
                all_inh_syns[sec_type].append(syn)

sim = QuickSim(duration)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['track_equilibrate'] = track_equilibrate
sim.parameters['global_theta_cycle_duration'] = global_theta_cycle_duration
sim.parameters['input_field_duration'] = input_field_duration
sim.parameters['track_length'] = track_length
sim.parameters['duration'] = duration
sim.parameters['stim_dt'] = stim_dt
sim.append_rec(cell, cell.tree.root, description='soma', loc=0.5)
sim.append_rec(cell, trunk, description='distal_trunk', loc=0.)
sim.append_rec(cell, trunk_bifurcation[0], description='proximal_trunk', loc=1.)
spike_output_vec = h.Vector()
cell.spike_detector.record(spike_output_vec)

# get the fraction of total spines contained in each sec_type

total_exc_syns = {sec_type: len(all_exc_syns[sec_type]) for sec_type in ['basal', 'trunk', 'apical', 'tuft']}
fraction_exc_syns = {sec_type: float(total_exc_syns[sec_type]) / float(np.sum(total_exc_syns.values())) for sec_type in
                 ['basal', 'trunk', 'apical', 'tuft']}
"""
total_exc_syns = {sec_type: len(all_exc_syns[sec_type]) for sec_type in ['basal', 'trunk', 'apical']}
fraction_exc_syns = {sec_type: float(total_exc_syns[sec_type]) / float(np.sum(total_exc_syns.values())) for sec_type in
                 ['basal', 'trunk', 'apical']}

total_exc_syns = {sec_type: len(all_exc_syns[sec_type]) for sec_type in ['tuft']}
fraction_exc_syns = {sec_type: float(total_exc_syns[sec_type]) / float(np.sum(total_exc_syns.values())) for sec_type in
                 ['tuft']}
"""
for sec_type in all_exc_syns:
    for i in local_random.sample(range(len(all_exc_syns[sec_type])), int(num_exc_syns*fraction_exc_syns[sec_type])):
        syn = all_exc_syns[sec_type][i]
        if sec_type == 'tuft':
            stim_exc_syns['ECIII'].append(syn)
        else:
            stim_exc_syns['CA3'].append(syn)

# get the fraction of inhibitory synapses contained in each sec_type

total_inh_syns = {sec_type: len(all_inh_syns[sec_type]) for sec_type in ['soma', 'basal', 'trunk', 'apical', 'tuft']}
fraction_inh_syns = {sec_type: float(total_inh_syns[sec_type]) / float(np.sum(total_inh_syns.values())) for sec_type in
                 ['soma', 'basal', 'trunk', 'apical', 'tuft']}
num_inh_syns = min(num_inh_syns, int(np.sum(total_inh_syns.values())))
"""
total_inh_syns = {sec_type: len(all_inh_syns[sec_type]) for sec_type in ['soma', 'basal', 'trunk', 'apical']}
fraction_inh_syns = {sec_type: float(total_inh_syns[sec_type]) / float(np.sum(total_inh_syns.values())) for sec_type in
                 ['soma', 'basal', 'trunk', 'apical']}

total_inh_syns = {sec_type: len(all_inh_syns[sec_type]) for sec_type in ['tuft']}
fraction_inh_syns = {sec_type: float(total_inh_syns[sec_type]) / float(np.sum(total_inh_syns.values())) for sec_type in
                 ['tuft']}
"""
for sec_type in all_inh_syns:
    for i in local_random.sample(range(len(all_inh_syns[sec_type])), int(num_inh_syns*fraction_inh_syns[sec_type])):
        syn = all_inh_syns[sec_type][i]
        if syn.node.type == 'tuft':
            if cell.is_terminal(syn.node):
                # GABAergic synapses on terminal tuft branches are about 25% feedforward
                group = local_random.choice(['tuft feedforward', 'tuft feedback', 'tuft feedback', 'tuft feedback'])
            else:
                # GABAergic synapses on intermediate tuft branches are about 50% feedforward
                group = local_random.choice(['tuft feedforward', 'tuft feedback'])
        else:
            distance = cell.get_distance_to_node(cell.tree.root, syn.node, syn.loc)
            group = 'perisomatic' if distance <= 75. else 'apical dendritic'
        stim_inh_syns[group].append(syn)

stim_t = np.arange(-track_equilibrate, track_duration, dt)

excitatory_theta_amp = excitatory_theta_modulation_depth / 2.
excitatory_theta_offset = 1. - excitatory_theta_amp
gauss_sigma = global_theta_cycle_duration * input_field_width / 3. / np.sqrt(2.)  # contains 99.7% gaussian area

rand_exc_seq_locs = {}
for group in stim_exc_syns.keys():
    rand_exc_seq_locs[group] = []
    if stim_exc_syns[group]:
        peak_locs[group] = np.arange(-0.75 * input_field_duration, (0.75 + track_length) * input_field_duration,
                          (1.5 + track_length) * input_field_duration / int(len(stim_exc_syns[group])))
        peak_locs[group] = peak_locs[group][:len(stim_exc_syns[group])]
    local_random.shuffle(peak_locs[group])
    peak_locs[group] = list(peak_locs[group])

for group in stim_exc_syns.keys():
    for syn in stim_exc_syns[group]:
        #peak_loc = local_random.uniform(-0.75 * input_field_duration, (0.75 + track_length) * input_field_duration)
        #peak_locs.append(peak_loc)
        if excitatory_stochastic:
            success_vec = h.Vector()
            stim_successes.append(success_vec)
            syn.netcon('AMPA_KIN').record(success_vec)
            rand_exc_seq_locs[group].append(syn.randObj.seq())
        #if syn.node.parent.parent not in [rec['node'] for rec in sim.rec_list]:
        #    sim.append_rec(cell, syn.node.parent.parent)
        # remove this synapse from the pool, so that additional "modulated" inputs can be selected from those that remain
        sim.append_rec(cell, syn.node, object=syn.target('AMPA_KIN'), param='_ref_i', description='i_AMPA')
        sim.append_rec(cell, syn.node, object=syn.target(NMDA_type), param='_ref_i', description='i_NMDA')
        all_exc_syns[syn.node.parent.parent.type].remove(syn)

# rand_inh_seq_locs = [] will need this when inhibitory synapses become stochastic
# stim_inh_successes = [] will need this when inhibitory synapses become stochastic

# modulate the weights of inputs that have peak_locs along this stretch of the track
#modulated_num_exc_syn = 100
modulated_field_center = track_duration * 0.6
gauss_mod_amp = {}
#modulated_start_index = len(stim_exc_syns)
for group in stim_exc_syns.keys():
    gauss_mod_amp[group] = 1.5 * np.exp(-((np.array(peak_locs[group]) - modulated_field_center) /
                                          (gauss_sigma * 1.4)) ** 2.) + 1.
    for i, syn in enumerate(stim_exc_syns[group]):
        syn.netcon('AMPA_KIN').weight[0] = gauss_mod_amp[group][i]
        """
        if (modulated_field_center - input_field_duration * 0.75 <= peak_locs[group][i] <
                    modulated_field_center - input_field_duration * 0.5) or (modulated_field_center +
                    input_field_duration * 0.5 < peak_locs[group][i] <= modulated_field_center +
                    input_field_duration * 0.75):
            syn.netcon('AMPA_KIN').weight[0] = 0.75
        elif (modulated_field_center - input_field_duration * 0.25 <= peak_locs[group][i] <= modulated_field_center +
                    input_field_duration * 0.25):
            syn.netcon('AMPA_KIN').weight[0] = 3.5
        """
"""
for sec_type in all_exc_syns:
    for i in local_random.sample(range(len(all_exc_syns[sec_type])), min(len(all_exc_syns[sec_type]),
                                                            int(modulated_num_exc_syn * fraction_exc_syns[sec_type]))):
        syn = all_exc_syns[sec_type][i]
        stim_exc_syns.append(syn)

for syn in stim_exc_syns[modulated_start_index:]:
    peak_loc = modulated_field_center
    peak_locs.append(peak_loc)
    if excitatory_stochastic:
        success_vec = h.Vector()
        stim_successes.append(success_vec)
        syn.netcon('AMPA_KIN').record(success_vec)
        rand_exc_seq_locs.append(syn.randObj.seq())
    # sim.append_rec(cell, syn.node, object=syn.target('AMPA_KIN'), param='_ref_i', description='i_AMPA')
    # sim.append_rec(cell, syn.node, object=syn.target(NMDA_type), param='_ref_i', description='i_NMDA')
    #if syn.node.parent.parent not in [rec['node'] for rec in sim.rec_list]:
    #    sim.append_rec(cell, syn.node.parent.parent)
    # remove this synapse from the pool, so that additional "modulated" inputs can be selected from those that remain
    all_exc_syns[syn.node.parent.parent.type].remove(syn)
"""

if trial_seed is None:
    trials = 0
    run_n_trials(1)
else:
    trials = trial_seed
    run_n_trials(1)

"""
global_phase_offset = 0.
end_baseline = np.where(stim_t >= input_field_duration / 2.)[0][0]
stim_forces = {'CA3': [], 'ECIII': []}
force_sum = {}
envelope = {}
mean_amp = {}
global_cos = {}
peak_times = {}
peak_phases = {}
for group in stim_exc_syns.keys():
    for i, syn in enumerate(stim_exc_syns[group]):
        gauss_force = excitatory_peak_rate * np.exp(-((stim_t - peak_locs[group][i]) / gauss_sigma)**2.)
        gauss_force *= gauss_mod_amp[group][i]
        if group == 'ECIII':
            theta_force = excitatory_theta_offset + excitatory_theta_amp * np.cos(2. * np.pi /
                                    global_theta_cycle_duration * stim_t - global_phase_offset -
                                    excitatory_theta_phase_offset['ECIII'])
        else:
            unit_phase_offset = peak_locs[group][i] * theta_compression_factor
            theta_force = excitatory_theta_offset + excitatory_theta_amp * np.cos(2. * np.pi /
                            unit_theta_cycle_duration * (stim_t - unit_phase_offset) - global_phase_offset -
                            excitatory_theta_phase_offset['CA3'])
        stim_force = np.multiply(gauss_force, theta_force)
        stim_forces[group].append(stim_force)
    force_sum[group] = np.sum(stim_forces[group], 0)
    mean_amp[group] = np.mean(force_sum[group][:end_baseline])
    force_sum[group] -= mean_amp[group]
    envelope[group] = np.mean(np.abs(signal.hilbert(force_sum[group][:end_baseline])))
    force_sum[group] += mean_amp[group]
    global_cos[group] = envelope[group] * np.cos(2. * np.pi / global_theta_cycle_duration * stim_t) + mean_amp[group]
fig, axes = plt.subplots(2, 3)
force_sum['Total'] = np.add(force_sum['CA3'], force_sum['ECIII'])
mean_amp['Total'] = np.mean(force_sum['Total'][:end_baseline])
force_sum['Total'] -= mean_amp['Total']
envelope['Total'] = np.mean(np.abs(signal.hilbert(force_sum['Total'][:end_baseline])))
force_sum['Total'] += mean_amp['Total']
global_cos['Total'] = envelope['Total'] * np.cos(2. * np.pi / global_theta_cycle_duration * stim_t) + mean_amp['Total']
axes[0][0].plot(stim_t, force_sum['CA3'])
axes[0][0].plot(stim_t, global_cos['CA3'])
axes[0][1].plot(stim_t, force_sum['ECIII'])
axes[0][1].plot(stim_t, global_cos['ECIII'])
axes[0][2].plot(stim_t, force_sum['Total'])
axes[0][2].plot(stim_t, global_cos['Total'])
for group in force_sum.keys():
    peak_times[group], peak_phases[group] = plot_waveform_phase_vs_time(stim_t, force_sum[group],
                                            time_offset=global_phase_offset / 2. / np.pi * global_theta_cycle_duration)
axes[1][0].scatter(peak_times['CA3'], peak_phases['CA3'])
axes[1][1].scatter(peak_times['CA3'], peak_phases['ECIII'])
axes[1][2].scatter(peak_times['CA3'], peak_phases['Total'])
plt.show()
plt.close()
"""