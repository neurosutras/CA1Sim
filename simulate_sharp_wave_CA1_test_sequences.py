__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import random
import sys

"""
Phase precession of CA3 inputs is implemented using the method from Chadwick et al., Elife, 2015, which uses a circular
gaussian with a phase sensitivity factor that effectively compresses the range of phases within each theta cycle that
each input is active.

Track wraps around, so that place field inputs with peak locations near the track edges are active on both ends.

The peak location, determined by a unimodal spatial weight distribution, is set from the command line.
Random seeds work as follows:
    gid: pick different synapse numbers and locations for each cell (0-10)
    synapse_seed: shuffle peak_locs uniquely for each cell (gid+100)
    trial_seed_offset: set theta phase offset for all cells for one trial
        trial_seed = gid * 1000 + trial_seed_offset. a more realistic approach would be to precompute spike trains for
        many CA3 inputs, and sample from them for each cell. for now we are just using statistically independent spike
        trains for the inputs to each cell.
"""

morph_filename = 'EB2-late-bifurcation.swc'

mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'

num_exc_syns = 1600
num_inh_syns = 600

# in range(10), used for synapse locations, field_location, and offsets for synapses_seed and trial_seed
if len(sys.argv) > 1:
    gid = int(sys.argv[1])
else:
    gid = 0
# allows parallel computation of multiple trials for the same spines with the same peak_locs, but with different
# input spike trains and stochastic synapses for each trial
if len(sys.argv) > 2:
    trial_seed_offset = int(sys.argv[2])
else:
    trial_seed_offset = 0

field_peak_loc = np.arange(60., 120., 6.)[gid]

synapses_seed = gid + 100

trial_seed = trial_seed_offset + gid * 1000

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())+'-gid'+\
               str(gid)+'-e'+str(num_exc_syns)+'-i'+str(num_inh_syns)+'-peak_loc'+str(int(field_peak_loc))+'_spw'+\
               str(trial_seed_offset)


def run_trial(simiter, run_sim=True):
    """

    :param simiter: int
    :param global_phase_offset: float
    :param run_sim: bool
    """
    local_random.seed(simiter)
    if run_sim:
        with h5py.File(data_dir + rec_filename + '-working.hdf5', 'a') as f:
            f.create_group(str(simiter))
            f[str(simiter)].create_group('train')
            f[str(simiter)].create_group('inh_train')
    rate_maps = {}
    index = 0
    extended_t = np.arange(-200., 400., dt)
    for group in stim_exc_syns:
        rate_maps[group] = []
        for i, syn in enumerate(stim_exc_syns[group]):
            # the stochastic sequence used for each synapse is unique for each trial,
            # up to 10000 input spikes per spine
            if excitatory_stochastic:
                syn.randObj.seq(rand_exc_seq_locs[group][i] + int(simiter * 1e4))
            if group == 'CA3':
                peak_loc = peak_locs[group][i]
                peak_time = peak_loc / 187. * 200.
                start = int((200. - spw_rate_peak_delay + peak_time) / dt)
                stim_force = np.zeros_like(extended_t)
                stim_force[start:start+len(spw_rate_shape)] = (excitatory_peak_rate[group] -
                                                               excitatory_background_rate[group]) * spw_rate_shape
                before = np.array(stim_force[:len(spw_rate_shape)])
                after = np.array(stim_force[2*len(spw_rate_shape):])
                stim_force[len(spw_rate_shape):len(spw_rate_shape) + len(before)] += before
                stim_force[len(spw_rate_shape):len(spw_rate_shape)+len(after)] += after
                stim_force = np.array(stim_force[len(spw_rate_shape):len(spw_rate_shape)+len(spw_rate_shape)])
                rate_map = np.zeros_like(stim_t)
                rate_map[int(track_equilibrate/dt):int(track_equilibrate/dt)+len(stim_force)] = stim_force
                rate_map = np.multiply(rate_map, spw_filter)
            elif group == 'ECIII':
                rate_map = np.zeros_like(stim_t)
            rate_map += excitatory_background_rate[group]
            if run_sim:
                train = get_inhom_poisson_spike_times_by_thinning(rate_map, stim_t, dt=dt, generator=local_random)
                syn.source.play(h.Vector(np.add(train, equilibrate + track_equilibrate)))
                with h5py.File(data_dir + rec_filename + '-working.hdf5', 'a') as f:
                    f[str(simiter)]['train'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                            data=train)
                    f[str(simiter)]['train'][str(index)].attrs['group'] = group
                    f[str(simiter)]['train'][str(index)].attrs['index'] = syn.node.index
                    f[str(simiter)]['train'][str(index)].attrs['type'] = syn.node.parent.parent.type
                    f[str(simiter)]['train'][str(index)].attrs['peak_loc'] = peak_locs[group][i]
                    f[str(simiter)]['train'][str(index)].attrs['weight'] = weights[group][i]
            else:
                rate_maps[group].append(rate_map)
            index += 1
    index = 0
    for group in stim_inh_syns:
        rate_maps[group] = []
        if group in inhibitory_ripple_modulation_depth:
            inh_peak_rate = 2. * inhibitory_mean_rate[group] / (2. - inhibitory_ripple_modulation_depth[group])
            inhibitory_ripple_force = np.cos(2. * np.pi * spw_rate_shape_t / ripple_cycle_duration)
            inhibitory_ripple_force -= np.min(inhibitory_ripple_force)
            inhibitory_ripple_force /= np.max(inhibitory_ripple_force)
            inhibitory_ripple_force *= inhibitory_ripple_modulation_depth[group]
            inhibitory_ripple_force += 1. - inhibitory_ripple_modulation_depth[group]
        else:
            inh_peak_rate = inhibitory_mean_rate[group]
            inhibitory_ripple_force = np.ones_like(stim_t)
        stim_force = np.exp(inhibitory_spw_tuning_factor) - np.exp(inhibitory_spw_tuning_factor *
                                                                   np.cos(2.*np.pi*spw_rate_shape_t/spw_duration))
        stim_force /= np.max(stim_force)
        stim_force *= (inhibitory_spw_modulation[group] - 1.) * inh_peak_rate
        rate_map = np.zeros_like(stim_t)
        rate_map[int((track_equilibrate + spw_delay) / dt):int((track_equilibrate + spw_delay + spw_duration) / dt)] = \
            stim_force[:int(spw_duration / dt)]
        rate_map += inhibitory_mean_rate[group]
        rate_map[int((track_equilibrate + spw_delay) / dt):int((track_equilibrate + spw_delay + spw_duration) / dt)] = \
            np.multiply(
                rate_map[int((track_equilibrate + spw_delay) / dt):int((track_equilibrate + spw_delay + spw_duration) / dt)],
                inhibitory_ripple_force[:int(spw_duration / dt)])
        for syn in stim_inh_syns[group]:
            if run_sim:
                train = get_inhom_poisson_spike_times_by_thinning(np.array(rate_map), stim_t, dt=dt,
                                                                  generator=local_random)
                syn.source.play(h.Vector(np.add(train, equilibrate + track_equilibrate)))
                with h5py.File(data_dir + rec_filename + '-working.hdf5', 'a') as f:
                    f[str(simiter)]['inh_train'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                                data=train)
                    f[str(simiter)]['inh_train'][str(index)].attrs['group'] = group
                    f[str(simiter)]['inh_train'][str(index)].attrs['index'] = syn.node.index
                    f[str(simiter)]['inh_train'][str(index)].attrs['loc'] = syn.loc
                    f[str(simiter)]['inh_train'][str(index)].attrs['type'] = syn.node.type
            else:
                rate_maps[group].append(rate_map)
            index += 1
    if run_sim:
        sim.run(v_init)
        with h5py.File(data_dir + rec_filename + '-working.hdf5', 'a') as f:
            sim.export_to_file(f, simiter)
            if excitatory_stochastic:
                f[str(simiter)].create_group('successes')
                index = 0
                for group in stim_exc_syns:
                    for syn in stim_exc_syns[group]:
                        f[str(simiter)]['successes'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                                    data=np.subtract(syn.netcon(
                                                                        'AMPA_KIN').get_recordvec().to_python(),
                                                                                     equilibrate + track_equilibrate))
                        index += 1
            # save the spike output of the cell, removing the equilibration offset
            f[str(simiter)].create_dataset('output', compression='gzip', compression_opts=9,
                                           data=np.subtract(cell.spike_detector.get_recordvec().to_python(),
                                                            equilibrate + track_equilibrate))
    else:
        return rate_maps


NMDA_type = 'NMDA_KIN5'

dt = 0.02  # ms
equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # (ms)
input_field_width = 90.  # cm
track_length = 187.  # cm
default_run_vel = 30.  # cm/s

dx = dt * default_run_vel / 1000.
x = np.arange(0., track_length, dx)

extended_x = np.concatenate([x - track_length, x, x + track_length])

track_equilibrate = 2. * global_theta_cycle_duration

excitatory_peak_rate = {'CA3': 70., 'ECIII': 40.}
excitatory_background_rate = {'CA3': 2., 'ECIII': 2.}
excitatory_stochastic = 1
inhibitory_mean_rate = {'perisomatic': 10., 'axo-axonic': 10., 'apical dendritic': 10., 'distal apical dendritic': 10.,
                        'tuft feedforward': 10., 'tuft feedback': 10.}  # Hz
inhibitory_spw_modulation = {'perisomatic': 6., 'axo-axonic': 0.5, 'apical dendritic': 2.,
                             'distal apical dendritic': 0.5, 'tuft feedforward': 2., 'tuft feedback': 2.}
inhibitory_ripple_modulation_depth = {'perisomatic': 0.8}
inhibitory_spw_tuning_factor = 2.

spw_duration = 100.
spw_delay = 50.
spw_rate_rise = 5.
spw_rate_decay = 10.
spw_rate_shape_t = np.arange(0., 200., dt)
spw_rate_shape = np.exp(-spw_rate_shape_t/spw_rate_decay) - np.exp(-spw_rate_shape_t/spw_rate_rise)
spw_rate_shape /= np.max(spw_rate_shape)
spw_rate_peak_delay = dt * np.where(spw_rate_shape == np.max(spw_rate_shape))[0][0]

stim_duration = 300.
duration = equilibrate + track_equilibrate + stim_duration

stim_t = np.append(np.arange(-track_equilibrate, 0., dt), np.arange(0., stim_duration, dt))

spw_filter = np.zeros_like(stim_t)
spw_filter[int((track_equilibrate+spw_delay)/dt):int((track_equilibrate+spw_delay+spw_duration)/dt)] = \
    np.array(np.exp(5.) -
             np.exp(5. * np.cos(2.*np.pi*stim_t / 100.)))[int(track_equilibrate/dt):int((track_equilibrate+
                                                                                         spw_duration)/dt)]
spw_filter /= np.max(spw_filter)
ripple_freq = 150.  # Hz
ripple_cycle_duration = 1000. / ripple_freq  # ms

v_init = -67.

syn_types = ['AMPA_KIN', NMDA_type]

local_random = random.Random()

local_random.seed(trial_seed_offset)
global_phase_offset = local_random.uniform(-np.pi, np.pi)

# choose a subset of synapses to stimulate with inhomogeneous poisson rates
local_random.seed(synapses_seed)

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True, gid=gid)
cell.set_terminal_branch_na_gradient()
# cell.zero_na()
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
all_inh_syns = {sec_type: [] for sec_type in ['soma', 'ais', 'basal', 'trunk', 'apical', 'tuft']}
stim_exc_syns = {'CA3': [], 'ECIII': []}
stim_inh_syns = {'perisomatic': [], 'axo-axonic': [], 'apical dendritic': [], 'distal apical dendritic': [],
                 'tuft feedforward': [], 'tuft feedback': []}
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

sim_dt = 0.01
sim = QuickSim(duration, cvode=0, dt=sim_dt)
# sim = QuickSim(duration, cvode=0, dt=0.02)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['track_equilibrate'] = track_equilibrate
sim.parameters['global_theta_cycle_duration'] = global_theta_cycle_duration
sim.parameters['track_length'] = track_length
sim.parameters['duration'] = duration
sim.parameters['stim_dt'] = dt
sim.parameters['dt'] = sim_dt
sim.parameters['field_peak_loc'] = field_peak_loc
sim.parameters['gid'] = gid
sim.parameters['trial_seed'] = trial_seed_offset
sim.parameters['run_vel'] = default_run_vel
sim.parameters['spw_duration'] = spw_duration
sim.parameters['spw_delay'] = spw_delay
sim.append_rec(cell, cell.tree.root, description='soma', loc=0.)
sim.append_rec(cell, trunk_bifurcation[0], description='proximal_trunk', loc=1.)
sim.append_rec(cell, trunk, description='distal_trunk', loc=1.)
spike_output_vec = h.Vector()
cell.spike_detector.record(spike_output_vec)

# get the fraction of total spines contained in each sec_type
total_exc_syns = {sec_type: len(all_exc_syns[sec_type]) for sec_type in ['basal', 'trunk', 'apical', 'tuft']}
fraction_exc_syns = {sec_type: float(total_exc_syns[sec_type]) / float(np.sum(total_exc_syns.values())) for sec_type in
                 ['basal', 'trunk', 'apical', 'tuft']}

for sec_type in all_exc_syns:
    for i in local_random.sample(range(len(all_exc_syns[sec_type])), int(num_exc_syns*fraction_exc_syns[sec_type])):
        syn = all_exc_syns[sec_type][i]
        if sec_type == 'tuft':
            stim_exc_syns['ECIII'].append(syn)
        else:
            stim_exc_syns['CA3'].append(syn)

# get the fraction of inhibitory synapses contained in each sec_type
total_inh_syns = {sec_type: len(all_inh_syns[sec_type]) for sec_type in ['soma', 'ais', 'basal', 'trunk', 'apical',
                                                                         'tuft']}
fraction_inh_syns = {sec_type: float(total_inh_syns[sec_type]) / float(np.sum(total_inh_syns.values())) for sec_type in
                 ['soma', 'ais', 'basal', 'trunk', 'apical', 'tuft']}
num_inh_syns = min(num_inh_syns, int(np.sum(total_inh_syns.values())))

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
        elif syn.node.type == 'trunk':
            distance = cell.get_distance_to_node(cell.tree.root, syn.node, syn.loc)
            if distance <= 50.:
                group = 'perisomatic'
            elif distance <= 150.:
                group = local_random.choice(['apical dendritic', 'apical dendritic', 'distal apical dendritic'])
            else:
                group = local_random.choice(['apical dendritic', 'distal apical dendritic', 'distal apical dendritic'])
        elif syn.node.type == 'basal':
            distance = cell.get_distance_to_node(cell.tree.root, syn.node, syn.loc)
            group = 'perisomatic' if distance <= 50. and not cell.is_terminal(syn.node) else 'apical dendritic'
        elif syn.node.type == 'soma':
            group = 'perisomatic'
        elif syn.node.type == 'apical':
            distance = cell.get_distance_to_node(cell.tree.root, cell.get_dendrite_origin(syn.node), loc=1.)
            if distance <= 150.:
                group = local_random.choice(['apical dendritic', 'apical dendritic', 'distal apical dendritic'])
            else:
                group = local_random.choice(['apical dendritic', 'distal apical dendritic', 'distal apical dendritic'])
        elif syn.node.type == 'ais':
            group = 'axo-axonic'
        stim_inh_syns[group].append(syn)

rand_exc_seq_locs = {}

for group in stim_exc_syns:
    rand_exc_seq_locs[group] = []
    for syn in stim_exc_syns[group]:
        if excitatory_stochastic:
            success_vec = h.Vector()
            stim_successes.append(success_vec)
            syn.netcon('AMPA_KIN').record(success_vec)
            rand_exc_seq_locs[group].append(syn.randObj.seq())
        all_exc_syns[syn.node.parent.parent.type].remove(syn)

weights = {}

delta_peak_locs = {}
extended_peak_locs = {}
within_track = {}
peak_mod_weight = 2.5
tuning_amp = (peak_mod_weight - 1.)
start_loc = field_peak_loc - input_field_width * 1.2 / 2.
end_loc = field_peak_loc + input_field_width * 1.2 / 2.
for group in stim_exc_syns:
    if stim_exc_syns[group]:
        delta_peak_locs[group] = track_length / len(stim_exc_syns[group])
        extended_peak_locs[group] = np.arange(-track_length, 2. * track_length, delta_peak_locs[group])
        within_track[group] = np.where((extended_peak_locs[group] >= 0.) &
                                       (extended_peak_locs[group] < track_length))[0]
        peak_locs[group] = extended_peak_locs[group][within_track[group]]
        extended_weights = np.zeros_like(extended_peak_locs[group])
        within_field = np.where((extended_peak_locs[group] >= start_loc) &
                                (extended_peak_locs[group] <= end_loc))[0]
        extended_weights[within_field] = tuning_amp / 2. * np.cos(2. * np.pi / (input_field_width * 1.2) *
                                                                  (extended_peak_locs[group][within_field] -
                                                                   field_peak_loc)) + tuning_amp / 2.
        this_weights = extended_weights[within_track[group]]
        before_track_indexes = np.where(extended_peak_locs[group] < 0.)[0]
        before_track = np.array(extended_weights[before_track_indexes[-min(len(before_track_indexes),
                                                                           len(this_weights)):]])
        after_track_indexes = np.where(extended_peak_locs[group] > track_length)[0]
        after_track = np.array(extended_weights[after_track_indexes[:min(len(after_track_indexes),
                                                                         len(this_weights))]])
        this_weights[-len(before_track):] += before_track
        this_weights[:len(after_track)] += after_track
        weights[group] = 1. + this_weights

for group in stim_exc_syns:
    peak_locs[group] = list(peak_locs[group])
    weights[group] = list(weights[group])
    indexes = range(len(peak_locs[group]))
    local_random.shuffle(indexes)
    peak_locs[group] = map(peak_locs[group].__getitem__, indexes)
    weights[group] = map(weights[group].__getitem__, indexes)
    for i, syn in enumerate(stim_exc_syns[group]):
        syn.netcon('AMPA_KIN').weight[0] = weights[group][i]

# rate_maps = run_trial(trial_seed, run_sim=False)

run_trial(trial_seed, global_phase_offset)
if os.path.isfile(data_dir+rec_filename+'-working.hdf5'):
    os.rename(data_dir+rec_filename+'-working.hdf5', data_dir+rec_filename+'.hdf5')

