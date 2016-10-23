__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import random
import sys
"""
In this version of the simulation, phase precession of CA3 inputs is implemented using the method from Chadwick et al.,
Elife, 2015, which uses a circular gaussian with a phase sensitivity factor that effectively compresses the range of
phases within each theta cycle that each input is active, which will reduce jitter across within-cycle input sequences.

"""
morph_filename = 'EB2-late-bifurcation.swc'
#mech_filename = '020516 altered km2 rinp - ampa nmda_kin5'
mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'


if len(sys.argv) > 1:
    synapses_seed = int(sys.argv[1])
else:
    synapses_seed = 0
if len(sys.argv) > 2:
    num_exc_syns = int(sys.argv[2])
else:
    num_exc_syns = 3200
if len(sys.argv) > 3:
    num_inh_syns = int(sys.argv[3])
else:
    num_inh_syns = 600
# whether to modulate the firing rate of all inhibitory inputs (0 = no, 1 = out of field at track start, 2 = in field, 
# 3 = entire length of track)
if len(sys.argv) > 4:
    mod_inh = int(sys.argv[4])
else:
    mod_inh = 0
# the synaptic AMPAR conductances at in-field inputs are multiplied by a factor with this value at the peak of the
# field, and decays with cosine spatial modulation away from the field
if len(sys.argv) > 5:
    mod_weights = float(sys.argv[5])
else:
    mod_weights = 2.5
# allows parallel computation of multiple trials for the same spines with the same peak_locs, but with different
# input spike trains and stochastic synapses for each trial
if len(sys.argv) > 6:
    trial_seed = int(sys.argv[6])
else:
    trial_seed = 0

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())+'-seed'+\
               str(synapses_seed)+'-e'+str(num_exc_syns)+'-i'+str(num_inh_syns)+'-subtr_mod_inh'+str(mod_inh)+\
               '-no_precess_'+str(mod_weights)+'_'+str(trial_seed)


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
    with h5py.File(data_dir+rec_filename+'-working.hdf5', 'a') as f:
        f.create_group(str(simiter))
        f[str(simiter)].create_group('train')
        f[str(simiter)].create_group('inh_train')
        f[str(simiter)].attrs['phase_offset'] = global_phase_offset / 2. / np.pi * global_theta_cycle_duration
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
            train = get_inhom_poisson_spike_times(stim_force, stim_t, dt=stim_dt, generator=local_random)
            syn.source.play(h.Vector(np.add(train, equilibrate + track_equilibrate)))
            with h5py.File(data_dir+rec_filename+'-working.hdf5', 'a') as f:
                f[str(simiter)]['train'].create_dataset(str(index), compression='gzip', compression_opts=9, data=train)
                f[str(simiter)]['train'][str(index)].attrs['group'] = group
                f[str(simiter)]['train'][str(index)].attrs['index'] = syn.node.index
                f[str(simiter)]['train'][str(index)].attrs['type'] = syn.node.parent.parent.type
                f[str(simiter)]['train'][str(index)].attrs['peak_loc'] = peak_locs[group][i]
            index += 1
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
            train = get_inhom_poisson_spike_times(stim_force, stim_t, dt=stim_dt, generator=local_random)
            syn.source.play(h.Vector(np.add(train, equilibrate + track_equilibrate)))
            with h5py.File(data_dir+rec_filename+'-working.hdf5', 'a') as f:
                f[str(simiter)]['inh_train'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                            data=train)
                f[str(simiter)]['inh_train'][str(index)].attrs['group'] = group
                f[str(simiter)]['inh_train'][str(index)].attrs['index'] = syn.node.index
                f[str(simiter)]['inh_train'][str(index)].attrs['loc'] = syn.loc
                f[str(simiter)]['inh_train'][str(index)].attrs['type'] = syn.node.type
            index += 1
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'-working.hdf5', 'a') as f:
        sim.export_to_file(f, simiter)
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


NMDA_type = 'NMDA_KIN5'

equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # (ms)
input_field_width = 20  # (theta cycles per 6 standard deviations)
input_field_duration = input_field_width * global_theta_cycle_duration
track_length = 2.5  # field widths
track_duration = track_length * input_field_duration
track_equilibrate = 2. * global_theta_cycle_duration
duration = equilibrate + track_equilibrate + track_duration  # input_field_duration
excitatory_peak_rate = {'CA3': 40., 'ECIII': 40.}
excitatory_theta_modulation_depth = {'CA3': 0.7, 'ECIII': 0.7}
# From Chadwick et al., ELife 2015
excitatory_theta_phase_tuning_factor = {'CA3': 0.8, 'ECIII': 0.8}
excitatory_precession_range = {}  # (ms, degrees)
# excitatory_precession_range['CA3'] = [(-input_field_duration*0.5, 180.), (-input_field_duration*0.35, 180.),
#                                      (input_field_duration*0.35, -180.), (input_field_duration*0.5, -180.)]
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

stim_dt = 0.02
dt = 0.02
v_init = -67.

syn_types = ['AMPA_KIN', NMDA_type]

local_random = random.Random()

# choose a subset of synapses to stimulate with inhomogeneous poisson rates
local_random.seed(synapses_seed)

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.set_terminal_branch_na_gradient()
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

sim = QuickSim(duration, cvode=0, dt=0.01)
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
        #peak_loc = local_random.uniform(-0.75 * input_field_duration, (0.75 + track_length) * input_field_duration)
        #peak_locs.append(peak_loc)
        if excitatory_stochastic:
            success_vec = h.Vector()
            stim_successes.append(success_vec)
            syn.netcon('AMPA_KIN').record(success_vec)
            rand_exc_seq_locs[group].append(syn.randObj.seq())
        # if syn.node.parent.parent not in [rec['node'] for rec in sim.rec_list]:
        #    sim.append_rec(cell, syn.node.parent.parent)
        # sim.append_rec(cell, syn.node, object=syn.target('AMPA_KIN'), param='_ref_i', description='i_AMPA')
        # sim.append_rec(cell, syn.node, object=syn.target(NMDA_type), param='_ref_i', description='i_NMDA')
        # remove this synapse from the pool, so that additional "modulated" inputs
        # can be selected from those that remain
        all_exc_syns[syn.node.parent.parent.type].remove(syn)

# rand_inh_seq_locs = [] will need this when inhibitory synapses become stochastic
# stim_inh_successes = [] will need this when inhibitory synapses become stochastic

# modulate the weights of inputs with peak_locs along this stretch of the track
modulated_field_center = track_duration * 0.6
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
    indexes = range(len(peak_locs[group]))
    local_random.shuffle(indexes)
    peak_locs[group] = map(peak_locs[group].__getitem__, indexes)
    cos_mod_weight[group] = map(cos_mod_weight[group].__getitem__, indexes)
    for i, syn in enumerate(stim_exc_syns[group]):
        syn.netcon('AMPA_KIN').weight[0] = cos_mod_weight[group][i]

run_trial(trial_seed)
if os.path.isfile(data_dir+rec_filename+'-working.hdf5'):
    os.rename(data_dir+rec_filename+'-working.hdf5', data_dir+rec_filename+'.hdf5')
