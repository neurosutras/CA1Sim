__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import random
import sys
"""

"""
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '020516 altered km2 rinp - ampa nmda_kin5'


if len(sys.argv) > 1:
    synapses_seed = int(sys.argv[1])
else:
    synapses_seed = 1
if len(sys.argv) > 2:
    num_exc_syns = int(sys.argv[2])
else:
    num_exc_syns = 3000
if len(sys.argv) > 3:
    num_inh_syns = int(sys.argv[3])
else:
    num_inh_syns = 500
# whether to modulate the peak rate of all inhibitory inputs (0 = no, 1 = out of field at track start, 2 = in field)
# input_field_width)
if len(sys.argv) > 4:
    mod_inh = int(sys.argv[4])
else:
    mod_inh = 0
# the synaptic AMPAR conductances at in-field inputs will be multiplied by a factor that follows an asymmetric kernel
# with this value at the peak of the field
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
               str(synapses_seed)+'-e'+str(num_exc_syns)+'-i'+str(num_inh_syns)+'-mod_inh'+str(mod_inh)+\
               '-weights_'+str(mod_weights)+'_'+str(trial_seed)


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
        elif mod_inh == 2:
            mod_inh_start = int((track_equilibrate + modulated_field_center - 0.3 * input_field_duration) / dt)
        sim.parameters['mod_inh_start'] = stim_t[mod_inh_start]
        mod_inh_stop = mod_inh_start + int(inhibitory_manipulation_duration * input_field_duration / dt)
        sim.parameters['mod_inh_stop'] = stim_t[mod_inh_stop]
    index = 0
    for group in stim_exc_syns:
        excitatory_theta_amp = excitatory_theta_modulation_depth[group] / 2.
        excitatory_theta_offset = 1. - excitatory_theta_amp
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
            train = get_inhom_poisson_spike_times_by_thinning(stim_force, stim_t, dt=stim_dt, generator=local_random)
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
        inhibitory_theta_amp = inhibitory_peak_rate[group] * inhibitory_theta_modulation_depth[group] / 2.
        inhibitory_theta_offset = inhibitory_peak_rate[group] - inhibitory_theta_amp
        inhibitory_phase_offset = inhibitory_theta_phase_offset[group]
        for syn in stim_inh_syns[group]:
            inhibitory_theta_force = inhibitory_theta_offset + inhibitory_theta_amp * np.cos(2. * np.pi /
                                                global_theta_cycle_duration * stim_t - global_phase_offset -
                                                inhibitory_phase_offset)
            if mod_inh > 0 and group in inhibitory_manipulation_fraction and syn in manipulated_inh_syns[group]:
                inhibitory_theta_force[mod_inh_start:mod_inh_stop] = 0.
            train = get_inhom_poisson_spike_times_by_thinning(inhibitory_theta_force, stim_t, dt=stim_dt,
                                                  generator=local_random)
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
excitatory_theta_modulation_depth = {'CA3': 0.75, 'ECIII': 0.7}
theta_compression_factor = 1. - unit_theta_cycle_duration / global_theta_cycle_duration
excitatory_theta_phase_offset = {}
excitatory_theta_phase_offset['CA3'] = 165. / 360. * 2. * np.pi  # radians
excitatory_theta_phase_offset['ECIII'] = 0. / 360. * 2. * np.pi  # radians
excitatory_stochastic = 1
inhibitory_manipulation_fraction = {'perisomatic': 0.35, 'axo-axonic': 0.35, 'apical dendritic': 0.35,
                                    'tuft feedback': 0.35}
inhibitory_manipulation_duration = 0.6  # Ratio of input_field_duration
inhibitory_peak_rate = {'perisomatic': 40., 'axo-axonic': 40., 'apical dendritic': 40., 'distal apical dendritic': 40.,
                        'tuft feedforward': 40., 'tuft feedback': 40.}
inhibitory_theta_modulation_depth = {'perisomatic': 0.5, 'axo-axonic': 0.5, 'apical dendritic': 0.5,
                                     'distal apical dendritic': 0.5, 'tuft feedforward': 0.5, 'tuft feedback': 0.5}
inhibitory_theta_phase_offset = {}
inhibitory_theta_phase_offset['perisomatic'] = 145. / 360. * 2. * np.pi  # Like PV+ Basket
inhibitory_theta_phase_offset['axo-axonic'] = 70. / 360. * 2. * np.pi  # Vargas et al., ELife, 2014
inhibitory_theta_phase_offset['apical dendritic'] = 210. / 360. * 2. * np.pi  # Like PYR-layer Bistratified
inhibitory_theta_phase_offset['distal apical dendritic'] = 180. / 360. * 2. * np.pi  # Like SR/SLM Border Cells
inhibitory_theta_phase_offset['tuft feedforward'] = 340. / 360. * 2. * np.pi  # Like Neurogliaform
inhibitory_theta_phase_offset['tuft feedback'] = 210. / 360. * 2. * np.pi  # Like SST+ O-LM

stim_dt = 0.02
dt = 0.02
v_init = -67.

syn_types = ['AMPA_KIN', NMDA_type]

local_random = random.Random()

# choose a subset of synapses to stimulate with inhomogeneous poisson rates
local_random.seed(synapses_seed)

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.set_terminal_branch_nas_gradient()
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
sim.append_rec(cell, cell.tree.root, description='soma', loc=0.5)
sim.append_rec(cell, trunk, description='distal_trunk', loc=0.)
sim.append_rec(cell, trunk_bifurcation[0], description='proximal_trunk', loc=1.)
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
                group = 'apical dendritic'
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
                group = 'apical dendritic'
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
trough_mod_weight = 0.8
tuning_amp1 = (peak_mod_weight - trough_mod_weight) / 2.
tuning_offset1 = tuning_amp1 + trough_mod_weight
tuning_amp2 = (peak_mod_weight - 1.) / 2.
tuning_offset2 = tuning_amp2 + 1.

pre_weights = {}
post_weights = {}
for group in stim_exc_syns:
    cos_mod_weight1 = tuning_amp1 * np.cos(2. * np.pi / (input_field_duration * 1.2) * (peak_locs[group] -
                                                                                       modulated_field_center)) + \
                      tuning_offset1
    cos_mod_weight2 = tuning_amp2 * np.cos(2. * np.pi / (input_field_duration * 1.2) * (peak_locs[group] -
                                                                                        modulated_field_center)) + \
                      tuning_offset2
    left = np.where((cos_mod_weight1 < 1.) & (peak_locs[group] < modulated_field_center -
                    input_field_duration * 1.2 / 2.) & (peak_locs[group] > modulated_field_center -
                    input_field_duration))[0][0]
    right = np.where(peak_locs[group] > modulated_field_center + input_field_duration * 1.2 / 2.)[0][0]
    center = np.where(peak_locs[group] > modulated_field_center)[0][0]
    cos_mod_weight[group] = np.array(cos_mod_weight1)
    cos_mod_weight[group][:left] = 1.
    cos_mod_weight[group][center:right] = cos_mod_weight2[center:right]
    cos_mod_weight[group][right:] = 1.
    peak_locs[group] = list(peak_locs[group])
    cos_mod_weight[group] = list(cos_mod_weight[group])
    indexes = range(len(peak_locs[group]))
    local_random.shuffle(indexes)
    peak_locs[group] = map(peak_locs[group].__getitem__, indexes)
    cos_mod_weight[group] = map(cos_mod_weight[group].__getitem__, indexes)
    pre_weights[group] = []
    post_weights[group] = []
    for i, syn in enumerate(stim_exc_syns[group]):
        pre_weight = syn.target('AMPA_KIN').gmax * 1000. * 0.3252  # scale by peak AMPAR PO and convert to nS
        pre_weights[group].append(pre_weight)
        syn.netcon('AMPA_KIN').weight[0] = cos_mod_weight[group][i]
        post_weights[group].append(pre_weight * cos_mod_weight[group][i])

manipulated_inh_syns = {}
for group in inhibitory_manipulation_fraction:
    num_syns = int(len(stim_inh_syns[group]) * inhibitory_manipulation_fraction[group])
    manipulated_inh_syns[group] = local_random.sample(stim_inh_syns[group], num_syns)
"""
run_trial(trial_seed)
if os.path.isfile(data_dir+rec_filename+'-working.hdf5'):
    os.rename(data_dir+rec_filename+'-working.hdf5', data_dir+rec_filename+'.hdf5')
"""
"""
svg_title = '030916 - cell68'
collapsed_peak_locs = []
collapsed_pre_weights = []
collapsed_post_weights = []
for group in peak_locs:
    collapsed_peak_locs.extend(peak_locs[group])
    collapsed_pre_weights.extend(pre_weights[group])
    collapsed_post_weights.extend(post_weights[group])
pre_bins, pre_density, pre_weight_mean = sliding_window(collapsed_peak_locs, collapsed_pre_weights)
post_bins, post_density, post_weight_mean = sliding_window(collapsed_peak_locs, collapsed_post_weights)
fig, axes = plt.subplots(1)
axes.scatter(pre_bins, pre_weight_mean, color='k', label='Postsynaptic Input Weight Distribution')
axes.set_xlim(0., 7500.)
axes.set_ylim(0., 2.)
axes.set_xlabel('Time (ms)', fontsize=20)
axes.set_ylabel('Peak Synaptic AMPAR\nConductances (nS)', fontsize=20)
axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=20)
clean_axes(axes)
plt.savefig(data_dir+svg_title+' - premod weights.svg', format='svg')
plt.show()
plt.close()
fig, axes = plt.subplots(1)
axes.scatter(post_bins, post_weight_mean, color='k', label='Postsynaptic Input Weight Distribution')
axes.set_xlim(0., 7500.)
axes.set_ylim(0., 2.)
axes.set_xlabel('Time (ms)', fontsize=20)
axes.set_ylabel('Peak Synaptic AMPAR\nConductances (nS)', fontsize=20)
axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=20)
clean_axes(axes)
plt.savefig(data_dir+svg_title+' - postmod weights.svg', format='svg')
plt.show()
plt.close()
fig, axes = plt.subplots(1)
axes.scatter(post_bins, post_density, color='r', label='Presynaptic Input Density Distribution')
axes.set_xlim(0., 7500.)
axes.set_ylim(0., 600.)
axes.set_xlabel('Time (ms)', fontsize=20)
axes.set_ylabel('Input Density (/s)', fontsize=20)
axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=20)
clean_axes(axes)
plt.savefig(data_dir+svg_title+' - input density.svg', format='svg')
plt.show()
plt.close()
"""