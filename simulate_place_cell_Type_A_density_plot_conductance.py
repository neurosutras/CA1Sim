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
    num_exc_syns = 800
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
# the inhibitory conductances in-field are multiplied by a factor with this value at the peak of the
# field, and decays with cosine spatial modulation away from the field
if len(sys.argv) > 5:
    shape_inh = float(sys.argv[5])
else:
    shape_inh = 1.
# allows parallel computation of multiple trials for the same spines with the same peak_locs, but with different
# input spike trains and stochastic synapses for each trial
if len(sys.argv) > 6:
    trial_seed = int(sys.argv[6])
else:
    trial_seed = 0


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
            this_range_piece = np.ones(del_t / dt) * phase0
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
    global_phase_offset = 0.  # local_random.uniform(-np.pi, np.pi)
    if mod_inh > 0:
        if mod_inh == 1:
            mod_inh_start = int(track_equilibrate / dt)
            mod_inh_stop = mod_inh_start + int(inhibitory_manipulation_duration * input_field_duration / dt)
        elif mod_inh == 2:
            mod_inh_start = int((track_equilibrate + modulated_field_center - 0.3 * input_field_duration) / dt)
            mod_inh_stop = mod_inh_start + int(inhibitory_manipulation_duration * input_field_duration / dt)
        elif mod_inh == 3:
            mod_inh_start = 0
            mod_inh_stop = len(stim_t)  # - 1
    index = 0
    for group in stim_exc_syns:
        exc_stim_forces[group] = []
        for i, syn in enumerate(stim_exc_syns[group]):
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
            exc_stim_forces[group].append(stim_force)
            index += 1
    index = 0
    for group in stim_inh_syns:
        inh_stim_forces[group] = []
        for syn in stim_inh_syns[group]:
            inhibitory_theta_force = np.exp(inhibitory_theta_phase_tuning_factor[group] *
                                 np.cos(inhibitory_theta_phase_offset[group] - 2. * np.pi * stim_t /
                                        global_theta_cycle_duration + global_phase_offset))
            inhibitory_theta_force -= np.min(inhibitory_theta_force)
            inhibitory_theta_force /= np.max(inhibitory_theta_force)
            inhibitory_theta_force *= inhibitory_theta_modulation_depth[group]
            inhibitory_theta_force += 1. - inhibitory_theta_modulation_depth[group]
            inhibitory_theta_force *= inhibitory_peak_rate[group]
            stim_force = cos_mod_inh * np.mean(inhibitory_theta_force)  # mean of theta modulation
            if mod_inh > 0 and group in inhibitory_manipulation_fraction and syn in manipulated_inh_syns[group]:
                if not mod_inh == 3:
                    stim_force[mod_inh_start:mod_inh_stop] = 0.
                    inh_stim_forces[group].append(stim_force)
            else:
                inh_stim_forces[group].append(stim_force)
            index += 1


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
excitatory_precession_range = {}
excitatory_precession_range['CA3'] = [(-input_field_duration*0.5, 180.), (-input_field_duration*0.35, 180.),
                                      (input_field_duration*0.35, -180.), (input_field_duration*0.5, -180.)]  # (ms, degrees)
excitatory_theta_phase_offset = {}
excitatory_theta_phase_offset['CA3'] = 165. / 360. * 2. * np.pi  # radians
excitatory_theta_phase_offset['ECIII'] = 0. / 360. * 2. * np.pi  # radians
excitatory_stochastic = 1
inhibitory_manipulation_fraction = {'perisomatic': 0.325, 'axo-axonic': 0.325, 'apical dendritic': 0.325,
                                    'distal apical dendritic': 0.325, 'tuft feedback': 0.325}
inhibitory_manipulation_duration = 0.6  # Ratio of input_field_duration
inhibitory_peak_rate = {'perisomatic': 40., 'axo-axonic': 40., 'apical dendritic': 40., 'distal apical dendritic': 40.,
                        'tuft feedforward': 40., 'tuft feedback': 40.}
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

stim_dt = 0.1
dt = 0.1
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

stim_t = np.arange(-track_equilibrate, track_duration + track_equilibrate, dt)

gauss_sigma = global_theta_cycle_duration * input_field_width / 3. / np.sqrt(2.)  # contains 99.7% gaussian area

for group in stim_exc_syns:
    if stim_exc_syns[group]:
        peak_locs[group] = np.arange(-0.75 * input_field_duration, (0.75 + track_length) * input_field_duration,
                          (1.5 + track_length) * input_field_duration / int(len(stim_exc_syns[group])))
        peak_locs[group] = peak_locs[group][:len(stim_exc_syns[group])]

for group in stim_exc_syns:
    for syn in stim_exc_syns[group]:
        all_exc_syns[syn.node.parent.parent.type].remove(syn)

# modulate the number and density of inputs with peak_locs along this stretch of the track
modulated_field_center = track_duration * 0.6
mod_density = 1.98
compression_factor = 0.15
peak_loc_choices = {}
peak_loc_probabilities = {}
pre_modulation_num_exc_syns = {}
total_modulated_num_exc_syns = 0
np_local_random = np.random.RandomState()

gauss_mod_probability = np.exp(-((stim_t - modulated_field_center) / (gauss_sigma * compression_factor)) ** 2.)
indexes = np.where(gauss_mod_probability > 0.01)[0]
start = stim_t[indexes[0]]
end = stim_t[indexes[-1]]

for group in stim_exc_syns:
    pre_modulation_num_exc_syns[group] = len(stim_exc_syns[group])
    peak_loc_t = np.arange(-0.75 * input_field_duration, (0.75 + track_length) * input_field_duration,
                          (1.5 + track_length) * input_field_duration / int(len(stim_exc_syns[group])) / 4.)
    in_field_indexes = np.where((peak_loc_t >= start) & (peak_loc_t <= end))[0]
    in_field_peak_locs = peak_loc_t[in_field_indexes]
    peak_loc_probabilities = np.exp(-((in_field_peak_locs - modulated_field_center) /
                                      (gauss_sigma * compression_factor)) ** 2.)
    peak_loc_probabilities /= np.sum(peak_loc_probabilities)  # sum of probabilities must equal 1
    baseline_num_exc_syns = len(np.where((peak_locs[group] >= modulated_field_center - input_field_duration / 2.) &
                                np.array(peak_locs[group] <= modulated_field_center + input_field_duration / 2.))[0])
    this_modulated_num_exc_syns = int(baseline_num_exc_syns * (mod_density - 1.))
    peak_loc_choices[group] = np_local_random.choice(in_field_peak_locs, this_modulated_num_exc_syns,
                                                     p=peak_loc_probabilities)
    total_modulated_num_exc_syns += len(peak_loc_choices[group])

modulated_num_exc_syns = {group: 0 for group in stim_exc_syns}
for sec_type in all_exc_syns:
    group = 'ECIII' if sec_type == 'tuft' else 'CA3'
    for i in local_random.sample(range(len(all_exc_syns[sec_type])),
                                 min(int(total_modulated_num_exc_syns*fraction_exc_syns[sec_type]),
                                     len(all_exc_syns[sec_type]))):
        syn = all_exc_syns[sec_type][i]
        if modulated_num_exc_syns[group] < len(peak_loc_choices[group]):
            stim_exc_syns[group].append(syn)
            modulated_num_exc_syns[group] += 1

for group in stim_exc_syns:
    peak_locs[group] = np.append(peak_locs[group], peak_loc_choices[group][:modulated_num_exc_syns[group]])
    for syn in stim_exc_syns[group][pre_modulation_num_exc_syns[group]:]:
        all_exc_syns[syn.node.parent.parent.type].remove(syn)

# modulate the weights of inputs with peak_locs along this stretch of the track
cos_mod_weight = {}
peak_mod_weight = 2.5
tuning_amp = (peak_mod_weight - 1.) / 2.
tuning_offset = tuning_amp + 1.

for group in stim_exc_syns:
    cos_mod_weight[group] = tuning_amp * np.cos(2. * np.pi / (input_field_duration * 1.2) * (peak_locs[group] -
                                                                        modulated_field_center)) + tuning_offset
    before = np.where(peak_locs[group] < modulated_field_center - input_field_duration * 1.2 / 2.)[0]
    after = np.where(peak_locs[group] > modulated_field_center + input_field_duration * 1.2 / 2.)[0]
    cos_mod_weight[group][before] = 1.
    cos_mod_weight[group][after] = 1.

manipulated_inh_syns = {}
for group in inhibitory_manipulation_fraction:
    num_syns = int(len(stim_inh_syns[group]) * inhibitory_manipulation_fraction[group])
    manipulated_inh_syns[group] = local_random.sample(stim_inh_syns[group], num_syns)

inh_tuning_amp = (shape_inh - 1.) / 2.
inh_tuning_offset = inh_tuning_amp + 1.
left = np.where(stim_t >= modulated_field_center - input_field_duration * 1.2 / 2.)[0][0]
right = np.where(stim_t > modulated_field_center + input_field_duration * 1.2 / 2.)[0][0]
cos_mod_inh = inh_tuning_amp * np.cos(2. * np.pi / (input_field_duration * 1.2) * (stim_t - modulated_field_center)) \
              + inh_tuning_offset
cos_mod_inh[:left] = 1.
cos_mod_inh[right:] = 1.

exc_stim_forces = {}
inh_stim_forces = {}
run_trial(trial_seed)
filter_t = np.arange(0., 150., stim_dt)
t = np.arange(0., track_duration, stim_dt)
start = int(track_equilibrate/stim_dt)
ampa_forces = {}
ampa_forces_sum = {}
gaba_forces = {}
gaba_forces_sum = {}

ampa_filter = np.exp(-filter_t/3.) - np.exp(-filter_t/.1)
ampa_filter /= np.sum(ampa_filter)
ampa_conversion_factor = 0.000002
gaba_filter = np.exp(-filter_t/9.) - np.exp(-filter_t/.2)
gaba_filter /= np.sum(gaba_filter)
gaba_conversion_factor = 0.000002

for group in exc_stim_forces:
    ampa_forces[group] = []
    for i, stim_force in enumerate(exc_stim_forces[group]):
        this_force = stim_force * cos_mod_weight[group][i] * ampa_conversion_factor
        this_force = np.convolve(this_force, ampa_filter)[:len(stim_t)]
        ampa_forces[group].append(this_force)

for group in inh_stim_forces:
    gaba_forces[group] = []
    for i, stim_force in enumerate(inh_stim_forces[group]):
        this_force = stim_force * gaba_conversion_factor
        this_force = np.convolve(this_force, gaba_filter)[:len(stim_t)]
        gaba_forces[group].append(this_force)

ampa_forces_sum['all'] = np.zeros_like(stim_t)
for group in ampa_forces:
    this_sum = np.sum(ampa_forces[group], axis=0)
    filtered = low_pass_filter(this_sum[start:-start], 2., track_duration, stim_dt)
    ampa_forces_sum[group] = filtered
    ampa_forces_sum['all'] = np.add(ampa_forces_sum['all'], this_sum)
filtered = low_pass_filter(ampa_forces_sum['all'][start:-start], 2., track_duration, stim_dt)
ampa_forces_sum['all'] = filtered

gaba_forces_sum['all'] = np.zeros_like(stim_t)
for group in gaba_forces:
    gaba_forces_sum[group] = np.sum(gaba_forces[group], axis=0)
    gaba_forces_sum['all'] = np.add(gaba_forces_sum['all'], gaba_forces_sum[group])
filtered = low_pass_filter(gaba_forces_sum['all'][start:-start], 2., track_duration, stim_dt)
gaba_forces_sum['all'] = filtered

plt.plot(t, ampa_forces_sum['all'])
plt.show()
plt.close()

plt.plot(t, gaba_forces_sum['all'])
plt.show()
plt.close()