__author__ = 'milsteina'
from specify_cells import *
from scipy import signal
import random
import sys
"""

"""
morph_filename = 'EB2-late-bifurcation.swc'
#mech_filename = '052915 pas_exp_scale kdr ka_scale ih_sig_scale ampar_exp_scale nmda - EB2'
#mech_filename = '071515 rebalanced nax kap kdr pas h - EB2 - spines'
#mech_filename = '071715 rebalanced na_kap_kdr_pas_h - EB2 - spines'
#mech_filename = '072515 optimized basal ka_scale dend_sh_ar_nas - EB2'
mech_filename = '080615 rebalanced na_ka ampa nmda - EB2'

if len(sys.argv) > 1:
    seed = sys.argv[1]
else:
    seed = 1

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())+'-seed'+str(seed)


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
        stim_trains = []
        global_phase_offset = local_random.uniform(0., global_theta_cycle_duration)
        for i, syn in enumerate(stim_syns):
            if syn.node.parent.parent.type == 'tuft':
                theta_force = 1. + np.sin(2. * np.pi / global_theta_cycle_duration * (stim_t - global_phase_offset +
                                                                                      tuft_phase_offset))
            else:
                unit_phase_offset = peak_locs[i] * theta_compression_factor
                theta_force = 1. + np.sin(2. * np.pi / unit_theta_cycle_duration * (stim_t - global_phase_offset -
                                                                              unit_phase_offset))
            start = int(((0.75 + track_length) * input_field_duration - peak_locs[i]) / dt)
            buffer = int((0.75 * input_field_duration - track_equilibrate) / dt)
            start += buffer
            end = start + len(stim_t)
            stim_force = gauss_force[start:end]
            stim_force = np.multiply(stim_force, theta_force)
            if simiter == 1:
                stim_forces.append(stim_force)
            train = get_inhom_poisson_spike_times(stim_force, stim_t, generator=local_random)
            stim_trains.append(train)
            syn.source.play(h.Vector(np.add(train, equilibrate + track_equilibrate)))
        stim_iterations.append(stim_trains)
        sim.run(v_init)
        with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
            sim.export_to_file(f, simiter)
            f[str(simiter)].create_group('train')
            f[str(simiter)].create_group('successes')
            f[str(simiter)]['train'].attrs['phase_offset'] = global_phase_offset
            for index, train in enumerate(stim_trains):
                f[str(simiter)]['train'].create_dataset(str(index), compression='gzip', compression_opts=9, data=train)
                f[str(simiter)]['train'][str(index)].attrs['index'] = stim_syns[index].node.index
                f[str(simiter)]['train'][str(index)].attrs['type'] = stim_syns[index].node.parent.parent.type
                f[str(simiter)]['successes'].create_dataset(str(index), compression='gzip', compression_opts=9,
                        data=np.subtract(stim_syns[index].netcon('AMPA_KIN').get_recordvec().to_python(),
                                         equilibrate + track_equilibrate))
                f[str(simiter)]['train'][str(index)].attrs['peak_loc'] = peak_locs[index]
            # save the spike output of the cell, removing the equilibration offset
            f[str(simiter)].create_dataset('output', compression='gzip', compression_opts=9,
                                        data=np.subtract(cell.spike_detector.get_recordvec().to_python(),
                                                         equilibrate + track_equilibrate))
    trials += n


NMDA_type = 'NMDA_KIN2'

equilibrate = 250.  # time to steady-state
global_theta_cycle_duration = 150.  # (ms)
input_field_width = 10  # (theta cycles per 7 standard deviations)
# Geissler...Buzsaki, PNAS 2010
unit_theta_cycle_duration = global_theta_cycle_duration * input_field_width / (input_field_width + 1.)
input_field_duration = input_field_width * global_theta_cycle_duration
track_length = 3  # field widths
track_duration = track_length * input_field_duration
track_equilibrate = 2. * global_theta_cycle_duration
duration = equilibrate + track_equilibrate + track_duration
gaussian_modulation_strength = 25.
theta_compression_factor = unit_theta_cycle_duration / input_field_duration
tuft_phase_offset = 45. /360. * global_theta_cycle_duration

stim_dt = 0.1
dt = 0.02
v_init = -67.

num_syns = 2000  # evenly distributed across sec_types
syn_types = ['AMPA_KIN', NMDA_type]

local_random = random.Random()

# choose a subset of synapses to stimulate with inhomogeneous poisson rates
# cell1 and cell2
local_random.seed(seed)

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
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

all_syns = {sec_type: [] for sec_type in ['basal', 'trunk', 'apical', 'tuft']}
stim_syns = []
stim_forces = []
stim_successes = []
peak_locs = []

# place synapses in every spine
for sec_type in all_syns:
    for node in cell.get_nodes_of_subtype(sec_type):
        for spine in node.spines:
            syn = Synapse(cell, spine, syn_types, stochastic=1)
            all_syns[sec_type].append(syn)
cell.init_synaptic_mechanisms()

sim = QuickSim(duration)
sim.parameters['equilibrate'] = equilibrate
sim.parameters['track_equilibrate'] = track_equilibrate
sim.parameters['input_field_duration'] = input_field_duration
sim.parameters['track_length'] = track_length
sim.parameters['duration'] = duration
sim.parameters['stim_dt'] = stim_dt
sim.append_rec(cell, cell.tree.root, description='soma', loc=0.5)
sim.append_rec(cell, trunk, description='trunk', loc=0.)
spike_output_vec = h.Vector()
cell.spike_detector.record(spike_output_vec)

# get the fraction of total spines contained in each dendritic type
total_syns = {sec_type: len(all_syns[sec_type]) for sec_type in ['basal', 'trunk', 'apical', 'tuft']}
fraction_syns = {sec_type: float(total_syns[sec_type]) / float(np.sum(total_syns.values())) for sec_type in
                 ['basal', 'trunk', 'apical', 'tuft']}

for sec_type in all_syns:
    for i in local_random.sample(range(len(all_syns[sec_type])), int(num_syns*fraction_syns[sec_type])):
        syn = all_syns[sec_type][i]
        stim_syns.append(syn)

stim_t = np.arange(-track_equilibrate, track_duration, dt)

gauss_sigma = global_theta_cycle_duration * input_field_width / 6.  # contains 99.7% gaussian area
gauss_force = gaussian_modulation_strength * signal.gaussian(int((2 * (track_length + 1.5) *
                                                    input_field_duration) / dt), gauss_sigma / dt)
for syn in stim_syns:
    peak_loc = local_random.uniform(-0.75 * input_field_duration, (0.75 + track_length) * input_field_duration)
    peak_locs.append(peak_loc)
    success_vec = h.Vector()
    stim_successes.append(success_vec)
    syn.netcon('AMPA_KIN').record(success_vec)
    if syn.node.parent.parent not in [rec['node'] for rec in sim.rec_list]:
        sim.append_rec(cell, syn.node.parent.parent)

stim_iterations = []
trials = 0
print 'Getting started with', num_syns, 'inputs, seed:', seed
run_n_trials(10)

"""
stim_forces = []
theta_forces = []
peak_locs = []
global_phase_offset = 0.  # local_random.uniform(0., global_theta_cycle_duration)
for i in range(500):
    peak_loc = local_random.uniform(-0.75 * input_field_duration, (0.75 + track_length) * input_field_duration)
    peak_locs.append(peak_loc)
    unit_phase_offset = peak_locs[i] * theta_compression_factor
    theta_force = 1. + np.sin(2. * np.pi / unit_theta_cycle_duration * (stim_t - global_phase_offset -
                                                                  unit_phase_offset))
    theta_forces.append(theta_force)
    start = int(((0.75 + track_length) * input_field_duration - peak_locs[i])/dt)
    buffer = int((0.75 * input_field_duration - track_equilibrate) / dt)
    start += buffer
    end = start + len(stim_t)
    stim_force = gauss_force[start:end]
    stim_force = np.multiply(stim_force, theta_force)
    stim_forces.append(stim_force)
    plt.plot(stim_t, stim_force)


for i, syn in enumerate(stim_syns):
    if syn.node.parent.parent.type == 'tuft':
        theta_force = 1. + np.sin(2. * np.pi / global_theta_cycle_duration * (stim_t - global_phase_offset +
                                                                              tuft_phase_offset))
    else:
        unit_phase_offset = peak_locs[i] * theta_compression_factor
        theta_force = 1. + np.sin(2. * np.pi / unit_theta_cycle_duration * (stim_t - global_phase_offset -
                                                                      unit_phase_offset))
    start = int(((0.5 + track_length) * input_field_duration - peak_locs[i])/dt)
    end = start + len(stim_t)
    stim_force = gauss_force[start:end]
    stim_force = np.multiply(stim_force, theta_force)
    if syn.node.parent.parent.type == 'tuft':
        tuft_forces.append(stim_force)
    else:
        stim_forces.append(stim_force)
"""