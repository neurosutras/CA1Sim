__author__ = 'Aaron D. Milstein'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
import scipy.integrate as integrate
import random
"""
This simulation uses scipy.optimize to iterate through NMDAR peak conductance to fit target NMDA\AMPA ratio by area.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
#mech_filename = '031815 calibrate nmda gmax.pkl'
mech_filename = '040815 kap_kad_ampar_scale low mg kd pas no_ih no_na.pkl'
rec_filename = '041315 calibrate_nmda_gmax - apical - EB2'


def nmda_ampa_area_error(x, plot=0):
    """
    :param x: list of parameters
    :param plot: int or bool: method can be called manually to compare actual to target and fit waveforms
    :return: float: Error
    """
    print('%s.gmax: %.6f' % (syn_types[1], x[0]))
    for syn in stim_syn_list:
        syn.target(syn_types[1]).gmax = x[0]
    sim.parameters['description'] = syn_types[0]+' + '+syn_types[1]
    sim.run(v_init)
    areas = []
    t = np.array(sim.tvec)
    left, right = time2index(t, equilibrate-2.0, equilibrate)
    vm = np.array(sim.rec_list[0]['vec'])
    baseline = np.average(vm[left:right])
    vm -= baseline
    depolarized = np.where(vm[right:] > 0)[0]
    start = depolarized[0] + right
    end = depolarized[-1] + right
    areas.append(integrate.trapz(vm[start:end], t[start:end]))
    if plot:
        with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
            sim.export_to_file(f, 0)
    sim.parameters['description'] = syn_types[0]+' alone'
    for syn in stim_syn_list:
        syn.target('NMDA_KIN').gmax = 0.
    sim.run(v_init)
    t = np.array(sim.tvec)
    left, right = time2index(t, equilibrate-2.0, equilibrate)
    vm = np.array(sim.rec_list[0]['vec'])
    baseline = np.average(vm[left:right])
    vm -= baseline
    depolarized = np.where(vm[right:] > 0)[0]
    start = depolarized[0] + right
    end = depolarized[-1] + right
    areas.append(integrate.trapz(vm[start:end], t[start:end]))
    result = {'%Area': (areas[0]-areas[1])/areas[0]}
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print 'Error:', Err, ', %Area (NMDA):', result['%Area']
    if plot:
        with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
            sim.export_to_file(f, 1)
        plot_superimpose_conditions(rec_filename)
    else:
        return Err


equilibrate = 150.  # time to steady-state
duration = 550.
v_init = -65.
num_stim = 40

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

syn_types = ['AMPA_KIN', 'NMDA_KIN']
spike_times = h.Vector([equilibrate])

# record from the trunk near where the tuft begins. to calibrate distal SR, stimulate only synapses on the trunk or
# alone apical obliques that originate <100 micron from the border of SR and SLM. to calibrate SLM, stimulate only
# synapses in tuft branches that originate <100 uM from the border.

sim = QuickSim(duration)  # , verbose=0)
trunk = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type in ['trunk','tuft'] and
                                            trunk.children[1].type in ['trunk', 'tuft']][0]  # tuft bifurcation
sim.append_rec(cell, trunk, description='trunk', loc=0.)
# Holding current calibrated to return Vm to -65
i_inj_delay = 25.
sim.append_stim(cell, trunk, loc=0., amp=0.223, delay=i_inj_delay, dur=duration-i_inj_delay,
                description='Holding I_inj')

syn_list = []
border_loc = cell.get_distance_to_node(cell.tree.root, trunk, loc=1.)

for branch in cell.trunk+cell.apical:  # cell.tuft:
    for spine in branch.spines:
        syn = Synapse(cell, spine, syn_types, stochastic=0)
        if (branch.type == 'trunk' and
            border_loc - cell.get_distance_to_node(cell.tree.root, branch, loc=syn.loc) <= 100.) or \
            (branch.type == 'apical' and
            border_loc - cell.get_distance_to_node(cell.tree.root, cell.get_dendrite_origin(branch), 1.) <= 100.) or \
            (branch.type == 'tuft' and
            cell.get_distance_to_node(cell.get_dendrite_origin(branch), branch, loc=syn.loc) <= 100.):
                #if branch.type == 'tuft':  # remove this line when probing trunk + apical synapses
                syn_list.append(syn)
cell.init_synaptic_mechanisms()

random.seed(0)
stim_syn_list = [syn_list[i] for i in random.sample(range(len(syn_list)), num_stim)]
for syn in stim_syn_list:
    syn.source.play(spike_times)

#the target values and acceptable ranges
target_val = {'%Area': 0.35}  # 0.35 for apical. 0.65 for tuft
target_range = {'%Area': 0.01}  #

#the initial guess and bounds
# x = [gmax]
x0 = [2.9e-4]
xmin = [1e-6]  # first-pass bounds
xmax = [0.1]

#mech_dict: before g_pas optimization
#x0 = [0.000145]  # result from L-BFGS-B for apical
#x0 = [0.000464]  # result from L-BFGS-B for tuft
#mech_dict: '040815 kap_kad_ampar_scale low mg kd pas no_ih no_na.pkl'
#x0 = [2.30e-4]  # result from L-BFGS-B for distal apical+trunk
#x0 = [1.39e-3]  # result from L-BFGS-B for proximal tuft

# rewrite the bounds in the way required by optimize.minimize
xbounds = [(low, high) for low, high in zip(xmin, xmax)]
"""
result = optimize.minimize(nmda_ampa_area_error, x0, method='L-BFGS-B', bounds=xbounds,
                           options={'ftol': 1e-3, 'eps': 1e-6})

print('%s.gmax: %.6f' % (syn_types[1], result.x[0]))
"""
result = optimize.minimize(nmda_ampa_area_error, x0, method='Nelder-Mead', options={'ftol': 1e-3, 'disp': True})
print result
nmda_ampa_area_error(result.x, plot=1)
"""
nmda_ampa_area_error(x0, plot=1)
"""