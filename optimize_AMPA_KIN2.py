__author__ = 'Aaron D. Milstein'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
import random
"""
This simulation uses scipy.optimize to iterate through AMPA_KIN mechanism parameters to fit target EPSP kinetics.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'

#mech_filename = '043015 pas_exp_scale kdr ka_scale ih_sig_scale - EB2'
#mech_filename = '072515 optimized basal ka_scale dend_sh_ar_nas - EB2'
mech_filename = '103015 interim dendritic excitability ampa'


def synaptic_kinetics_error(x, plot=0):
    """
    :param x: list of parameters
    :param plot: int or bool: method can be called manually to compare actual to target and fit waveforms
    :return: float: Error
    """
    if x[0] > 100. or x[3] > 1.:
        return 1e9
    for i, syn in enumerate(stim_syn_list):
        syn.target(syn_type).kon = x[0]
        syn.target(syn_type).koff = x[1]
        syn.target(syn_type).Beta = x[2]
        syn.target(syn_type).Alpha = x[3]
    sim.run(v_init)
    t = np.array(sim.tvec)
    g = np.array(sim.rec_list[0]['vec'])
    interp_t = np.arange(0, duration, 0.001)
    interp_g = np.interp(interp_t, t, g)
    Rb = np.interp(interp_t, t, np.array(sim.rec_list[1]['vec']))
    Ro = np.interp(interp_t, t, np.array(sim.rec_list[2]['vec']))
    left, right = time2index(interp_t, equilibrate-3., equilibrate-1.)
    baseline = np.average(interp_g[left:right])
    interp_g -= baseline
    start, end = time2index(interp_t, equilibrate, duration)
    y = interp_g[start:end]
    interp_t = interp_t[start:end]
    interp_t -= interp_t[0]
    amp = np.max(y)
    t_peak = np.where(y == amp)[0][0]
    y /= amp
    rise_10 = np.where(y[0:t_peak] >= 0.1)[0][0]
    rise_90 = np.where(y[0:t_peak] >= 0.9)[0][0]
    rise_tau = interp_t[rise_90] - interp_t[rise_10]
    decay_90 = np.where(y[t_peak:] <= 0.9)[0][0]
    decay_10 = np.where(y[t_peak:] <= 0.1)[0]
    if decay_10.any():
        decay_tau = interp_t[decay_10[0]] - interp_t[decay_90]
    else:
        decay_tau = 1000.  # large error if trace has not decayed to 10% in 1 second
    Ro_peak = np.max(Ro)
    Ro_peak_loc = np.where(Ro == Ro_peak)[0][0]
    Rc_max = Ro_peak + Rb[Ro_peak_loc]
    result = {'rise_tau': rise_tau, 'decay_tau': decay_tau, 'Rc_max': Rc_max}
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print('[kon, koff, Beta, Alpha]: [%.3f, %.3f, %.3f, %.3f], Error: %.4E, Rise: %.3f, Decay: %.3f, '
        'Rc_max: %.3f' % (x[0], x[1], x[2], x[3], Err, rise_tau, decay_tau, Rc_max))
    if plot:
        plt.plot(interp_t, y)
        plt.show()
        plt.close()
    else:
        return Err


equilibrate = 250.  # time to steady-state
duration = 1250.
v_init = -67.
num_syns = 1
spike_times = h.Vector([equilibrate])

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

syn_type = 'AMPA_KIN2'

sim = QuickSim(duration)

# look for a trunk bifurcation
trunk_bifurcation = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type == 'trunk' and
                     trunk.children[1].type == 'trunk']

# get where the thickest trunk branch gives rise to the tuft
if trunk_bifurcation:  # follow the thicker trunk
    trunk = max(trunk_bifurcation[0].children[:2], key=lambda node: node.sec(0.).diam)
    trunk = (node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type for child in
                                                                                             node.children)).next()
else:
    trunk = (node for node in cell.trunk if 'tuft' in (child.type for child in node.children)).next()
tuft = (child for child in trunk.children if child.type == 'tuft').next()
trunk = trunk_bifurcation[0]

#sim.append_rec(cell, trunk, loc=1., description='trunk vm')

spine_list = []
spine_list.extend(trunk.spines)
for spine in spine_list:
    syn = Synapse(cell, spine, [syn_type], stochastic=0)
local_random = random.Random()
local_random.seed(0)
stim_syn_list = [spine_list[i].synapses[0] for i in local_random.sample(range(len(spine_list)), num_syns)]

for i, syn in enumerate(stim_syn_list):
    syn.source.play(spike_times)

sim.append_rec(cell, syn.node, object=syn.target(syn_type), param='_ref_g')
sim.append_rec(cell, syn.node, object=syn.target(syn_type), param='_ref_Rb')
sim.append_rec(cell, syn.node, object=syn.target(syn_type), param='_ref_Ro')

#the target values and acceptable ranges
target_val = {'rise_tau': .1, 'decay_tau': 7., 'Rc_max': 0.9}  # extrapolating from Chen...Murphy and Harnett...Magee
target_range = {'rise_tau': 0.01, 'decay_tau': 0.1, 'Rc_max': 0.01}

#the initial guess and bounds
#x = [kon, koff, Beta, Alpha)
#x0 = [12.88, 6.47, 10., 0.14]
x0 = [64.20, 16.19, 19.99, 0.72]
x1 = [62.88, 16.63, 20.53, 0.71]  # Error: 7.7753E+01, Rise: 0.138, Decay: 7.000, Rc_max: 0.980
xmin = [10., 1., 1., 0.1]
xmax = [100., 30., 30., 1.]

mytakestep = Normalized_Step(x0, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)
"""
result = optimize.basinhopping(synaptic_kinetics_error, x0, niter=720, niter_success=200, disp=True, interval=20,
                                                            minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
synaptic_kinetics_error(result.x, plot=1)
"""
polished_result = optimize.minimize(synaptic_kinetics_error, x0, method='Nelder-Mead', options={'ftol': 1e-3,
                                                                                        'xtol': 1e-3, 'disp': True})
synaptic_kinetics_error(polished_result.x, plot=1)

#synaptic_kinetics_error(x1, plot=1)
