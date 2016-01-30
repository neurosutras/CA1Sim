__author__ = 'Aaron D. Milstein'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
import random
"""
This simulation uses scipy.optimize to iterate through GABA_A_KIN mechanism parameters to fit target IPSG kinetics.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'

mech_filename = '080615 rebalanced na_ka ampa nmda - EB2'



def synaptic_kinetics_error(x, plot=0):
    """
    :param x: list of parameters
    :param plot: int or bool: method can be called manually to compare actual to target and fit waveforms
    :return: float: Error
    """
    syn.target(syn_type).kon = x[0]
    syn.target(syn_type).koff = x[1]
    syn.target(syn_type).CC = x[2]
    syn.target(syn_type).CO = x[3]
    syn.target(syn_type).Beta = x[4]
    syn.target(syn_type).Alpha = x[5]
    sim.run(v_init)

    t = np.array(sim.tvec)
    vm = np.array(sim.rec_list[0]['vec'])
    interp_t = np.arange(0., duration, 0.01)
    interp_vm = np.interp(interp_t, t, vm)
    left, right = time2index(interp_t, equilibrate-3., equilibrate-1.)
    baseline = np.average(interp_vm[left:right])
    interp_vm -= baseline
    start, end = time2index(interp_t, equilibrate, duration)
    y = interp_vm[start:end]
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
    Ro = np.array(sim.rec_list[3]['vec'])
    Rc_max = np.max(np.array(sim.rec_list[1]['vec'])+np.array(sim.rec_list[2]['vec'])+Ro)
    """
    if 4. * decay_tau > duration - equilibrate:
        steady_state = Ro[-1]
    else:
        t_steady = time2index(t, equilibrate, equilibrate + 4. * decay_tau)[1]
        steady_state = Ro[t_steady]
    if steady_state < target_val['steady_state']:
        steady_state = target_val['steady_state']  # don't penalize decay to less than target
    rise_tau = optimize.curve_fit(model_exp_rise, interp_t[1:t_peak], y[1:t_peak], p0=target_val['rise_tau'])[0]
    decay_tau = optimize.curve_fit(model_exp_decay, interp_t[t_peak+1:]-interp_t[t_peak], y[t_peak+1:],
                                     p0=target_val['decay_tau'])[0]
    """
    result = {'rise_tau': rise_tau, 'decay_tau': decay_tau, 'Rc_max': Rc_max} #  , 'steady_state': steady_state}

    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print('kon: %.3f, koff: %.3f, CC: %.3f, CO: %.3f, Beta: %.3f, Alpha: %.3f, Error: %.4E, Rise: %.3f, Decay: %.3f, '
        'Rc_max: %.3f' % (x[0], x[1], x[2], x[3], x[4], x[5], Err, rise_tau, decay_tau, Rc_max))
    if plot:
        #fit_rise = model_exp_rise(interp_t[:t_peak], rise_tau)
        #fit_decay = model_exp_decay(interp_t[:-t_peak], decay_tau)
        #target_rise = model_exp_rise(interp_t[:t_peak], target_val['rise_tau'])
        #target_decay = model_exp_decay(interp_t[:-t_peak], target_val['decay_tau'])
        plt.plot(interp_t, y, label="actual", color='b')
        #plt.plot(interp_t[:t_peak], fit_rise, label="fit", color='r')
        #plt.plot(interp_t[:-t_peak]+interp_t[t_peak], fit_decay, color='r')
        #plt.plot(interp_t[:t_peak], target_rise, label="target", color='g')
        #plt.plot(interp_t[:-t_peak]+interp_t[t_peak], target_decay, color='g')
        plt.legend(loc='best')
        plt.show()
        plt.close()
    else:
        return Err


equilibrate = 250.  # time to steady-state
duration = 1000.
v_init = -67.
num_syns = 1
spike_times = h.Vector([equilibrate])

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

cell.zero_na()

syn_type = 'GABA_A_KIN'

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


syn = Synapse(cell, cell.tree.root, [syn_type], stochastic=0)
syn.target(syn_type).gmax *= 10.
syn.source.play(spike_times)

sim.append_rec(cell, syn.node, object=syn.target(syn_type), param='_ref_g')
sim.append_rec(cell, syn.node, object=syn.target(syn_type), param='_ref_Rb')
sim.append_rec(cell, syn.node, object=syn.target(syn_type), param='_ref_Rc')
sim.append_rec(cell, syn.node, object=syn.target(syn_type), param='_ref_Ro')
sim.append_rec(cell, cell.tree.root, loc=0.)

sim.append_stim(cell, cell.tree.root, loc=0., amp=0.2, delay=0., dur=duration)

#the target values and acceptable ranges
target_val = {'rise_tau': .5, 'decay_tau': 22., 'Rc_max': 0.9}  # to hit target for exponential rise and decay of
                                                                # 0.3 ms and 10 ms, respectively
target_range = {'rise_tau': 0.01, 'decay_tau': 0.1, 'Rc_max': 0.01}

#the initial guess and bounds
#x = [kon, koff, CC, CO, Beta, Alpha)
#x0 = [10., .05, 25., 31., 2.5, 1.5]
#x0 = [60., 10., 60., 5., 100., 60.]
#x0 = [12.88, 6.47, 69.97, 6.16, 100.63, 173.04]
x0 = [5.655, 1.276, 126.608, 15.053, 105.914, 234.470] # Error: 9.1212E+02, Rise: 0.230, Decay: 22.990, Rc_max: 0.808
xmin = [.1, .1, 10., 1., 10., 10.]
xmax = [100., 100., 500., 300., 500., 500.]

take_step = Normalized_Step(x0, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)
"""
result = optimize.basinhopping(synaptic_kinetics_error, x0, niter=1200, niter_success=600, disp=True, interval=30,
                                                            minimizer_kwargs=minimizer_kwargs, take_step=take_step)
print result.x

#synaptic_kinetics_error(result.x, plot=1)

polished_result = optimize.minimize(synaptic_kinetics_error, x1, method='Nelder-Mead', options={'ftol': 1e-3,
                                                                                        'xtol': 1e-3, 'disp': True})
synaptic_kinetics_error(polished_result.x, plot=1)
"""
#synaptic_kinetics_error(x0, plot=1)
