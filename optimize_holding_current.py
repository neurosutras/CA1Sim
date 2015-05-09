__author__ = 'Aaron D. Milstein'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
import scipy.integrate as integrate
import random
"""
This simulation uses scipy.optimize to fit a 'holding current' amplitude to compensate for hyperpolarization in the
absence of ih.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
#mech_filename = '031815 calibrate nmda gmax.pkl'
#mech_filename = '040815 kap_kad_ampar_scale low mg kd pas no_ih no_na'
mech_filename = '050715 pas_exp_scale kdr ka_scale no_ih ampar_exp_scale - EB2'
#rec_filename = '041315 calibrate holding i_inj - EB2'


def holding_current_error(x, plot=0):
    """
    :param x: list of parameters
    :param plot: int or bool: method can be called manually to compare actual to target and fit waveforms
    :return: float: Error
    """
    print('Holding I_inj: %.6f' % (x[0]))
    sim.modify_stim(0, amp=x[0])
    sim.run(v_init)
    t = np.array(sim.tvec)
    left, right = time2index(t, duration-3., duration-1.)
    vm = np.array(sim.rec_list[0]['vec'])
    result = {'Vm_Rest': np.average(vm[left:right])}
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print 'Error:', Err, ', Vm_Rest:', result['Vm_Rest']
    if plot:
        sim.plot()
    else:
        return Err


delay = 25.
duration = 300.
amp = 0.16
v_init = -67.
num_stim = 40

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

sim = QuickSim(duration, verbose=0)
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
sim.append_rec(cell, trunk, description='trunk', loc=0.)
sim.append_stim(cell, trunk, loc=0., amp=amp, delay=delay, dur=duration-delay, description='Holding I_inj')

#the target values and acceptable ranges
target_val = {'Vm_Rest': -67.}
target_range = {'Vm_Rest': 1.}

#the initial guess and bounds
x0 = [amp]
xmin = [0.1]  # first-pass bounds
xmax = [1.]

# rewrite the bounds in the way required by optimize.minimize
xbounds = [(low, high) for low, high in zip(xmin, xmax)]

result = optimize.minimize(holding_current_error, x0, method='L-BFGS-B', bounds=xbounds,
                           options={'ftol': 1e-3, 'eps': 0.001, 'disp': True})
holding_current_error(result.x, plot=1)