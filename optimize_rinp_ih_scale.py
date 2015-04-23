__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
"""
Optimizes input resistance by injecting hyperpolarizing current at soma. Modifies ghbar to fit target values from
Magee 1995, Bittner et al., 2012 in the presence of Ih. Implements an linear gradient of ghbar along trunk, which
is inherited by apical and tuft.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'

#mech_filename = '042115 pas_ka_ih_scale kdr - EB2.pkl'
mech_filename = '042215 pas_exp_scale kdr ka_scale ih - EB2.pkl'


def rinp_error(x, plot=0):
    """
    :param x: array
    :param plot: int
    :return: float
    """
    # target = {'r_soma', 'v_rest_trunk', 'r_trunk', 'sag'}
    # x = [soma.ghbar, trunk.ghbar slope]

    start_time = time.time()

    cell.modify_mech_param('soma', 'h', 'ghbar', x[0])
    cell.modify_mech_param('trunk', 'h', 'ghbar', origin='soma', slope=x[1])
    #cell.modify_mech_param('trunk', 'h', 'ghbar', origin='soma', slope=x[1], tau=x[2])

    for sec_type in ['basal', 'apical', 'tuft']:
        cell.reinitialize_subset_mechanisms(sec_type, 'h')

    result = {}
    sim.modify_stim(node=cell.tree.root)
    sim.modify_rec(node=cell.tree.root)
    sim.run(v_init)
    v_rest, rinp_peak, rinp_steady = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate,
                                              duration, amp)
    result['r_soma'] = rinp_steady
    #result['v_rest'] = v_rest
    sim.modify_stim(node=trunk)
    sim.modify_rec(node=trunk)
    sim.run(v_init)
    v_rest, rinp_peak, rinp_steady = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate,
                                              duration, amp)
    result['v_rest_trunk'] = v_rest
    result['r_trunk'] = rinp_steady
    result['sag'] = 100*((rinp_peak-rinp_steady)/rinp_peak)
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print('Simulation took %.3f s' % (time.time()-start_time))
    print('soma ghbar: %.4E, trunk slope: %.4E, Error: %.4E, soma R_inp: %.3f, V_rest: %.3f, '
          'trunk R_inp: %.3f, sag: %.3f' % (x[0], x[1], Err, result['r_soma'], result['v_rest_trunk'],
                                            result['r_trunk'], result['sag']))
    if plot:
        sim.plot()
    else:
        return Err


equilibrate = 200.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
amp = -0.15
v_init = -67.

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

#cell.modify_mech_param('trunk', 'h', 'ghbar', origin='soma', slope=8.57e-7, tau=100.)

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

sim = QuickSim(duration, verbose=False)
sim.parameters['description'] = 'RInp'
sim.append_stim(cell, cell.tree.root, 0., amp, equilibrate, stim_dur)
sim.append_rec(cell, cell.tree.root, loc=0.)

#the target values and acceptable ranges
target_val = {'r_soma': 66., 'v_rest_trunk': -65., 'r_trunk': 39., 'sag': 27.}
target_range = {'r_soma': 1., 'v_rest_trunk': 1., 'r_trunk': 1., 'sag': 1.}

#the initial guess and bounds
# x = [soma ghbar, trunk.ghbar_h slope]
x0 = [5e-5, 7.56e-8]
xmin = [1e-10, 1e-9]  # first-pass bounds
xmax = [1e-4, 1e-6]
#x1 =  # following first pass basinhopping, building on exp pas_scale
#x1 =  # following polishing by simplex

# rewrite the bounds in the way required by optimize.minimize
xbounds = [(low, high) for low, high in zip(xmin, xmax)]

blocksize = 0.5  # defines the fraction of the xrange that will be explored at each step
                 #  basinhopping starts with this value and reduces it by 10% every 'interval' iterations

mytakestep = MyTakeStep(blocksize, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)

result = optimize.basinhopping(rinp_error, x0, niter= 720, niter_success=100, disp=True, interval=20,
                                                    minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
#rinp_error(result.x, plot=1)

polished_result = optimize.minimize(rinp_error, result.x, method='Nelder-Mead', options={'ftol': 1e-3, 'disp': True})
rinp_error(polished_result.x, plot=1)
