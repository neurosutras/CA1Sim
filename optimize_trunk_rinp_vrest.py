__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
"""
Optimizes input resistance by injecting hyperpolarizing current at soma. Modifies g_pas to fit target values from
Magee 1995, Bittner et al., 2012 in the absence of Ih. Modifies e_pas to set resting vm at -77 mV without ih.
Implements a linear gradient of g_pas along trunk, which is inherited by apical and tuft. Spines inherit from their
parent branches.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '042015 soma_pas kdr ka_scale - adjusted - EB2.pkl'
#rec_filename = 'calibrate_rinp'


def rinp_error(x, plot=0):
    """
    :param x: array
    :param plot: int
    :return: float
    """
    start_time = time.time()

    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.modify_mech_param('soma', 'pas', 'e', -x[1])
    cell.modify_mech_param('trunk', 'pas', 'g', origin='soma', slope=x[2])

    for sec_type in ['basal', 'axon', 'ais', 'trunk', 'apical', 'tuft', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')

    result = {}
    sim.modify_stim(node=cell.tree.root)
    sim.modify_rec(node=cell.tree.root)
    sim.run(v_init)
    v_rest, rinp_peak, rinp_steady = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate,
                                              duration, amp)
    result['soma'] = rinp_steady
    result['v_rest'] = v_rest
    sim.modify_stim(node=trunk)
    sim.modify_rec(node=trunk)
    sim.run(v_init)
    result['trunk'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[2]
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print('Simulation took %.3f s' % (time.time()-start_time))
    print('soma g_pas: %.4E, e_pas: %.3f, trunk slope: %.4E, Error: %.4E, soma R_inp: %i, V_rest: %.3f, trunk R_inp: %i'
          % (x[0], -x[1], x[2], Err, result['soma'], result['v_rest'], result['trunk']))
    if plot:
        sim.plot()
    else:
        return Err


equilibrate = 200.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
amp = -0.1
v_init = -77.

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
#cell = CA1_Pyr(morph_filename, full_spines=False)

cell.modify_mech_param('apical', 'pas', 'g', origin='trunk')
cell.modify_mech_param('tuft', 'pas', 'g', origin='trunk')
cell.modify_mech_param('trunk', 'pas', 'g', origin='soma', slope=1e-7)

trunk = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type in ['trunk','tuft'] and
                                            trunk.children[1].type in ['trunk', 'tuft']][0]  # tuft bifurcation

sim = QuickSim(duration, verbose=False)
sim.parameters['description'] = 'RInp (Steady-state)'
sim.append_stim(cell, cell.tree.root, 0., amp, equilibrate, stim_dur)
sim.append_rec(cell, cell.tree.root, loc=0.)

#the target values and acceptable ranges
target_val = {'soma': 110., 'v_rest': -77., 'trunk': 97.}
target_range = {'soma': 1., 'v_rest': 1., 'trunk': 1.}

#the initial guess and bounds
# x = [soma.g_pas, e_pas, trunk.g_pas slope]
x0 = [1e-5, 77., 1e-7]
xmin = [1e-7, 60., 1e-9]  # first-pass bounds
xmax = [1e-4, 90., 1e-5]
#x0 = [  2.37074536e-05,   7.51790072e+01,   2.04897295e-07]  # following first pass basinhopping
x1 = [  1.79e-05,   75.1,   3.04e-07]  # following polishing by simplex
#x1 = [1.7906E-05, 75.093, 3.0412E-07]

# rewrite the bounds in the way required by optimize.minimize
xbounds = [(low, high) for low, high in zip(xmin, xmax)]

blocksize = 0.5  # defines the fraction of the xrange that will be explored at each step
                 #  basinhopping starts with this value and reduces it by 10% every 'interval' iterations

mytakestep = MyTakeStep(blocksize, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)

result = optimize.basinhopping(rinp_error, x0, niter= 720, niter_success=100, disp=True, interval=20,
                                                    minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
rinp_error(result.x, plot=1)
#polished_result = optimize.minimize(rinp_error, result.x, method='Nelder-Mead', options={'ftol': 1e-3, 'disp': True})
#rinp_error(polished_result.x, plot=1)

#rinp_error(x1, plot=1)