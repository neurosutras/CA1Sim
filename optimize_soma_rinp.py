__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
"""
Optimizes input resistance by injecting hyperpolarizing current at soma. Modifies g_pas to fit target
values from Magee 1995, Bittner et al., 2012 in the absence of Ih.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
#mech_filename = '032315 kap_kad_ampar_scale nmda kd pas no_ih no_na.pkl'
#mech_filename = '032315 kap_kad_ih_ampar_scale nmda kd pas no_na.pkl'
mech_filename = '042015 soma_pas - EB2.pkl'
#rec_filename = 'calibrate_rinp'


def rinp_error(x, plot=0):
    """
    :param x: array
    :param plot: int
    :return: float
    """
    start_time = time.time()
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    cell.reinit_mechanisms()
    result = {}
    sim.run(v_init)
    result['soma'] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[2]
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print('Simulation took %.3f s' % (time.time()-start_time))
    print('g_pas: %.4E, Error: %.4E, R_Inp: soma: %i' % (x[0], Err, result['soma']))
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
#cell = CA1_Pyr(morph_filename, mech_filename, full_spines=False)

sim = QuickSim(duration, verbose=False)
sim.parameters['description'] = 'RInp (Steady-state)'
sim.append_stim(cell, cell.tree.root, 0., amp, equilibrate, stim_dur)
sim.append_rec(cell, cell.tree.root, loc=0.)

#the target values and acceptable ranges
target_val = {'soma': 110.}
target_range = {'soma': 1.}

#the initial guess and bounds
# x = [soma.g_pas]
x0 = [1e-6]
xmin = [1e-7]  # first-pass bounds
xmax = [1e-4]
#x0 = [8.95E-06]  # following first pass basinhopping

# rewrite the bounds in the way required by optimize.minimize
xbounds = [(low, high) for low, high in zip(xmin, xmax)]

blocksize = 0.5  # defines the fraction of the xrange that will be explored at each step
                 #  basinhopping starts with this value and reduces it by 10% every 'interval' iterations

mytakestep = MyTakeStep(blocksize, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)

"""
result = optimize.basinhopping(rinp_error, x0, niter= 720, niter_success=100, disp=True, interval=20,
                                                    minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
"""
result = optimize.minimize(rinp_error, x0, method='Nelder-Mead', options={'ftol': 1e-3, 'disp': True})
rinp_error(result.x, plot=1)