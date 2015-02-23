__author__ = 'Aaron D. Milstein'
import time
from plot_results import *
from specify_cells import *
from function_lib import *
from neuron import h
import scipy as sp
"""
This simulation uses scipy.optimize to iterate through NMDAR gating parameters to fit target EPSP kinetics.
"""

morph_filename = 'MC120914100xC3-scaled.swc'
mech_filename = '022315 kap_scale kd ih_scale no_na.pkl'


def diff_of_exp_model_func(t, amp, rise_tau, decay_tau):
    """
    Returns a difference of exponentials waveform with specified rise and decay kinetics.
    :param t: :class:'np.array'
    :param amp: float
    :param rise_tau: float
    :param decay_tau: float
    :return: :class:'np.array'
    """
    return np.round(amp*(np.exp(t/decay_tau)-np.exp(t/rise_tau)), 10)


def fit_exp_nonlinear(t, y, rise, decay):
    """
    Fits the input vectors to a difference of exponentials and returns the fit parameters.
    :param t:
    :param y:
    :param rise:
    :param decay:
    :return:
    """
    opt_parms, parm_cov = sp.optimize.curve_fit(diff_of_exp_model_func, t, y, p0=[1., rise, decay], maxfev=2000)
    A1, tau1, tau2 = opt_parms
    return A1, tau1, tau2


def null_minimizer(fun, x0, args, **options):
    """
    Rather than allow basinhopping to pass each local mimimum to a gradient descent algorithm for polishing, this method
    just catches and passes all local minima so basinhopping can proceed.
    """
    return sp.optimize.OptimizeResult(x=x0, fun=fun(x0, *args), success=True, nfev=1)


class MyBounds(object):
    """
    Implementation of parameters bounds used by some sp.optimize minimizers.
    """
    def __init__(self, xmin, xmax):
        self.xmin = np.array(xmin)
        self.xmax = np.array(xmax)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin


class MyTakeStep(object):
    """
    Converts basinhopping absolute stepsize into different stepsizes for each parameter such that the stepsizes are
    some fraction of the ranges specified by xmin and xmax.
    """
    def __init__(self, stepsize, blocksize, xmin, xmax):
        self.stepsize = stepsize
        self.blocksize = blocksize
        self.xmin = xmin
        self.xmax = xmax
        self.xrange = []
        for i in range(len(xmin)):
            self.xrange.append(self.xmax[i] - self.xmin[i])

    def __call__(self, x):
        for i in range(len(x)):
            snew = self.stepsize * self.blocksize * self.xrange[i]/0.5
            sinc = min(xmax[i] - x[i], snew)
            sdec = min(x[i]-xmin[i], snew)
            x[i] += np.random.uniform(-sdec, sinc)
        return x


def nmda_kinetics_error(x, syn_type, equilibrate, target_val, target_range, plot=0):
    """
    :param x: list of parameters
    :param syn_type: str: name of point_process
    :param equilibrate: float: time from start of simulation to steady-state
    :param target_val: dict: contains target values for features extracted from simulation
    :param target_range: dict: contains acceptable ranges around target values
    :param plot: int or bool: method can be called manually to compare actual to target and fit waveforms
    :return: float: error
    """
    syn_type = syn_type
    equilibrate = equilibrate
    target_val = target_val
    target_range = target_range

    syn.target(syn_type).kon = x[0]
    syn.target(syn_type).koff = x[1]
    syn.target(syn_type).Beta = x[2]
    syn.target(syn_type).Alpha = x[3]
    syn.target(syn_type).Delta = x[4]
    syn.target(syn_type).Gamma = x[5]
    sim.run()
    t = np.array(sim.tvec)
    start = np.where(t >= equilibrate)[0][0]
    right = start - 1
    left = right - 1
    while t[right]-t[left] < 2:
        left -= 1
    vm = np.array(sim.rec_list[0]['vec'])
    baseline = np.average(vm[left:right])
    Rovec = np.array(sim.rec_list[2]['vec'])
    Ro = np.max(Rovec)
    Rbvec = np.array(sim.rec_list[1]['vec']) + Rovec + np.array(sim.rec_list[3]['vec'])
    Rb = np.max(Rbvec)
    y = vm[start:] - baseline
    amp = np.max(y)
    t = t[start:]
    t -= t[0]
    exp_amp, rise_tau, decay_tau = fit_exp_nonlinear(t, y, target_val['rise_tau'], target_val['decay_tau'])
    print('kon, koff, Beta, Alpha, Delta, Gamma: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f' % (x[0], x[1], x[2], x[3], x[4],
                                                                                        x[5]))
    print 'Rise tau: ', rise_tau,' , Decay tau: ', decay_tau, ', 1-Ru: ', Rb, ', Ro: ', Ro, ', Amp: ', amp
    result = {'Ro': Ro, 'rise_tau': rise_tau, 'decay_tau': decay_tau}  # '1-Ru': Rb,
    if plot:
        y /= amp
        fit_y = diff_of_exp_model_func(t, exp_amp, rise_tau, decay_tau)
        fit_y /= np.max(fit_y)
        target_y = diff_of_exp_model_func(t, 1., target_val['rise_tau'], target_val['decay_tau'])
        target_y /= np.max(target_y)
        plt.plot(t, y, label="actual")
        plt.plot(t, fit_y, label="fit")
        plt.plot(t, target_y, label="target")
        plt.plot(t, Rbvec[start:], label="1-Ru")
        plt.plot(t, Rovec[start:], label="Ro")
        plt.legend(loc='best')
        plt.show()
        plt.close()
    else:
        Err = 0.
        for target in result:
            Err += round(((target_val[target] - result[target])/target_range[target])**2., 10)
        return Err


equilibrate = 200.  # time to steady-state
duration = 600.

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=False)
syn_type = 'NMDA_D'
for node in cell.basal+cell.trunk+cell.apical+cell.tuft:
    cell.insert_spine(node, 0.5)
cell._reinit_mech(cell.spine, 1)

sim = QuickSim(duration)
branch = cell.get_node_by_distance_to_soma(200, 'apical')
head = branch.spines[0]
sim.append_rec(cell, head, 0.5, description='spine_Vm')
syn = Synapse(cell, head, [syn_type], stochastic=0)
syn.target(syn_type).mg = 0.1
spike_times = h.Vector([equilibrate])
syn.source.play(spike_times)
sim.append_rec(cell, head, param='_ref_Rb', object=syn.target(syn_type), description='Rb')
sim.append_rec(cell, head, param='_ref_Ro', object=syn.target(syn_type), description='Ro')
sim.append_rec(cell, head, param='_ref_Rd', object=syn.target(syn_type), description='Rd')

#the target values and acceptable ranges
target_val = {'Ro': 0.25, 'rise_tau': -7., 'decay_tau': -70.} # '1-Ru': 0.9,
target_range = {'Ro': 0.025, 'rise_tau': 0.5, 'decay_tau': 5}  # '1-Ru': 0.025,

#the initial guess
# x = [kon, koff, Beta, Alpha, Delta, Gamma]
x0 = [7.5, .01, .0838, .0838, .0152, .0094]  # Lester and Jahr
#x0 = [5.5355, 0.0220, 0.0680, 0.0986, 0.0343, 0.7161]  # result from blocksize = 0.5, niter = 720, interval = 20

# the bounds
stepsize = 0.5  # basinhopping absolute stepsize
blocksize = 0.5  # MyTakeStep relative stepsize
xmin = [1., .005, .005, .005, .0005, .001]
xmax = [10., .1, .2, .5, .1, 2.]
mytakestep = MyTakeStep(stepsize, blocksize, xmin, xmax)

# rewrite the bounds in the way required by L-BFGS-B
bounds = [(low, high) for low, high in zip(xmin, xmax)]

# use method L-BFGS-B
#minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds, tol=1.)

minimizer_kwargs = dict(method=null_minimizer, args=(syn_type, equilibrate, target_val, target_range))
#mybounds = MyBounds(xmin, xmax)

result = sp.optimize.basinhopping(nmda_kinetics_error, x0, niter= 720, minimizer_kwargs=minimizer_kwargs, disp=True,
                                  interval=20, take_step=mytakestep)
print('kon, koff, Beta, Alpha, Delta, Gamma: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f' % (result.x[0], result.x[1],
                                                                                    result.x[2], result.x[3],
                                                                                    result.x[4], result.x[5]))
nmda_kinetics_error(result.x, syn_type, equilibrate, target_val, target_range, plot=1)
"""
result2 = sp.optimize.minimize(nmda_kinetics_error, result.x, method='L-BFGS-B', args=(syn_type, equilibrate,
                                                                                       target_val, target_range),
                               bounds=bounds)
"""