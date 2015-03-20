__author__ = 'Aaron D. Milstein'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
"""
This simulation uses scipy.optimize to iterate through AMPAR gating parameters to fit target mEPSP kinetics.
"""
#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '031215 kap_kad_ih_ampar_scale nmda kd no_na.pkl'


def null_minimizer(fun, x0, args, **options):
    """
    Rather than allow basinhopping to pass each local mimimum to a gradient descent algorithm for polishing, this method
    just catches and passes all local minima so basinhopping can proceed.
    """
    return optimize.OptimizeResult(x=x0, fun=fun(x0, *args), success=True, nfev=1)


class MyTakeStep(object):
    """
    Converts basinhopping absolute stepsize into different stepsizes for each parameter such that the stepsizes are
    some fraction of the ranges specified by xmin and xmax. Also enforces bounds for x.
    """
    def __init__(self, blocksize, xmin, xmax, stepsize=0.5):
        self.stepsize = stepsize
        self.blocksize = blocksize
        self.xmin = xmin
        self.xmax = xmax
        self.xrange = []
        for i in range(len(self.xmin)):
            self.xrange.append(self.xmax[i] - self.xmin[i])

    def __call__(self, x):
        for i in range(len(x)):
            if x[i] < self.xmin[i]:
                x[i] = self.xmin[i]
            if x[i] > self.xmax[i]:
                x[i] = self.xmax[i]
            snew = self.stepsize / 0.5 * self.blocksize * self.xrange[i] / 2.
            sinc = min(self.xmax[i] - x[i], snew)
            sdec = min(x[i]-self.xmin[i], snew)
            x[i] += np.random.uniform(-sdec, sinc)
        return x


def synaptic_kinetics_error(x, plot=0):
    """
    :param x: list of parameters
    :param plot: int or bool: method can be called manually to compare actual to target and fit waveforms
    :return: float: Error
    """
    print('kon: %.2f, koff: %.2f, CC: %.2f, CO: %.2f, Beta: %.2f, Alpha: %.2f' % (x[0], x[1], x[2], x[3], x[4], x[5]))

    syn.target(syn_type).kon = x[0]
    syn.target(syn_type).koff = x[1]
    syn.target(syn_type).CC = x[2]
    syn.target(syn_type).CO = x[3]
    syn.target(syn_type).Beta = x[4]
    syn.target(syn_type).Alpha = x[5]
    sim.run(v_init)

    interp_t = np.arange(0, duration, 0.001)
    t = np.array(sim.tvec)
    left, right = time2index(t, equilibrate-2.0, equilibrate)
    vm = np.array(sim.rec_list[0]['vec'])
    baseline = np.average(vm[left:right])
    vm -= baseline
    interp_vm = np.interp(interp_t, t, vm)
    Rovec = np.array(sim.rec_list[2]['vec'])
    Ro = np.max(Rovec)
    Rbvec = np.array(sim.rec_list[1]['vec']) + Rovec
    Rb = np.max(Rbvec)
    start, end = time2index(interp_t, equilibrate, duration)
    y = interp_vm[start:end]
    amp = np.max(y)
    t_peak = np.where(y == amp)[0][0]
    y /= amp
    interp_t = interp_t[start:end]
    interp_t -= interp_t[0]
    rise_tau = optimize.curve_fit(model_exp_rise, interp_t[1:t_peak], y[1:t_peak], p0=target_val['rise_tau'])[0]
    decay_tau = optimize.curve_fit(model_exp_decay, interp_t[t_peak+1:]-interp_t[t_peak], y[t_peak+1:],
                                     p0=target_val['decay_tau'])[0]
    result = {'Ro': Ro, 'rise_tau': rise_tau, 'decay_tau': decay_tau}

    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print 'Error:', Err, 'Rise:', rise_tau, 'Decay:', decay_tau, '1-Ru:', Rb, 'Ro:', Ro, 'Amp:', amp
    if plot:
        fit_rise = model_exp_rise(interp_t[:t_peak], rise_tau)
        fit_decay = model_exp_decay(interp_t[:-t_peak], decay_tau)
        target_rise = model_exp_rise(interp_t[:t_peak], target_val['rise_tau'])
        target_decay = model_exp_decay(interp_t[:-t_peak], target_val['decay_tau'])
        plt.plot(interp_t, y, label="actual", color='b')
        plt.plot(interp_t[:t_peak], fit_rise, label="fit", color='r')
        plt.plot(interp_t[:-t_peak]+interp_t[t_peak], fit_decay, color='r')
        plt.plot(interp_t[:t_peak], target_rise, label="target", color='g')
        plt.plot(interp_t[:-t_peak]+interp_t[t_peak], target_decay, color='g')
        plt.plot(t[right:]-equilibrate, Rbvec[right:], label="1-Ru", color='m')
        plt.plot(t[right:]-equilibrate, Rovec[right:], label="Ro", color='c')
        plt.legend(loc='best')
        plt.show()
        plt.close()
    else:
        return Err


equilibrate = 150.  # time to steady-state
duration = 225.
v_init = -65.

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=False)
cell.insert_spines_in_subset(['trunk', 'apical'])
syn_type = 'AMPA_KIN'

sim = QuickSim(duration)
branch = cell.get_node_by_distance_to_soma(100, 'trunk')
head = branch.spines[0]
sim.append_rec(cell, branch, 0.5, description='trunk_Vm')
syn = Synapse(cell, head, [syn_type], stochastic=0)
spike_times = h.Vector([equilibrate])
syn.source.play(spike_times)

sim.append_rec(cell, head, param='_ref_Rb', object=syn.target(syn_type), description='Rb')
sim.append_rec(cell, head, param='_ref_Ro', object=syn.target(syn_type), description='Ro')

#the target values and acceptable ranges
target_val = {'Ro': 0.6, 'rise_tau': 0.2, 'decay_tau': 5.}
target_range = {'Ro': 0.06, 'rise_tau': 0.02, 'decay_tau': 0.5}

#the initial guess and bounds
# x = [kon, koff, Beta, Alpha]
#x0 = [10., 50., 25., 2.5, 10., 6.]  # Milstein 2007
#xmin = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]  # first-pass bounds
#xmax = [100., 100., 100., 100., 100., 100.]
x0 = [64.86, 31.92, 93.75, 8.33, 53.76, 32.99]  # result from blocksize = 0.5
xmin = [10., 10., 10., 0.1, 10., 10.]  # second-pass bounds did not change result after 100 iterations.
xmax = [200., 200., 200., 100., 200., 200.]

# rewrite the bounds in the way required by optimize.minimize
xbounds = [(low, high) for low, high in zip(xmin, xmax)]

blocksize = 0.5  # defines the fraction of the xrange that will be explored at each step
                 #  basinhopping starts with this value and reduces it by 10% every 'interval' iterations

mytakestep = MyTakeStep(blocksize, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)

result = optimize.basinhopping(synaptic_kinetics_error, x0, niter= 720, niter_success=100, disp=True, interval=20,
                                                                minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)

print('kon, koff, CC, CO, Beta, Alpha: %.2f, %.2f, %.2f, %.2f, %.2f, %.2f' % (result.x[0], result.x[1], result.x[2],
                                                                              result.x[3], result.x[4], result.x[5]))
synaptic_kinetics_error(result.x, plot=1)
"""
synaptic_kinetics_error(x0, plot=1)

result2 = optimize.minimize(synaptic_kinetics_error, result.x, method='L-BFGS-B', bounds=xbounds)
"""