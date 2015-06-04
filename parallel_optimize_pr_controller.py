__author__ = 'Aaron D. Milstein'
from IPython.parallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
import parallel_optimize_pr_engine
"""
This simulation uses scipy.optimize.basinhopping to fit parameters for a stochastic release probability mechanism to
data collected by stimulating synapses at a range of inter-stimulus intervals.
Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""


def release_dynamics_error(x, plot=0):
    """

    :param x: list
    :return: float
    """
    start_time = time.time()
    ISI_list = [300, 300, 300, 100, 100, 100, 50, 50, 50, 25, 25, 25, 10, 10, 10]
    #x[0] = int(1000 * x[0])
    #print dv.block, len(v)
    dv['x'] = x
    #map_result = v.map_async(parallel_optimize_pr_engine.sim_stim_train, [300, 100])
    map_result = v.map_async(parallel_optimize_pr_engine.sim_stim_train, ISI_list)
    while not map_result.ready():
        time.sleep(30)
        clear_output()
        for stdout in [stdout for stdout in map_result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
        sys.stdout.flush()
    dv.execute('restore_random_sequence_locations()')
    #v.map_sync(parallel_optimize_pr_engine.restore_random_sequence_locations, range(len(c)))
    results = {}
    for result in map_result.get():
        for ISI in result:
            if not ISI in results:
                results[ISI] = result[ISI]
            else:
                results[ISI] += result[ISI]
    for ISI in results:
        results[ISI] /= 3.
    #return results
    Err = 0.
    for ISI in results:
        Err += round(((target_val[ISI] - np.max(results[ISI]))/target_range[ISI])**2., 10)
    print 'Parallel simulation took:', time.time()-start_time, 's, Error:', Err
    print ('Num synapses: %i, P0: %.3f, f: %.3f, tau_F: %.3f, D: %.3f, tau_D: %.3f' % (int(x[0]*1000), x[1], x[2], x[3],
                                                                                       x[4], x[5]))
    interp_t = {}
    if plot:
        interp_dt = parallel_optimize_pr_engine.interp_dt
        fig, axes = plt.subplots(5, 1)
        for i, ISI in enumerate([300, 100, 50, 25, 10]):
            interp_t[ISI] = np.arange(0., interp_dt*len(results[ISI]), interp_dt)
            interp_t[ISI] -= 2.
            axes[i].plot(interp_t[ISI], results[ISI])
            axes[i].set_title(str(ISI)+' ms ISI')
            axes[i].set_ylabel('EPSP Amp (mV)')
        axes[4].set_xlabel('Time (ms)')
        fig.subplots_adjust(hspace=0.55, wspace=0.3, left=0.05, right=0.95, top=0.95, bottom=0.08)
        plt.show()
        plt.close()
    else:
        return Err

#the target values and acceptable ranges
target_val = {300: 3.8, 100: 7.1, 50: 9.0, 25: 11.6, 10: 15.1}
target_range = {300: 0.6, 100: 1.0, 50: 1.4, 25: 1.7, 10: 1.8}

param_names = ['n', 'P0', 'f', 'tau_F', 'd1', 'tau_D1']

#the initial guess
#x0 = [0.1, 0.1, 0.2, 100.0, 0.9, 1000.0]  # n will be filtered by f(n) = int(n * 1000)
#x0 = [0.10, 0.19, 1.61, 161.6, 0.92, 3.8]
x0 = [0.09, 0.17, 1.31, 180.7, 0.80, 0.6]

# the bounds
#xmin = [0.01, 0.1, 0.01, 1., 0.01, 1.]  # n will be filtered by f(n) = int(n * 1000)
#xmax = [0.3, 0.9, 50., 1e4, 1.0, 1e5]
xmin = [0.05, 0.1, 0.1, 100., 0.5, .1]  # n will be filtered by f(n) = int(n * 1000)
xmax = [0.2, 0.4, 5., 300., 1.0, 100.]

#x1 = [0.101, 0.19, 1.61, 162.3, 0.93, 4.0]  # first pass basinhopping
x1 = [0.09, 0.17, 1.31, 180.7, 0.80, 0.6]  # second pass basinhopping

c = Client()
dv = c[:]
dv.clear()
dv.block = True
global_start_time = time.time()
dv.execute('from parallel_optimize_pr_engine import *')
v = c.load_balanced_view()

blocksize = 0.5  # defines the fraction of the xrange that will be explored at each step
                 #  basinhopping starts with this value and reduces it by 10% every 'interval' iterations

mytakestep = MyTakeStep(blocksize, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)

result = optimize.basinhopping(release_dynamics_error, x0, niter= 720, niter_success=100, disp=True, interval=20,
                                                            minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
#release_dynamics_error(result.x, plot=1)
print result
"""
polished_result = optimize.minimize(release_dynamics_error, result.x, method='Nelder-Mead', options={'ftol': 1e-3, 'disp': True})
release_dynamics_error(polished_result.x, plot=1)
release_dynamics_error(x1, plot=1)
"""