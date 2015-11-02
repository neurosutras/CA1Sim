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
    repeat = 10
    ISI_list = [300, 100, 50, 25, 10]
    instructions = []
    for ISI in ISI_list:
        for i in range(repeat):
            instructions.append(ISI)
    dv['x'] = x
    map_result = v.map_async(parallel_optimize_pr_engine.sim_stim_train, instructions)
    while not map_result.ready():
        time.sleep(30)
        clear_output()
        for stdout in [stdout for stdout in map_result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
        sys.stdout.flush()
    dv.execute('restore_random_sequence_locations()')
    results = {}
    unit_amps = []
    interp_dt = parallel_optimize_pr_engine.interp_dt
    num_stims = parallel_optimize_pr_engine.num_stims
    for result in map_result.get():
        for ISI in result:
            if ISI == 300:
                interp_t = np.arange(0., interp_dt*len(result[ISI]), interp_dt)
                these_unit_amps = []
                for i in range(num_stims):
                    left, right = time2index(interp_t, 2.0+ISI*i, 102.0+ISI*i)
                    these_unit_amps.append(np.max(result[ISI][left:right]))
                unit_amps.append(these_unit_amps)
            if not ISI in results:
                results[ISI] = result[ISI]
            else:
                results[ISI] += result[ISI]
    for ISI in results:
        results[ISI] /= float(repeat)

    # calculate slope across multiple stims at 300 ms ISI:
    xi = np.arange(0, num_stims)
    A = np.array([xi, np.ones(num_stims)])
    unit_slopes = []
    for these_unit_amps in unit_amps:
        this_slope = np.linalg.lstsq(A.T, these_unit_amps)[0][0]
        unit_slopes.append(this_slope)
    mean_unit_slope = np.mean(unit_slopes)
    mean_unit_amp = np.mean(unit_amps)
    Err = 0.
    for ISI in results:
        if ISI == 300:
            Err += round(((target_val[ISI] - mean_unit_amp)/target_range[ISI])**2., 10)
            Err += round(((target_val['unit_slope'] - mean_unit_slope)/target_range['unit_slope'])**2., 10)
        elif ISI == 10:
            # 3rd pulse in burst of 5
            left = 0
            right = int((2.+ISI*3)/interp_dt)
            amp = np.max(results[ISI][left:right])
            #print '3rd pulse amp:', amp
            Err += round(((target_val[ISI] - amp)/target_range[ISI])**2., 10)
            # 5th pulse in burst of 5
            left = right
            right += int(ISI*2/interp_dt)
            amp = np.max(results[ISI][left:right])
            #print '5th pulse amp:', amp
            Err += round(((target_val['5th pulse'] - amp)/target_range['5th pulse'])**2., 10)
            # 6th pulse after recovery
            left = right + int(100./interp_dt)
            right = left + int(ISI/interp_dt)
            amp = np.max(results[ISI][left:right])
            #print 'recovery pulse amp:', amp
            Err += round(((target_val['recovery'] - amp)/target_range['recovery'])**2., 10)
        else:
            Err += round(((target_val[ISI] - np.max(results[ISI]))/target_range[ISI])**2., 10)
    print 'Parallel simulation took %i s, Error: %.4E' % (time.time()-start_time, Err)
    print '[Num synapses, P0, f, tau_F, D, tau_D]: [%i, %.3f, %.3f, %.3f, %.3f, %.3f], unit amp: %.3f, unit slope: ' \
          '%.2E, recovery unit amp: %.3f' % (int(x[0]*1000), x[1], x[2], x[3], x[4], x[5], mean_unit_amp,
                                                          mean_unit_slope, amp)
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
target_val = {300: 3.8, 100: 7.1, 50: 9.0, 25: 11.6, 10: 15.1, '5th pulse': 21.7, 'recovery': 5.8, 'unit_slope': 0.}
#target_range = {300: 0.6, 100: 1.0, 50: 1.4, 25: 1.7, 10: 1.8, '5th pulse': 2., 'recovery': 0.5, 'unit_slope': .01}
target_range = {300: 0.1, 100: 0.1, 50: 0.1, 25: 0.1, 10: 0.1, '5th pulse': 0.1, 'recovery': 0.1, 'unit_slope': .01}
#target_val = {300: 3.8, 100: 7.1, 50: 9.0, 25: 11.6, 10: 15.1, 'unit_slope': 0.}
#target_range = {300: 0.6, 100: 1.0, 50: 1.4, 25: 1.7, 10: 1.8, 'unit_slope': .01}

param_names = ['n', 'P0', 'f', 'tau_F', 'd1', 'tau_D1']

#the initial guess
#x0 = [0.067, 0.18, 0.92, 105.5, 0.64, 10.]  # first pass after calibrating NMDA_KIN2.gmax for cooperativity
# with 103115 interim dendritic excitability ampa nmda_kin3.pkl as mech_dict:
#x0 = [0.082, 0.233, 1.603, 46.732, 0.830, 133.609]  # Err 1.0880E+03
x0 = [0.094, 0.190, 1.625, 65.775, 0.785, 113.151]
# unit amp: 4.180, unit slope: -6.69E-02, recovery unit amp: 6.308, Error: 4.8017E+02

# the bounds
#xmin = [0.01, 0.1, 0.01, 1., 0.01, 1.]  # n will be filtered by f(n) = int(n * 1000)
#xmax = [0.3, 0.9, 50., 1e4, 1.0, 1e5]
xmin = [0.06, 0.15, 0.8, 25., 0.5, 50.]  # n will be filtered by f(n) = int(n * 1000)
xmax = [0.15, 0.3, 1.8, 150., 0.9, 300.]

if len(sys.argv) > 1:
    cluster_id = sys.argv[1]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

dv = c[:]
dv.clear()
dv.block = True
global_start_time = time.time()
dv.execute('from parallel_optimize_pr_engine import *')
v = c.load_balanced_view()

mytakestep = Normalized_Step(x0, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)
"""
result = optimize.basinhopping(release_dynamics_error, x0, niter=400, niter_success=100, disp=True, interval=30,
                                                            minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
#print result
release_dynamics_error(result.x, plot=1)

polished_result = optimize.minimize(release_dynamics_error, x0, method='Nelder-Mead',
                                    options={'xtol': 1e-3, 'ftol': 1e-3, 'maxiter': 200, 'disp': True})
#print polished_result

release_dynamics_error(polished_result.x, plot=1)
"""
release_dynamics_error(x0, plot=1)
