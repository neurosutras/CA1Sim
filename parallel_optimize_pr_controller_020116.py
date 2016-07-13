__author__ = 'Aaron D. Milstein'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
import parallel_optimize_pr_engine_020116
"""
This simulation uses scipy.optimize.basinhopping to fit parameters for a stochastic release probability mechanism to
data collected by stimulating synapses at a range of inter-stimulus intervals.
Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""


def basal_release_error(x, plot=0):
    """

    :param x: list
    :return: float
    """
    start_time = time.time()
    repeat = 20
    instructions = []
    for i in range(repeat):
        instructions.append(x)
    map_result = v.map_async(parallel_optimize_pr_engine_020116.sim_stim_train_basal, instructions)
    while not map_result.ready():
        time.sleep(30)
        clear_output()
        for stdout in [stdout for stdout in map_result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
        sys.stdout.flush()
    unit_amps = []
    interp_dt = parallel_optimize_pr_engine_020116.interp_dt
    num_stims = parallel_optimize_pr_engine_020116.num_stims
    ISI = 300.
    mean_trace = None
    for result in map_result.get():
        interp_t = np.arange(0., interp_dt*len(result), interp_dt)
        these_unit_amps = []
        for i in range(num_stims):
            left, right = time2index(interp_t, 2.0+ISI*i, 102.0+ISI*i)
            these_unit_amps.append(np.max(result[left:right]))
        unit_amps.append(these_unit_amps)
        if mean_trace is None:
            mean_trace = np.array(result)
        else:
            mean_trace += result
    mean_trace /= float(repeat)
    mean_unit_amp = np.mean(unit_amps)
    Err = 0.
    Err += round(((target_val[300] - mean_unit_amp)/target_range[300])**2., 10)
    P0 = parallel_optimize_pr_engine_020116.P0
    print 'Parallel simulation took %i s, Error: %.4E' % (time.time()-start_time, Err)
    print '[Num synapses, P0]: [%i, %.3f], unit amp: %.3f' % (int(x[0]*10000.), P0, mean_unit_amp)
    if plot:
        interp_dt = parallel_optimize_pr_engine_020116.interp_dt
        interp_t = np.arange(0., interp_dt*len(mean_trace), interp_dt)
        interp_t -= 2.
        plt.plot(interp_t, mean_trace)
        plt.title('300 ms ISI')
        plt.ylabel('EPSP Amp (mV)')
        plt.xlabel('Time (ms)')
        plt.show()
        plt.close()
    return Err


def release_dynamics_error(x, plot=0):
    """

    :param x: list
    :return: float
    """
    start_time = time.time()
    repeat = 20
    ISI_list = [300, 100, 50, 25, 10]
    instructions = []
    for ISI in ISI_list:
        for i in range(repeat):
            instructions.append(ISI)
    dv['x'] = x
    map_result = v.map_async(parallel_optimize_pr_engine_020116.sim_stim_train_dynamics, instructions)
    while not map_result.ready():
        time.sleep(30)
        clear_output()
        for stdout in [stdout for stdout in map_result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
        sys.stdout.flush()
    results = {}
    unit_amps = []
    interp_dt = parallel_optimize_pr_engine_020116.interp_dt
    num_stims = parallel_optimize_pr_engine_020116.num_stims
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
            print 'ISI:', ISI, 'Amp:', mean_unit_amp
            Err += round(((target_val[ISI] - mean_unit_amp)/target_range[ISI])**2., 10)
            Err += round(((target_val['unit_slope'] - mean_unit_slope)/target_range['unit_slope'])**2., 10)
        elif ISI == 10:
            # 3rd pulse in burst of 5
            left = 0
            right = int((2.+ISI*3)/interp_dt)
            amp = np.max(results[ISI][left:right])
            print 'ISI:', ISI, 'Amp:', amp
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
            amp = np.max(results[ISI])
            print 'ISI:', ISI, 'Amp:', amp
            Err += round(((target_val[ISI] - amp)/target_range[ISI])**2., 10)
    N = int(x0['basal'][0] * 10000.)
    P0 = parallel_optimize_pr_engine_020116.P0
    print 'Parallel simulation took %i s, Error: %.4E' % (time.time()-start_time, Err)
    print '[Num synapses, P0, f, tau_F, d, tau_D]: [%i, %.3f, %.3f, %.3f, %.3f, %.3f], unit amp: %.3f, unit slope: ' \
          '%.3E, recovery unit amp: %.3f' % (N, P0, x[0], x[1], x[2], x[3], mean_unit_amp, mean_unit_slope, amp)
    interp_t = {}
    if plot:
        interp_dt = parallel_optimize_pr_engine_020116.interp_dt
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
    return Err, results


#the target values and acceptable ranges
target_val = {300: 3.8, 100: 7.1, 50: 9.0, 25: 11.6, 10: 15.1, '5th pulse': 21.7, 'recovery': 5.8, 'unit_slope': 0.}
#target_range = {300: 0.6, 100: 1.0, 50: 1.4, 25: 1.7, 10: 1.8, '5th pulse': 2., 'recovery': 0.5, 'unit_slope': .01}
target_range = {300: 0.1, 100: 0.1, 50: 0.1, 25: 0.1, 10: 0.1, '5th pulse': 0.1, 'recovery': 0.1, 'unit_slope': .01}
#target_val = {300: 3.8, 100: 7.1, 50: 9.0, 25: 11.6, 10: 15.1, 'unit_slope': 0.}
#target_range = {300: 0.6, 100: 1.0, 50: 1.4, 25: 1.7, 10: 1.8, 'unit_slope': .01}

#the initial guess
x0 = {}
xmin = {}
xmax = {}

# x0['basal'] = [num_syns_to_stim]
x0['basal'] = [0.00301]  # n will be filtered by f(n) = int(n * 10000)
xmin['basal'] = [0.0030]
xmax['basal'] = [0.0040]

# x0['dynamics'] = ['f', 'tau_F', 'd1', 'tau_D1']
# x0['dynamics'] = [1.625, 65.775, 0.785, 113.151]
# unit amp: 4.180, unit slope: -6.69E-02, recovery unit amp: 6.308, Error: 4.8017E+02
# x0['dynamics'] = [1.749, 70.195, 0.875, 88.964]
# unit amp: 4.444, unit slope: 3.817E-02, recovery unit amp: 7.101, basinhopping step 273: Error: 1323.4
# x0['dynamics'] = [1.786, 66.242, 0.883, 105.858]
# unit amp: 3.794, unit slope: -3.023E-02, recovery unit amp: 7.664, basinhopping step 199: Error: 729.99
x0['dynamics'] = [1.769, 67.351, 0.878, 92.918]
# the bounds
xmin['dynamics'] = [0.8, 25., 0.5, 50.]
xmax['dynamics'] = [1.8, 150., 0.9, 300.]

if len(sys.argv) > 1:
    cluster_id = sys.argv[1]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()

dv = c[:]
dv.clear()
dv.block = True
global_start_time = time.time()
result = dv.execute('from parallel_optimize_pr_engine_020116 import *')
while not result.ready():
    time.sleep(30)
v = c.load_balanced_view()

"""
# first optimize basal_release_error once in order to determine the ideal number of synapses to produce the target
# single pulse EPSP with the specified basal release probability
result = optimize.minimize(basal_release_error, x0['basal'], method='Nelder-Mead',
                                    options={'xtol': 1e-3, 'ftol': 1e-3, 'maxiter': 100, 'disp': True})
print result
"""


# then use the resulting number of synapses during optimization of the parameters governing release dynamics
N = int(x0['basal'][0] * 10000.)
result = dv.execute('syn_list.choose_syns_to_stim('+str(N)+')')
while not result.ready():
    time.sleep(30)


"""
mytakestep = Normalized_Step(x0['dynamics'], xmin['dynamics'], xmax['dynamics'])
minimizer_kwargs = dict(method=null_minimizer)

result = optimize.basinhopping(release_dynamics_error, x0['dynamics'], niter=720, niter_success=200, disp=True,
                                interval=30, minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
print result
"""
"""
# release_dynamics_error(result.x, plot=1)

polished_result = optimize.minimize(release_dynamics_error, x0['dynamics'], method='Nelder-Mead',
                                    options={'xtol': 1e-3, 'ftol': 1e-3, 'maxiter': 200, 'disp': True})
print polished_result

release_dynamics_error(polished_result.x, plot=1)
"""

Err, results = release_dynamics_error(x0['dynamics'], plot=0)

"""
020716:
[Num synapses, P0, f, tau_F, d, tau_D]: [30, 0.200, 1.786, 66.242, 0.883, 105.858],
unit amp: 3.794, unit slope: -3.023E-02, recovery unit amp: 7.664, basinhopping step 199: Error: 729.99

[Num synapses, P0, f, tau_F, d, tau_D]: [30, 0.200, 1.769, 67.351, 0.878, 92.918],
unit amp: 3.782, unit slope: 8.499E-02, recovery unit amp: 7.519, basinhopping step 313: Error: 687.832
"""