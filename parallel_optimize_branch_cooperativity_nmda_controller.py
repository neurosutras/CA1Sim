__author__ = 'Aaron D. Milstein'
from IPython.parallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
import parallel_optimize_branch_cooperativity_nmda_engine
import os, glob
"""
This simulation uses scipy.optimize.minimize to fit gmax_NMDA_KIN to branch cooperativity data from
Harnett et al., 2012 (performed in TTX). One apical oblique ~150 um from the soma is chosen, up to 50 spines within a
30 um long stretch of branch are stimulated. Actual and Expected EPSP amplitudes are compared, with a target peak
nonlinearity of 44 +/- 6%.

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
#new_rec_filename = '052215 apical oblique cooperativity'
new_rec_filename = '052815 apical oblique cooperativity - proximal - new_mg - no nmda'


def branch_cooperativity_error(x, plot=0):
    """

    :param x: list
    :return: float
    """
    start_time = time.time()
    print 'gmax: %.3E, Kd: %.2f, gamma: %.3f' % (x[0], x[1], x[2])
    dv['gmax'] = x[0]
    dv['Kd'] = x[1]
    dv['gamma'] = x[2]
    num_spines = min(50, len(parallel_optimize_branch_cooperativity_nmda_engine.spine_list))
    result = v.map_async(parallel_optimize_branch_cooperativity_nmda_engine.stim_expected, range(num_spines))
    #result = v.map_async(parallel_optimize_branch_cooperativity_nmda_engine.stim_expected, range(len(c)))
    while not result.ready():
        #if time.time() % 60 < 1.:
        time.sleep(30)
        clear_output()
        for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
        sys.stdout.flush()
    rec_file_list = dv['rec_filename']
    combine_output_files(rec_file_list, new_rec_filename+'_expected')
    for filename in glob.glob(data_dir+'out*'):
        os.remove(filename)
    instructions = []
    for group in range(1, num_spines+1):
        instructions.append(range(group))
    result = v.map_async(parallel_optimize_branch_cooperativity_nmda_engine.stim_actual, instructions)
    #result = v.map_async(parallel_optimize_branch_cooperativity_nmda_engine.stim_actual, instructions[:len(c)])
    while not result.ready():
        #if time.time() % 60 < 1.:
        time.sleep(30)
        clear_output()
        for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
        sys.stdout.flush()
    combine_output_files(rec_file_list, new_rec_filename+'_actual')
    for filename in glob.glob(data_dir+'out*'):
        os.remove(filename)
    with h5py.File(data_dir+new_rec_filename+'_expected.hdf5', 'r') as expected_file:
        unit_with_nmda = get_expected_EPSP(expected_file,
                                            parallel_optimize_branch_cooperativity_nmda_engine.spine.index,
                                            parallel_optimize_branch_cooperativity_nmda_engine.equilibrate,
                                            parallel_optimize_branch_cooperativity_nmda_engine.duration)
    dv['gmax'] = 0.
    result = v.map_sync(parallel_optimize_branch_cooperativity_nmda_engine.stim_expected, [0])
    with h5py.File(data_dir+result[0]+'.hdf5', 'r') as unit_no_nmda_file:
        unit_no_nmda = get_expected_EPSP(unit_no_nmda_file,
                                            parallel_optimize_branch_cooperativity_nmda_engine.spine.index,
                                            parallel_optimize_branch_cooperativity_nmda_engine.equilibrate,
                                            parallel_optimize_branch_cooperativity_nmda_engine.duration)
    os.remove(data_dir+result[0]+'.hdf5')
    result = {'unitary_nmda_contribution': (np.max(unit_with_nmda['trunk']) - np.max(unit_no_nmda['trunk'])) /
                                           np.max(unit_no_nmda['trunk'])}
    with h5py.File(data_dir+new_rec_filename+'_expected.hdf5', 'r') as expected_file:
        with h5py.File(data_dir+new_rec_filename+'_actual.hdf5', 'r') as actual_file:
            sorted_sim_keys = actual_file.keys()
            sorted_sim_keys.sort(key=lambda x: len(actual_file[x].attrs['syn_indexes']))
            expected_dict, actual_dict = get_expected_vs_actual(expected_file, actual_file, sorted_sim_keys)
    expected = np.array(expected_dict['trunk'])
    actual = np.array(actual_dict['trunk'])
    supralinearity = (actual - expected) / expected * 100.
    peak_supralinearity = np.max(supralinearity)
    if peak_supralinearity < 0.:  # there is no gradient if integration is always sublinear
        peak_supralinearity = np.min(supralinearity)  # exaggerate error for sublinear integration
        min_supralinearity = np.min(supralinearity)
    else:
        peak_index = np.where(supralinearity==peak_supralinearity)[0][0]
        if peak_index == 0:
            min_supralinearity = supralinearity[0]
        else:
            min_supralinearity = np.min(supralinearity[:peak_index])
    result['peak_supralinearity'] = peak_supralinearity
    result['min_supralinearity'] = min_supralinearity
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print 'Parallel simulation took %i s, Error: %.4E' % (time.time()-start_time, Err)
    if plot:
        print result['peak_supralinearity']
        plt.plot(expected, actual)
        plt.xlabel('Expected EPSP (mV)')
        plt.ylabel('Actual EPSP (mV)')
        plt.title('Expected vs. Actual EPSP')
        plt.show()
        plt.close()
        plt.plot(expected, supralinearity)
        plt.xlabel('Expected EPSP (mV)')
        plt.ylabel('Supralinearity (%)')
        plt.title('Supralinearity')
        plt.show()
        plt.close()
        plt.plot(unit_with_nmda['time'], unit_with_nmda['trunk'], label='with NMDA')
        plt.plot(unit_no_nmda['time'], unit_no_nmda['trunk'], label='without NMDA')
        plt.xlabel('Time (ms)')
        plt.ylabel('EPSP Amplitude (mV)')
        plt.legend(loc='best')
        plt.title('NMDAR Contribution to Unitary EPSP')
        plt.show()
        plt.close()
    else:
        return Err

#the target values and acceptable ranges
target_val = {'peak_supralinearity': 44., 'min_supralinearity': 0., 'unitary_nmda_contribution': 0.}
target_range = {'peak_supralinearity': 6., 'min_supralinearity': 0.1,'unitary_nmda_contribution': 0.1}

#the initial guess
# x = ['gmax', 'Kd', 'gamma']
x0 = [1e-3, 3.57, 0.08]
xmin = [5e-4, 1.6, 0.06]
xmax = [5e-3, 10., 0.1]
#x1 = [3.11e-3, 7.18, 0.097]
#x1 = [3.59e-3, 8.44, 0.12]
x1 = [2.72e-3, 9.82, 0.089]  # Err: 4.18
x2 = [0., 9.82, 0.089]

blocksize = 0.5  # defines the fraction of the xrange that will be explored at each step
                 #  basinhopping starts with this value and reduces it by 10% every 'interval' iterations

mytakestep = MyTakeStep(blocksize, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)

c = Client()
dv = c[:]
dv.clear()
dv.block = True
global_start_time = time.time()
dv.execute('from parallel_optimize_branch_cooperativity_nmda_engine import *')
v = c.load_balanced_view()
"""
result = optimize.basinhopping(branch_cooperativity_error, x0, niter= 720, niter_success=100, disp=True, interval=20,
                                                            minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
branch_cooperativity_error(result.x, plot=1)


#branch_cooperativity_error(x0, plot=1)

result = optimize.minimize(branch_cooperativity_error, x1, method='Nelder-Mead', options={'xtol': 1e-5, 'ftol': 1e-3,
                                                                                          'disp': True})
branch_cooperativity_error(result.x, plot=1)
"""
branch_cooperativity_error(x2, 1)