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
new_rec_filename = '052115 apical oblique cooperativity'


def branch_cooperativity_error(x, plot=0):
    """

    :param x: list
    :return: float
    """
    start_time = time.time()
    dv['gmax'] = x[0]
    num_spines = min(50, len(parallel_optimize_branch_cooperativity_nmda_engine.spine_list))
    result = v.map_async(parallel_optimize_branch_cooperativity_nmda_engine.stim_expected, range(num_spines))
    #result = v.map_async(parallel_optimize_branch_cooperativity_nmda_engine.stim_expected, range(len(c)))
    while not result.ready():
        if time.time() % 60 < 1.:
            clear_output()
            for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
                lines = stdout.split('\n')
                if lines[-2]:
                    print lines[-2]
            sys.stdout.flush()
            time.sleep(1)
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
        if time.time() % 60 < 1.:
            clear_output()
            for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
                lines = stdout.split('\n')
                if lines[-2]:
                    print lines[-2]
            sys.stdout.flush()
            time.sleep(1)
    combine_output_files(rec_file_list, new_rec_filename+'_actual')
    for filename in glob.glob(data_dir+'out*'):
        os.remove(filename)
    with h5py.File(data_dir+new_rec_filename+'_expected.hdf5', 'r') as expected_file:
        with h5py.File(data_dir+new_rec_filename+'_actual.hdf5', 'r') as actual_file:
            sorted_sim_keys = actual_file.keys()
            sorted_sim_keys.sort(key=lambda x: len(actual_file[x].attrs['syn_indexes']))
            expected_dict, actual_dict = get_expected_vs_actual(expected_file, actual_file, sorted_sim_keys)
    expected = np.array(expected_dict['trunk'])
    actual = np.array(actual_dict['trunk'])
    supralinearity = (actual - expected) / expected * 100.
    result = {'peak_supralinearity': np.max(supralinearity)}
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print 'Parallel simulation took:', time.time()-start_time, 's, Error:', Err
    print ('gmax: %.3E' % (x[0]))
    if plot:
        plt.plot(expected, actual, label='actual')
        plt.show()
        plt.close()
        plt.plot(expected, supralinearity, label='supralinearity')
        plt.show()
        plt.close()
    else:
        return Err

#the target values and acceptable ranges
target_val = {'peak_supralinearity': 44.}
target_range = {'peak_supralinearity': 6.}

#the initial guess
x0 = [3e-4]

c = Client()
dv = c[:]
dv.clear()
dv.block = True
global_start_time = time.time()
dv.execute('from parallel_optimize_branch_cooperativity_nmda_engine import *')
v = c.load_balanced_view()
#branch_cooperativity_error(x0, plot=1)

result = optimize.minimize(branch_cooperativity_error, x0, method='Nelder-Mead', options={'xtol': 1e-5, 'ftol': 1e-3,
                                                                                          'disp': True})
branch_cooperativity_error(result.x, plot=1)
