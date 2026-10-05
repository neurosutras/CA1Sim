__author__ = 'Aaron D. Milstein'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
import parallel_optimize_branch_cooperativity_na_engine
import os
"""
This simulation uses scipy.optimize.minimize to fit dendritic nas.gmax to branch cooperativity data from
Losonczy et al., 2006. One apical oblique ~150 um from the soma is chosen, up to 50 spines within a
30 um long stretch of branch are stimulated. Actual and Expected EPSP amplitudes are compared, with a target peak
nonlinearity of 44 +/- 6%.

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
new_rec_filename = '110515 apical oblique cooperativity - na'


def branch_cooperativity_error(x, plot=0):
    """

    :param x: list
    :return: float
    """
    start_time = time.time()
    num_spines = min(30, len(parallel_optimize_branch_cooperativity_na_engine.spine_list))
    result = v.map_async(parallel_optimize_branch_cooperativity_na_engine.stim_expected, range(num_spines))
    while not result.ready():
        time.sleep(30)
        clear_output()
        for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
        sys.stdout.flush()
    rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir+filename+'.hdf5')]
    combine_output_files(rec_file_list, new_rec_filename+'_expected')
    for filename in rec_file_list:
        os.remove(data_dir+filename+'.hdf5')
    instructions = []
    for group in range(1, num_spines+1):
        instructions.append(range(group))
    result = v.map_async(parallel_optimize_branch_cooperativity_na_engine.stim_actual, instructions)
    while not result.ready():
        time.sleep(30)
        clear_output()
        for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
        sys.stdout.flush()
    rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir+filename+'.hdf5')]
    combine_output_files(rec_file_list, new_rec_filename+'_actual')
    for filename in rec_file_list:
        os.remove(data_dir+filename+'.hdf5')
    with h5py.File(data_dir+new_rec_filename+'_expected.hdf5', 'r') as expected_file:
        expected_index_map = get_expected_spine_index_map(expected_file).itervalues().next()
        with h5py.File(data_dir+new_rec_filename+'_actual.hdf5', 'r') as actual_file:
            sorted_sim_keys = actual_file.keys()
            sorted_sim_keys.sort(key=lambda x: len(actual_file[x].attrs['syn_indexes']))
            expected_dict, actual_dict = get_expected_vs_actual(expected_file, actual_file, expected_index_map,
                                                                sorted_sim_keys)
    expected = np.array(expected_dict['soma'])
    actual = np.array(actual_dict['soma'])
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
    result = {'peak_supralinearity': peak_supralinearity, 'min_supralinearity': min_supralinearity}
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print 'Peak Supralinearity: %.2f, Min Supralinearity: %.2f' % \
          (result['peak_supralinearity'], result['min_supralinearity'])
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
    else:
        return Err

#the target values and acceptable ranges
target_val = {'peak_supralinearity': 44., 'min_supralinearity': 0.}
target_range = {'peak_supralinearity': 1., 'min_supralinearity': 1.}

#the initial guess
x0 = [3.613E-03, 0.100, 7.51, 1.81]
xmin = [1e-4, 0.05, 3., 1.]
xmax = [5e-3, 0.1, 10., 4.]

mytakestep = Normalized_Step(x0, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)

c = Client()
dv = c[:]
dv.block = True
global_start_time = time.time()
dv.execute('from parallel_optimize_branch_cooperativity_na_engine import *')
v = c.load_balanced_view()

"""
result = optimize.basinhopping(branch_cooperativity_error, x0, niter=700, niter_success=200, disp=True, interval=20,
                                                            minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
branch_cooperativity_error(result.x, plot=1)

result = optimize.minimize(branch_cooperativity_error, x0, method='Nelder-Mead', options={'xtol': 1e-3, 'ftol': 1e-3,
                                                                                    'disp': True, 'maxiter': 200})
branch_cooperativity_error(result.x, plot=1)
"""
branch_cooperativity_error(x0, 1)
