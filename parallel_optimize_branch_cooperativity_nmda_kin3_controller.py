__author__ = 'Aaron D. Milstein'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import sys
import parallel_optimize_branch_cooperativity_nmda_kin3_engine
import os
"""
This simulation uses scipy.optimize.minimize to fit gmax_NMDA_KIN to branch cooperativity data from
Harnett et al., 2012 (performed in TTX). One apical oblique ~150 um from the soma is chosen, up to 50 spines within a
30 um long stretch of branch are stimulated. Actual and Expected EPSP amplitudes are compared, with a target peak
nonlinearity of 44 +/- 6%.

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
# new_rec_filename = '052215 apical oblique cooperativity'
# new_rec_filename = '052815 apical oblique cooperativity - proximal - new_mg - no nmda'
# new_rec_filename = '072815 apical oblique cooperativity - proximal'
# new_rec_filename = '103015 apical oblique cooperativity - NMDA_KIN3'
# new_rec_filename = '012816 altered intrinsic properties - NMDA_KIN3'
new_rec_filename = '020516 altered km2 rinp - tuning NMDA_KIN3'


def create_no_nmda_expected_file():
    """

    """
    start_time = time.time()
    dv['gmax'] = 0.
    num_spines = min(31, len(parallel_optimize_branch_cooperativity_nmda_kin3_engine.spine_list))
    result = v.map_async(parallel_optimize_branch_cooperativity_nmda_kin3_engine.stim_expected, range(num_spines))
    while not result.ready():
        time.sleep(30)
        clear_output()
        for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
        sys.stdout.flush()
    rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir+filename+'.hdf5')]
    combine_output_files(rec_file_list, new_rec_filename+'_no_nmda_expected')
    for filename in rec_file_list:
        os.remove(data_dir+filename+'.hdf5')


def branch_cooperativity_error(x, plot=0):
    """

    :param x: list
    :return: float
    """
    start_time = time.time()
    if x[0] < 0. or x[1] < 0.05 or x[1] > 0.12 or x[2] < 3. or x[2] > 10. or x[3] < 1.:
        return 1e9
    dv['gmax'] = x[0]
    dv['gamma'] = x[1]
    dv['Kd'] = x[2]
    dv['kin_scale'] = x[3]
    num_spines = min(31, len(parallel_optimize_branch_cooperativity_nmda_kin3_engine.spine_list))
    result = v.map_async(parallel_optimize_branch_cooperativity_nmda_kin3_engine.stim_expected, range(num_spines))
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
    result = v.map_async(parallel_optimize_branch_cooperativity_nmda_kin3_engine.stim_actual, instructions)
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
        units_with_nmda = []
        for group_index in expected_file:
            unit_with_nmda = get_expected_EPSP(expected_file, group_index,
                                            parallel_optimize_branch_cooperativity_nmda_kin3_engine.equilibrate,
                                            parallel_optimize_branch_cooperativity_nmda_kin3_engine.duration)
            units_with_nmda.append(unit_with_nmda['trunk'])
        unit_time = unit_with_nmda['time']
    with h5py.File(data_dir+new_rec_filename+'_no_nmda_expected.hdf5', 'r') as unit_no_nmda_file:
        units_no_nmda = []
        for group_index in unit_no_nmda_file:
            unit_no_nmda = get_expected_EPSP(unit_no_nmda_file, group_index,
                                            parallel_optimize_branch_cooperativity_nmda_kin3_engine.equilibrate,
                                            parallel_optimize_branch_cooperativity_nmda_kin3_engine.duration)
            units_no_nmda.append(unit_no_nmda['trunk'])
    unit_with_nmda = np.mean(units_with_nmda, 0)
    unit_no_nmda = np.mean(units_no_nmda, 0)
    result = {'unitary_nmda_contribution': (np.max(unit_with_nmda) - np.max(unit_no_nmda)) /
                                           np.max(unit_no_nmda) * 100.}
    with h5py.File(data_dir+new_rec_filename+'_expected.hdf5', 'r') as expected_file:
        expected_index_map = get_expected_spine_index_map(expected_file).itervalues().next()
        with h5py.File(data_dir+new_rec_filename+'_actual.hdf5', 'r') as actual_file:
            sorted_sim_keys = actual_file.keys()
            sorted_sim_keys.sort(key=lambda x: len(actual_file[x].attrs['syn_indexes']))
            expected_dict, actual_dict = get_expected_vs_actual(expected_file, actual_file, expected_index_map,
                                                                sorted_sim_keys)
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
    print '[gmax, gamma, Kd, kin_scale]: [%.3E, %.3f, %.2f, %.2f]' % (x[0], x[1], x[2], x[3])
    print 'Peak Supralinearity: %.2f, Min Supralinearity: %.2f, Unitary %% NMDA: %.3f' % \
          (result['peak_supralinearity'], result['min_supralinearity'], result['unitary_nmda_contribution'])
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
        plt.plot(unit_time, unit_with_nmda, label='with NMDA')
        plt.plot(unit_time, unit_no_nmda, label='without NMDA')
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
target_range = {'peak_supralinearity': 1., 'min_supralinearity': 0.1,'unitary_nmda_contribution': 0.5}

#the initial guess
# x = ['gmax', 'gamma', 'Kd', 'kin_scale']
# x0 = [3.13E-03, 0.091, 9.44, 2.]
# x0 = [3.992E-03, 0.100, 9.20, 1.29]
# x0 = [3.607E-03, 0.098, 7.36, 1.92]
x0 = [3.613E-03, 0.100, 7.51, 1.81]
# x1 = [4.550E-03, 0.113, 8.96, 2.63]
x1 = [3.026E-03, 0.101, 9.98, 1.83]
xmin = [5e-4, 0.05, 3., 1.]
xmax = [5e-3, 0.12, 10., 4.]

mytakestep = Normalized_Step(x0, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)

c = Client()
dv = c[:]
dv.clear()
dv.block = True
global_start_time = time.time()
dv.execute('from parallel_optimize_branch_cooperativity_nmda_kin3_engine import *')
#time.sleep(240)
v = c.load_balanced_view()
#create_no_nmda_expected_file()  # run once for each new mech_dict
"""
result = optimize.basinhopping(branch_cooperativity_error, x0, niter=720, niter_success=300, disp=True, interval=30,
                                                            minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
print result

polished_result = optimize.minimize(branch_cooperativity_error, x0, method='Nelder-Mead',
                                    options={'xtol': 1e-3, 'ftol': 1e-3, 'disp': True, 'maxiter': 400})

print polished_result
"""

#branch_cooperativity_error(result.x, plot=1)
branch_cooperativity_error(x1, 1)


"""
012916: Branch chosen with trunk origin between 75 - 100 um, spines between 30 - 60 um along ~90 um length.
[gmax, gamma, Kd, kin_scale]: [2.228E-03, 0.100, 7.65, 1.78]
Peak Supralinearity: 43.51, Min Supralinearity: -4.87, Unitary % NMDA: 0.092
Error: 1.0885E+02

020416: Branch chosen with trunk origin between ~100 um, spines between 30 - 60 um along ~90 um length.
Simplex, favoring lower gmax over min error:
[gmax, gamma, Kd, kin_scale]: [3.026E-03, 0.101, 9.98, 1.83]
Peak Supralinearity: 45.42, Min Supralinearity: -1.22, Unitary % NMDA: 12.850
Parallel simulation took 210 s, Error: 8.1046E+02
"""