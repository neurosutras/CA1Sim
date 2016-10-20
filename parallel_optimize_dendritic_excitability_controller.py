__author__ = 'Grace Ng'
from ipyparallel import Client
from IPython.display import clear_output
from plot_results import *
import scipy.optimize as optimize
import sys
import parallel_optimize_dendritic_excitability_engine
import os
"""
Aims for spike initiation at initial segment by increasing nax density and decreasing activation V1/2 relative to soma,
axon_hill, and axon compartments. Extend linear kap gradient into basals and obliques, aim for 60% spike attenuation
at bifurcation of trunk and tuft.

Hierarchical optimization:
I) optimize g_pas for target rinp at soma, trunk bifurcation, and tuft bifurcation [without h].
II) optimize ghbar_h for target rinp at soma, trunk bifurcation, and tuft bifurcation, while also optimizing for v_rest
offset between soma and tuft, and EPSP shape changes between proximal and distal synapses measured at the soma.
III) optimize gbar_nax/nas/sh/sha, gkabar_kap/d, gkdrbar for target na spike threshold, AHP amp, and vm stability

Parallel version dynamically submits jobs to available cores.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

new_rec_filename = '101116 altered km2 rinp - tuning pas'


class history(object):
    def __init__(self):
        self.x_values = []
        self.error_values = []
        self.Rinp_values = dict([('soma', []), ('trunk', []), ('tuft', [])])

    def find_lowest_error(self):
        lowest_Err = float(self.error_values[0])
        best_x = self.x_values[0]
        for index, error in enumerate(self.error_values):
            if error < lowest_Err:
                lowest_Err = float(error)
                best_x = self.x_values[index]
        return lowest_Err

"""
#write a method in the class that can look up the lowest error and return the corresponding x value that led to it
global lowest_Err
global best_x
if Err < lowest_Err:
    lowest_Err = float(Err)
    best_x = list(x)
"""

hist = history()


def pas_error(x, plot=0):
    """
    :param x: array (soma.g_pas, trunk.g_pas slope, trunk.g_pas tau)
    :param plot: int
    :return: float
    """
    if np.any(x < 0.):
        return 1e9
    start_time = time.time()
    dv['x'] = x
    #dv.push(xmin['pas'])
    #dv.push(xmax['pas'])

    hist.x_values.append(x)

    #parallel_optimize_dendritic_excitability_engine.node_pass_error
    sec_list = ['soma', 'trunk', 'tuft']
    print 'Process %i using current x: [soma g_pas, trunk slope, trunk tau]: [%.2E, %.2E, %.1f]' % (os.getpid(), x[0],
                                                                                                    x[1], x[2])
    result = v.map_async(parallel_optimize_dendritic_excitability_engine.section_pas_error, sec_list)

    while not result.ready():
        # time.sleep(5)
        clear_output()
        for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
            lines = stdout.split('\n')
            if lines[-2]:
                print lines[-2]
        sys.stdout.flush()
    result = result.get()

    Err = 0.

    final_result = {}
    for dict in result:
        final_result.update(dict)
    for section in final_result:
        Err += ((target_val['pas'][section] - final_result[section]) / target_range['pas'][section]) ** 2.
        hist.Rinp_values[section].append(final_result[section])
    hist.error_values.append(Err)

    print('Simulation took %.3f s' % (time.time() - start_time))
    print 'Process %i: [soma g_pas, trunk slope, trunk tau]: [%.2E, %.2E, %.1f], soma R_inp: %.1f, ' \
          'trunk R_inp: %.1f, tuft R_inp: %.1f' % (os.getpid(), x[0], x[1], x[2], final_result.get('soma'),
                                                   final_result.get('trunk'), final_result.get('tuft'))
    print 'Process %i: Error: %.4E' % (os.getpid(), Err)
    """
    if plot:
        sim.plot()
    """
    return Err


#variable values from optimize_dendritic_excitability_020416
v_init = -67.
soma_ek = -77.

#the target values and acceptable ranges
target_val = {}
target_range = {}
target_val['pas'] = {'soma': 120., 'trunk': 97., 'tuft': 121.}
target_range['pas'] = {'soma': 1., 'trunk': 1., 'tuft': 1.}
target_val['h'] = {'v_rest': v_init, 'soma': 75., 'trunk': 39., 'tuft': 42., 'v_rest_offset': 0., 'delta_EPSP_dur': 0.}
target_range['h'] = {'v_rest': 0.25, 'soma': 0.5, 'trunk': 1., 'tuft': 2., 'v_rest_offset': 0.1, 'delta_EPSP_dur': 0.05}
target_val['v_rest'] = {'soma': v_init, 'tuft_offset': 0.}
target_range['v_rest'] = {'soma': 0.25, 'tuft_offset': 0.1}
target_val['na_ka'] = {'v_rest': v_init, 'th_v': -51., 'soma_peak': 40., 'trunk_amp': 0.6, 'ADP': 0., 'AHP': 4.,
                       'stability': 0., 'ais_delay': 0., 'slow_depo': 25.}
target_range['na_ka'] = {'v_rest': 0.25, 'th_v': .2, 'soma_peak': 2., 'trunk_amp': 0.01, 'ADP': 0.01, 'AHP': .2,
                         'stability': 1., 'ais_delay': 0.001, 'slow_depo': 1.}


x0 = {}

#x1 = {}
xmin = {}
xmax = {}

# x0['pas'] = [1.43E-06, 1.58E-05, 121.9]
# x0['pas'] = [1.46E-06, 1.56E-05, 124.2]  # Error: 9.5227E+00
x0['pas'] = [1.52E-06, 1.63E-05, 121.9]  # Error: 2.6862E+00
xmin['pas'] = [1.0E-6, 1.0E-06, 25.]
xmax['pas'] = [2.0E-5, 2.0E-05, 200.]

# x0['h'] = [5.93E-08, 4.88E-03, 114.4, 459.0]
# x0['h'] = [2.27E-09, 2.14E-03, 265.3, 484.0]  # Error: 4.0003E+02
# x0['h'] = [2.35E-09, 1.76E-03, 381.8, 280.2]  # Error: 4.0260E+02
x0['h'] = [2.43E-09, 1.85E-03, 383.1, 260.7]  # Error: 3.5202E+02
xmin['h'] = [1.0E-9, 1.e-4, 25., 50.]
xmax['h'] = [1.0E-7, 1.e-2, 200., 500.]

# [soma.gkabar, soma.gkdrbar, axon.gkabar_kap factor, axon.gbar_nax factor]
# x0['na_ka_stability'] = [0.03555, 0.04975, 4.51, 2.18]
# x0['na_ka_stability'] = [0.0305, 0.0478, 2.97, 2.23]
x0['na_ka_stability'] = [0.0305, 0.0478, 2.97, 2.34]  # Error: 4.2416E+02
xmin['na_ka_stability'] = [0.01, 0.01, 1., 2.]
xmax['na_ka_stability'] = [0.05, 0.05, 5., 5.]

# [soma.sh_nas, trunk.ka factor]
# x0['na_ka_dend'] = [0., 2.9]
x0['na_ka_dend'] = [1.7, 2.8]
xmin['na_ka_dend'] = [0., 1.5]
xmax['na_ka_dend'] = [4., 5.]

# [ais.sha_nas, ais.gbar_nax factor]
x0['ais_delay'] = [-3.6, 4.4]

# [soma.e_pas, tuft.e_pas]
x0['v_rest'] = [-63., -77.]
xmin['v_rest'] = [v_init, soma_ek]
xmax['v_rest'] = [-63., -63.]

#x1 = dict(x0)
log = []

# maxfev = 200  # max number of iterations to run
maxfev = 1  # max number of iterations to run
take_step = Normalized_Step(x0['pas'], xmin['pas'], xmax['pas'])
minimizer_kwargs = dict(method=null_minimizer)


c = Client()
dv = c[:]
dv.clear()
dv.block = True
global_start_time = time.time()
dv.execute('from parallel_optimize_dendritic_excitability_engine import *')
#time.sleep(240)

v = c.load_balanced_view()

result = optimize.minimize(pas_error, x0['pas'], method='Nelder-Mead', options={'ftol': 1e-3,
                                                    'xtol': 1e-3, 'disp': True, 'maxiter': maxfev})
print result

"""
polished_result = optimize.basinhopping(pas_error, x0['pas'], niter=maxfev, niter_success=1,
                               disp=True, interval=20, minimizer_kwargs=minimizer_kwargs, take_step=take_step)

print polished_result
"""

#current_tuft_e_pas = x1['v_rest'][1]
#update_v_rest(x1['v_rest'])

# parallel_optimize_dendritic_excitability_engine.update_pas(result)

#update_pas(x1['pas'])
#update_h(x1['h'])

#current_ais_nax_factor = x1['ais_delay'][1]
#current_ka_trunk_factor = x1['na_ka_dend'][1]
#original_soma_kap = x1['na_ka_stability'][0]

#update_na_ka_stability(x1['na_ka_stability'])
#update_na_ka_dend(x1['na_ka_dend'])
#update_ais_delay(x1['ais_delay'])