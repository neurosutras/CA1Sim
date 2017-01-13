__author__ = 'milsteina'
import parallel_optimize_spine_cm_correction_engine
import sys
import os
from ipyparallel import Client
from IPython.display import clear_output
import scipy.optimize as optimize
from function_lib import *

"""
Parallel version: Iterates through every section, injecting hyperpolarizing current and measuring input resistance.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""

spines = False
if len(sys.argv) > 1:
    mech_filename = str(sys.argv[1])
else:
    mech_filename = '121516 DG_GC pas spines'
if len(sys.argv) > 2:
    cluster_id = sys.argv[2]
    c = Client(cluster_id=cluster_id)
else:
    c = Client()


dv = c[:]
dv.clear()
dv.block = True
dv.execute('run parallel_optimize_spine_cm_correction_engine \"%s\"' % mech_filename)
# time.sleep(120)
v = c.load_balanced_view()


class History(object):
    def __init__(self):
        """

        """
        self.x = []
        self.Err = []

    def report_best(self):
        """
        Report the input parameters and output values with the lowest error.
        """
        lowest_Err = min(self.Err)
        index = self.Err.index(lowest_Err)
        best_x = self.x[index]
        formatted_x = '[' + ', '.join(['%.2f' % xi for xi in best_x]) + ']'
        print 'best x: %s' % formatted_x
        print 'lowest Err: %.3E' % lowest_Err
        return best_x


hist = History()


def check_bounds(x):
    """
    Check that the current set of parameters for optimization are within the desired bounds.
    :param x: array
    :return: bool
    """
    for i in range(len(x)):
        if ((xmin[i] is not None and x[i] < xmin[i]) or
                (xmax[i] is not None and x[i] > xmax[i])):
            return False
    return True


def cm_correction_error(x):
    """

    :param x: array
    :return: float
    """
    if not check_bounds(x):
        print 'Aborting: Invalid parameter values.'
        return 1e9

    dv.map_sync(parallel_optimize_spine_cm_correction_engine.correct_for_spines_cell, [x[0] for i in range(len(c))])
    num_secs = len(parallel_optimize_spine_cm_correction_engine.nodes)
    result = v.map_async(parallel_optimize_spine_cm_correction_engine.test_single_section, range(num_secs))
    start_time = time.time()
    last = []
    while not result.ready():
        time.sleep(1.)
        clear_output()
        for i, stdout in enumerate([stdout for stdout in result.stdout if stdout][-num_secs:]):
            line = stdout.splitlines()[-1]
            if line not in last:
                print line
                last.append(line)
        if len(last) > num_secs:
            last = last[-num_secs:]
        sys.stdout.flush()
    print 'Parallel execution took: ', time.time()-start_time, ' s'
    rec_file_list = [filename for filename in dv['rec_filename'] if os.path.isfile(data_dir+filename+'.hdf5')]

    distances = []
    decay_90 = []
    for filename in rec_file_list:
        with h5py.File(data_dir + filename + '.hdf5', 'r') as r:
            for trial in r.itervalues():
                rec = trial['rec'].itervalues().next()
                distance = rec.attrs['soma_distance']
                distances.append(distance)
                this_decay_90 = trial.attrs['decay_90']
                decay_90.append(this_decay_90)
    for filename in rec_file_list:
        os.remove(data_dir+filename+'.hdf5')
    indexes = range(len(distances))
    indexes.sort(key=distances.__getitem__)
    decay_90 = map(decay_90.__getitem__, indexes)
    if len(decay_90) != len(target_decay_90):
        raise Exception('Mismatch in number of sections')
    Err = 0.
    for i in range(len(decay_90)):
        Err += ((target_decay_90[i] - decay_90[i])/sigma) ** 2.

    print 'cm_correction_fraction: %.2f; Error: %.4E' % (x[0], Err)
    hist.x.append(x)
    hist.Err.append(Err)
    return Err


with_spines_filename = '011120171404_Rinp_uniform_pas_spines'
target_distances = []
target_decay_90 = []
dend_sec_types = ['basal', 'trunk', 'apical', 'tuft']
with h5py.File(data_dir + with_spines_filename + '.hdf5', 'r') as f:
    for trial in f['Rinp_data'].itervalues():
        sec_type = trial.attrs['type']
        if sec_type in dend_sec_types:
            distance = trial.attrs['soma_distance']
            this_decay_90 = trial.attrs['decay_90']
            target_distances.append(distance)
            target_decay_90.append(this_decay_90)
indexes = range(len(target_distances))
indexes.sort(key=target_distances.__getitem__)
target_decay_90 = map(target_decay_90.__getitem__, indexes)

sigma = 1.

x0 = [0.5]
xmin = [0.]
xmax = [1.]

polish_niter = 200

result = optimize.minimize(cm_correction_error, x0, method='Nelder-Mead', options={'ftol': 1e-5,
                                                                                   'disp': True,
                                                                                   'maxiter': polish_niter})

# cm_correction_error(x0)