__author__ = 'milsteina'
from function_lib import *
import parallel_optimize_BTSP
from ipyparallel import Client
from IPython.display import clear_output
from ipyparallel import interactive

"""
These methods attempt to fit the experimental BTSP data with a simple 3-state markov-like process model with 
nonstationary rates dependent on the presence of two intracellular biochemical signals.

Assumptions:
1) Synaptic weights are all = 1 prior to induction 1.
"""

try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass


if len(sys.argv) > 1:
    cluster_id = sys.argv[1]
    c = Client(cluster_id=cluster_id)
else:
    cluster_id = None
    c = Client()

if len(sys.argv) > 2:
    sub_controller_id = int(sys.argv[2])
else:
    sub_controller_id = 0

if len(sys.argv) > 3:
    group_size = int(sys.argv[3])
else:
    group_size = 1


experimental_file_dir = data_dir
experimental_filenames = {'cell': '121216 magee lab first induction', 'spont_cell': '120216 magee lab spont'}

dv = c[sub_controller_id+1:sub_controller_id+group_size]
dv.block = True

cell_ids = []
labels = []

for label, experimental_filename in experimental_filenames.iteritems():
    with h5py.File(experimental_file_dir+experimental_filename+'.hdf5') as f:
        cell_ids.extend(f.keys())
        labels.extend([label for i in range(len(f))])

for i, id in enumerate(dv.targets):
    c[id].execute('run parallel_optimize_BTSP_kinetic_rule_engine_052017 %s %i %i %s %s' %
                  (cluster_id, sub_controller_id, group_size, cell_ids[i], labels[i]), block=True)


@interactive
def ramp_population_error(x, induction=None):
    """
    Push the specified parameters to the engines. Each engine processes a different cell from an experimental dataset.
    Calculate the mean error across engines.
    :param x: array [local_decay_tau, global_decay_tau, local_kon, global_kon, global_koff, saturated_delta_weights]
    :param xmin: array
    :param xmax: array
    :return: float
    """
    formatted_x = '[' + ', '.join(['%.3E' % xi for xi in x]) + ']'
    if induction is None:
        induction = 1
    start_time = time.time()
    result = dv.apply_async(parallel_optimize_BTSP.report, x, induction)
    # result = dv.apply_async(parallel_optimize_BTSP_kinetic_rule_engine_052017.ramp_error_parametric, x, induction)
    last_buffer_len = []
    ready = False
    while not ready:
        try:
            ready = result.ready()
            time.sleep(1.)
            clear_output()
            for i, stdout in enumerate(result.stdout):
                if (i + 1) > len(last_buffer_len):
                    last_buffer_len.append(0)
                if stdout:
                    lines = stdout.splitlines()
                    if len(lines) > last_buffer_len[i]:
                        for line in lines[last_buffer_len[i]:]:
                            print line
                        last_buffer_len[i] = len(lines)
            sys.stdout.flush()
        except:
            ready = True
    result = result.get()
    Err = np.mean(result)
    # print result
    print 'Sub_controller %i: Mean error: %.4E; calculation of error across %i cells took %i s' % \
          (sub_controller_id, Err, len(dv), time.time() - start_time)
    sys.stdout.flush()
    return Err


@interactive
def report(x, induction):
    result = dv.apply_sync(parallel_optimize_BTSP.report, x, induction)
    return sub_controller_id, os.getpid(), result