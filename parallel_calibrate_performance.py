__author__ = 'milsteina'
from specify_cells import *
import os
"""
Runs 4 iterations on each node of a simple simulation with a step current injection into the soma of a single cell.
Allows calibration of performance scaling with number of processes/cores/hyperthreads.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
Hyperthreading  is a feature of the machine, can be disabled through XCode Instruments on Mac.
"""


def split_simulation(piece_index):
    """

    :param piece_index: int
    :return: str
    """
    #morph_filename = 'EB1-early-bifurcation.swc'
    morph_filename = 'EB2-late-bifurcation.swc'
    mech_filename = '030515 kap_kad_ih_ampar_scale kd with_na.pkl'
    rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())

    cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

    nodes = cell.basal+cell.trunk+cell.apical+cell.tuft
    piece_size = 4
    start = piece_size * piece_index
    stop = start + 4

    equilibrate = 150.  # time to steady-state
    duration = 250.
    stim_dur = 100.
    amp = -0.1
    v_init = -65.
    sim = QuickSim(duration, verbose=False)
    sim.append_rec(cell, cell.tree.root, 0.5)
    sim.append_stim(cell, cell.tree.root, 0.5, amp, equilibrate, stim_dur)
    
    start_time = time.time()
    simiter = 0
    for node in nodes[start:stop]:
        iteration_time = time.time()
        sim.modify_rec(0, node=node)
        sim.modify_stim(0, node=node)
        sim.run(v_init)
        print 'Process:', os.getpid(), 'Iteration:', simiter, 'Section:', node.name, 'took', time.time() - \
                                                                            iteration_time, 's'
        simiter += 1
    print 'Process:', os.getpid(), 'completed', simiter, 'iterations in', time.time() - start_time, 's'
    return simiter

if __name__ == '__main__':
    from IPython.parallel import Client
    from IPython.display import clear_output
    from plot_results import *
    import sys

    c = Client()
    v = c[:]
    start_time = time.time()
    v.execute('from parallel_calibrate_performance import *')
    result = v.map_async(split_simulation, range(len(c)))
    while not result.ready():
        clear_output()
        for stdout in result.stdout:
            if stdout:
                lines = stdout.split('\n')
                if lines[-2]:
                    print lines[-2]
        sys.stdout.flush()
        time.sleep(5)
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    print 'Parallel execution took: ', time.time()-start_time, ' s'