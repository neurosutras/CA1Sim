__author__ = 'milsteina'
from specify_cells import *
import os
"""
Parallel version: split the sections into 4 groups, and run 1 group on each of 4 cores.
Iterate through every section, inject hyperpolarizing current to measure input resistance.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
Works best with hyperthreading off (a feature of the machine, can be disabled through XCode Instruments on Mac)
"""
num_cores = 4


def split_simulation(piece_index):
    """

    :param piece_index: int
    :return: str
    """
    morph_filename = 'Erik_Bloss_CA1_0215_Stitched_Proofread.swc'  # EB1
    #morph_filename = 'EB022715-stitched-proofread.swc'  # EB2
    mech_filename = '030515 kap_kad_ih_ampar_scale kd with_na.pkl'
    new_rec_filename = '030515 kap_kad_ih_scale kd na full_spines - EB2 - rinp'
    rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())

    #cell = HocCell(morph_filename, mech_filename)
    cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

    nodes = cell.basal+cell.trunk+cell.apical+cell.tuft
    piece_size = len(nodes) / num_cores
    pieces = []
    start = 0
    simiter = 0
    for i in range(1, num_cores):  # the first n-1 pieces contain equal number of items
        pieces.append(nodes[start:i*piece_size])
        start += piece_size
    pieces.append(nodes[start:])  # the last piece contains the rest of the items

    equilibrate = 150.  # time to steady-state
    duration = 250.
    stim_dur = 100.
    amp = -0.1
    v_init = -65.
    sim = QuickSim(duration)
    sim.append_rec(cell, cell.tree.root, 0.5)
    sim.append_stim(cell, cell.tree.root, 0.5, amp, equilibrate, stim_dur)

    for node in pieces[piece_index]:
        sim.modify_rec(0, node)
        sim.modify_stim(node=node)
        print 'Run:', simiter, 'on Process:', os.getpid(), ', Section:', node.name
        sim.run(v_init)
        with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
            sim.export_to_file(f, simiter)
        simiter += 1
    print 'Process:', os.getpid(), 'completed', simiter, 'iterations in', time.time() - start_time, 's'
    return rec_filename

if __name__ == '__main__':
    from IPython.parallel import Client
    from IPython.display import clear_output
    from plot_results import *
    import sys

    c = Client()
    v = c[:]
    start_time = time.time()
    v.execute('from parallel_rinp import *')
    result = v.map_async(split_simulation, range(num_cores))
    while not result.ready():
        clear_output()
        for stdout in result.stdout:
            if stdout:
                lines = stdout.split('\n')
                for line in lines[-4:-1]:
                    if line:
                        print line
        sys.stdout.flush()
        time.sleep(60)
    rec_file_list = result.get()
    print 'Parallel execution took: ', time.time()-start_time, ' s'
    #new_rec_filename = '030515 kap_kad_ih_scale kd na full_spines - EB1 - rinp'
    #combine_output_files(rec_file_list, new_rec_filename)
    #plot_Rinp(new_rec_filename)