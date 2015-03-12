__author__ = 'milsteina'
from specify_cells import *
import os
"""
Parallel version: split the sections into pieces, and run 1 piece in each of many processes (4 cores).
Iterate through every section, activate AMPA_KIN synapses to measure EPSP attenuation.

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
Works best with hyperthreading off (a feature of the machine, can be disabled through XCode Instruments on Mac)
"""
num_cores = 4
new_rec_filename = '031115 kap_kad_ih_ampar_scale kd no_na - EB1 - epsp_attenuation sample'


def split_simulation(piece_index):
    """

    :param piece_index: int
    :return: str
    """
    import random

    morph_filename = 'Erik_Bloss_CA1_0215_Stitched_Proofread.swc'  # EB1
    #morph_filename = 'EB022715-stitched-proofread.swc'  # EB2
    #mech_filename = '030515 kap_kad_ih_ampar_scale kd with_na.pkl'
    mech_filename = '031015 kap_kad_ih_ampar_scale kd no_na.pkl'
    rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())

    equilibrate = 150.  # time to steady-state
    duration = 200.
    v_init = -65.
    syn_type = 'AMPA_KIN'

    syn_list = []
    cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
    """
    for node in cell.get_nodes_of_subtype('spine_head'):
        syn = Synapse(cell, node, ['AMPA_KIN'], stochastic=0)
        syn_list.append(syn)
    """
    random.seed(num_cores)
    for branch in cell.basal+cell.trunk+cell.apical+cell.tuft:
        if len(branch.spines) > 1:
            if branch.sec.L <= 10.:
                node = branch.spines[random.sample(range(0, len(branch.spines)), 1)[0]]
                syn = Synapse(cell, node, [syn_type], stochastic=0)
                syn_list.append(syn)
            else:
                num_syns = min(len(branch.spines), int(branch.sec.L//10.))  # a random synapse every 10 um
                for i in random.sample(range(0, len(branch.spines)), num_syns):
                    node = branch.spines[i]
                    syn = Synapse(cell, node, [syn_type], stochastic=0)
                    syn_list.append(syn)
        elif branch.spines:
            node = branch.spines[0]
            syn = Synapse(cell, node, [syn_type], stochastic=0)
            syn_list.append(syn)
    cell.init_synaptic_mechanisms()

    piece_size = len(syn_list) / num_cores
    #piece_size = 40 / num_cores
    pieces = []
    start = 0
    for i in range(1, num_cores):  # the first n-1 pieces contain equal number of items
        pieces.append(syn_list[start:i*piece_size])
        start += piece_size
    pieces.append(syn_list[start:])  # the last piece contains the rest of the items
    #pieces.append(syn_list[start:40])  # the last piece contains the rest of the items

    sim = QuickSim(duration)
    sim.parameters['equilibrate'] = equilibrate
    sim.parameters['duration'] = duration
    sim.append_rec(cell, cell.tree.root, description='soma')
    trunk = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type == 'trunk' and
                                               trunk.children[1].type == 'trunk'][0]  # trunk bifurcation
    sim.append_rec(cell, trunk, description='trunk')
    sim.append_rec(cell, trunk, description='branch')  # placeholders for branch and spine
    sim.append_rec(cell, trunk, description='spine')

    spike_times = h.Vector([equilibrate])

    start_time = time.time()

    simiter = 0
    for syn in pieces[piece_index]:
        spine = syn.node
        branch = spine.parent.parent
        sim.modify_rec(2, branch)
        sim.parameters['input_loc'] = branch.type
        sim.modify_rec(3, spine)
        syn.source.play(spike_times)
        print 'Run:', simiter, 'on Process:', os.getpid(), ', Section:', branch.name, ', Spine:', spine.name
        sim.run(v_init)
        with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
            sim.export_to_file(f, simiter)
        syn.source.play(h.Vector())  # playing an empty vector turns this synapse off for future runs while keeping the
                                 # VecStim source object in existence so it can be activated again
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
    v.execute('from parallel_EPSP_attenuation import *')
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
    print 'Parallel execution took:', time.time()-start_time, 's'
    print rec_file_list
    combine_output_files(rec_file_list, new_rec_filename)
    plot_EPSP_attenuation(new_rec_filename)
    plot_EPSP_kinetics(new_rec_filename)