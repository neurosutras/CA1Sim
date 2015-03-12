__author__ = 'milsteina'
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which section to simulate (full sampling of sections).
"""
from specify_cells import *
import os

#morph_filename = 'Erik_Bloss_CA1_0215_Stitched_Proofread.swc'  # EB1
morph_filename = 'EB022715-stitched-proofread.swc'  # EB2
#mech_filename = '030515 kap_kad_ih_ampar_scale kd with_na.pkl'
mech_filename = '031215 kap_kad_ih_ampar_scale kd no_na.pkl'
rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())


def test_single_section(sec_index):
    """
    :param sec_index: int
    :return: str
    """
    start_time = time.time()
    node = nodes[sec_index]
    sim.modify_rec(0, node=node)
    sim.modify_stim(0, node=node)
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, sec_index)
    print 'Process:', os.getpid(), 'completed Iteration:', sec_index, 'Node:', node.name, 'in', \
        time.time() - start_time, 's'
    return rec_filename


equilibrate = 150.  # time to steady-state
duration = 250.
stim_dur = 100.
amp = -0.1
v_init = -65.

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

nodes = cell.basal+cell.trunk+cell.apical+cell.tuft

sim = QuickSim(duration, verbose=False)
sim.append_rec(cell, cell.tree.root, 0.5)
sim.append_stim(cell, cell.tree.root, 0.5, amp, equilibrate, stim_dur)