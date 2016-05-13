__author__ = 'milsteina'
from specify_cells import *
import os
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which section to simulate (full sampling of sections).
"""

morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())


def test_single_section(sec_index):
    """
    :param sec_index: int
    :return: str
    """
    start_time = time.time()
    node = nodes[sec_index]
    loc = 0. if node.type == 'soma' else 0.5
    sim.modify_rec(0, node=node, loc=loc)
    sim.modify_stim(0, node=node, loc=loc)
    sim.run(v_init)
    with h5py.File(data_dir+rec_filename+'.hdf5', 'a') as f:
        sim.export_to_file(f, sec_index)
    print 'Process:', os.getpid(), 'completed Iteration:', sec_index, 'Node:', node.name, 'in', \
        time.time() - start_time, 's'
    return rec_filename


equilibrate = 250.  # time to steady-state
stim_dur = 400.
duration = equilibrate + stim_dur + 100.
amp = -0.15
v_init = -67.

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.zero_na()

nodes = cell.soma+cell.basal+cell.trunk+cell.apical+cell.tuft

sim = QuickSim(duration, verbose=False)
sim.append_rec(cell, cell.tree.root, 0.)
sim.append_stim(cell, cell.tree.root, 0., amp, equilibrate, stim_dur)