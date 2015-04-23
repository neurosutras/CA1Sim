__author__ = 'milsteina'
from specify_cells import *
import os

"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which section to simulate (full sampling of sections).
"""

#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'

#mech_filename = '042215 soma_pas - EB2.pkl'
#mech_filename = '042215 soma_pas spines - EB2.pkl'
#mech_filename = '042215 soma_pas kdr ka_scale - EB2.pkl'
#mech_filename = '042215 soma_pas kdr ka_scale - adjusted - EB2.pkl'
#mech_filename = '042215 pas_exp_scale kdr ka_scale - EB2.pkl'
mech_filename = '042315 pas_ka_scale kdr - EB2.pkl'

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


equilibrate = 200.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
amp = -0.15
v_init = -80.

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
#cell = CA1_Pyr(morph_filename, mech_filename, full_spines=False)

nodes = cell.soma+cell.basal+cell.trunk+cell.apical+cell.tuft

sim = QuickSim(duration, verbose=False)
sim.append_rec(cell, cell.tree.root, 0.5)
sim.append_stim(cell, cell.tree.root, 0.5, amp, equilibrate, stim_dur)