__author__ = 'milsteina'
from specify_cells import *
import os
"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by an index corresponding to
which section to simulate (full sampling of sections).
"""

morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())+'_bAP'


def offset_vm(sec_type, vm_target=None):
    """

    :param sec_type: str
    :param vm_target: float
    """
    if vm_target is None:
        vm_target = v_init
    sim.modify_stim(0, amp=0.)
    if sec_type == 'soma':
        sim.modify_stim(1, node=cell.tree.root, loc=0., amp=0.)
        rec = sim.rec_list[0]['vec']
    """
    elif sec_type == 'trunk':
        sim.modify_stim(1, node=trunk, loc=1., amp=0.)
        rec = sim.rec_list[1]['vec']
    elif sec_type == 'tuft':
        sim.modify_stim(1, node=tuft, loc=1., amp=0.)
        rec = sim.rec_list[2]['vec']
    """
    offset = True
    sim.tstop = equilibrate
    sim.run(vm_target)
    t = np.arange(0., equilibrate, dt)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    sim.modify_stim(1, amp=i_holding[sec_type])
    sim.run(vm_target)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    if v_rest < vm_target - 0.5:
        i_holding[sec_type] += 0.01
        while offset:
            if sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (i_holding[sec_type], sec_type)
            sim.modify_stim(1, amp=i_holding[sec_type])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest < vm_target - 0.5:
                i_holding[sec_type] += 0.01
            else:
                offset = False
    elif v_rest > vm_target + 0.5:
        i_holding[sec_type] -= 0.01
        while offset:
            if sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (i_holding[sec_type], sec_type)
            sim.modify_stim(1, amp=i_holding[sec_type])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                i_holding[sec_type] -= 0.01
            else:
                offset = False
    sim.tstop = duration
    return initial_v_rest

def find_spike_threshold():
    """

    """
    t = np.arange(0., duration, dt)
    spike = False
    amp = 0.05
    while not spike:
        sim.modify_stim(0, amp=amp)
        if sim.verbose:
            print 'increasing amp to %.3f' % amp
        sim.run(v_init)
        vm = np.interp(t, sim.tvec, sim.rec_list[0]['vec'])
        if np.any(vm[int(equilibrate/dt):int((equilibrate+50.)/dt)] > -30.):
            spike = True
        else:
            amp += 0.05


equilibrate = 250.  # time to steady-state
stim_dur = 100.
duration = equilibrate + stim_dur + 100.
amp = 0.1
v_active = -61.
v_init = -67.
dt = 0.02

i_holding = {'soma': 0.12, 'trunk': 0., 'tuft': 0.}

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.set_terminal_branch_na_gradient()

nodes = cell.basal+cell.trunk+cell.apical+cell.tuft

sim = QuickSim(duration)
sim.append_rec(cell, cell.tree.root, 0.)
sim.append_stim(cell, cell.tree.root, 0., 0., equilibrate, stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)
offset_vm('soma', v_active)
find_spike_threshold()

for node in nodes:
    sim.append_rec(cell, node, 0.5)

sim.run(v_active)

with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
    sim.export_to_file(f, 0)
