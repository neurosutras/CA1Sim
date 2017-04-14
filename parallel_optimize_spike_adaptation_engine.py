__author__ = 'Grace Ng'
from specify_cells2 import *
import os
import sys
from ipyparallel import interactive
import mkl

"""
Builds a cell locally so each engine is ready to receive jobs one at a time, specified by a value for the amplitude of
a somatic current injection to test spike adaptation.
"""

mkl.set_num_threads(1)

neurotree_filename = '121516_DGC_trees.pkl'
neurotree_dict = read_from_pkl(morph_dir+neurotree_filename)

rec_filename = str(time.strftime('%m%d%Y', time.gmtime()))+'_'+str(time.strftime('%H%M%S', time.gmtime()))+\
               '_pid'+str(os.getpid())+'_sim_output'


# placeholder for optimization parameter, must be pushed to each engine on each iteration
# x: array [soma.gCa factor, soma.gCadepK factor, soma.gkmbar]
x = None  # Placeholder for parameters pushed from controller.

xmin = {}
xmax = {}

# [soma.gCa factor, soma.gCadepK factor, soma.gkmbar]
xmin['spike_adaptation'] = [0.5, 0.5, 0.0005]
xmax['spike_adaptation'] = [2., 2., 0.003]

check_bounds = CheckBounds(xmin, xmax)


if len(sys.argv) > 1:
    spines = bool(int(sys.argv[1]))
else:
    spines = False
if len(sys.argv) > 2:
    mech_filename = str(sys.argv[2])
else:
    mech_filename = '041417 GC optimizing excitability'


@interactive
def update_mech_dict():
    update_spike_adaptation(x)
    cell.export_mech_dict(cell.mech_filename)


@interactive
def update_spike_adaptation(x):
    """

    :param x: array [soma.gCa factor, soma.gCadepK factor, soma.gkmbar]
    """

    cell.modify_mech_param('soma', 'Ca', 'gcamult', x[0])
    cell.modify_mech_param('soma', 'CadepK', 'gcakmult', x[1])
    cell.modify_mech_param('soma', 'km3', 'gkmbar', x[2])
    for sec_type in ['axon_hill', 'ais', 'axon']:
        cell.reinitialize_subset_mechanisms(sec_type, 'km3')


@interactive
def offset_vm(description, vm_target=None):
    """

    :param description: str
    :param vm_target: float
    """
    if vm_target is None:
        vm_target = v_init
    sim.modify_stim(0, amp=0.)
    node = rec_nodes[description]
    loc = rec_locs[description]
    rec_dict = sim.get_rec(description)
    sim.modify_stim(1, node=node, loc=loc, amp=0.)
    rec = rec_dict['vec']
    offset = True
    sim.tstop = equilibrate
    t = np.arange(0., equilibrate, dt)
    sim.modify_stim(1, amp=i_holding[description])
    sim.run(vm_target)
    vm = np.interp(t, sim.tvec, rec)
    v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
    initial_v_rest = v_rest
    if v_rest < vm_target - 0.5:
        i_holding[description] += 0.01
        while offset:
            if sim.verbose:
                print 'increasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(1, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest < vm_target - 0.5:
                i_holding[description] += 0.01
            else:
                offset = False
    elif v_rest > vm_target + 0.5:
        i_holding[description] -= 0.01
        while offset:
            if sim.verbose:
                print 'decreasing i_holding to %.3f (%s)' % (i_holding[description], description)
            sim.modify_stim(1, amp=i_holding[description])
            sim.run(vm_target)
            vm = np.interp(t, sim.tvec, rec)
            v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
            if v_rest > vm_target + 0.5:
                i_holding[description] -= 0.01
            else:
                offset = False
    sim.tstop = duration
    return v_rest


@interactive
def get_rheobase(local_x=None, plot=False):
    """
    :param local_x: array [soma.gCa factor, soma.gCadepK factor, soma.gkmbar]
    :param plot: bool
    :return: float
    """
    if local_x is None:
        local_x = x
    if not check_bounds.within_bounds(local_x, 'spike_adaptation'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return None
    start_time = time.time()
    update_spike_adaptation(local_x)
    soma_vm = offset_vm('soma', v_active)
    print 'Process %i: Getting here - after offset_vm' % os.getpid()
    sim.modify_stim(0, dur=100.)
    duration = equilibrate + 100.
    sim.tstop = duration
    t = np.arange(0., duration, dt)
    spike = False
    d_amp = 0.01
    amp = max(0., i_th['soma'] - 0.05)
    while not spike:
        sim.modify_stim(0, amp=amp)
        sim.run(v_active)
        vm = np.interp(t, sim.tvec, sim.get_rec('soma')['vec'])
        if np.any(vm[:int(equilibrate/dt)] > -30.):
            print 'Process %i: Aborting - spontaneous firing' % (os.getpid())
            return None
        if np.any(vm[int(equilibrate/dt):int((equilibrate+50.)/dt)] > -30.):
            spike = True
        elif amp >= 0.4:
            print 'Process %i: Aborting - rheobase outside target range' % (os.getpid())
            return None
        else:
            amp += d_amp
            if sim.verbose:
                print 'increasing amp to %.3f' % amp
    i_th['soma'] = amp
    print 'Process %i took %i s to find spike rheobase at amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    if plot:
        sim.plot()
    return amp


@interactive
def sim_f_I(amp, local_x=None, plot=False):
    """
    
    :param amp: float 
    :param local_x: array
    :param plot: bool
    :return: dict
    """
    if local_x is None:
        local_x = x
    if not check_bounds.within_bounds(local_x, 'spike_adaptation'):
        print 'Process %i: Aborting - Parameters outside optimization bounds.' % (os.getpid())
        return None
    update_spike_adaptation(local_x)
    soma_vm = offset_vm('soma', v_active)
    print 'Process %i: Getting here - after offset_vm' % os.getpid()
    sim.parameters['amp'] = amp
    start_time = time.time()
    sim.modify_stim(0, dur=stim_dur, amp=amp)
    duration = equilibrate + stim_dur
    sim.tstop = duration
    sim.run(v_active)
    if plot:
        sim.plot()
    spike_times = np.subtract(cell.spike_detector.get_recordvec().to_python(), equilibrate)
    result = {}
    result['spike_times'] = spike_times
    result['amp'] = amp
    print 'Process %i took %i s to run simulation with I_inj amp: %.3f' % (os.getpid(), time.time() - start_time, amp)
    return result


@interactive
def export_sim_results():
    """
    Export the most recent time and recorded waveforms from the QuickSim object.
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
        sim.export_to_file(f)


equilibrate = 250.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
dt = 0.02
amp = 0.3
th_dvdt = 10.
v_init = -77.
v_active = -77.

cell = DG_GC(neurotree_dict=neurotree_dict[0], mech_filename=mech_filename, full_spines=spines)
if spines is False:
    cell.correct_for_spines()

cell.modify_mech_param('soma', 'Ca', 'gcamult', 1.)
cell.modify_mech_param('soma', 'CadepK', 'gcakmult', 1.)

rec_locs = {'soma': 0., 'axon': 1.}
rec_nodes = {'soma': cell.tree.root, 'axon': cell.axon[2]}

# sim = QuickSim(duration, verbose=False)
sim = QuickSim(duration, verbose=True)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)

for description, node in rec_nodes.iteritems():
    sim.append_rec(cell, node, loc=rec_locs[description], description=description)

i_holding = {'soma': 0.}
i_th = {'soma': 0.05}

spike_output_vec = h.Vector()
cell.spike_detector.record(spike_output_vec)

"""
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ina_nas',
               description='Soma nas_i')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ik_kap',
               description='Soma kap_i')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ik_kdr',
               description='Soma kdr_i')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ik_km3',
               description='Soma km_i')
sim.append_rec(cell, cell.axon[2], loc=0.5, description='Axon Vm')

sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_i_Ca',
               description='Soma Ca_i')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_ca_i_Ca',
               description='Soma Ca_intra_Ca')
sim.append_rec(cell, cell.tree.root, loc=0.5, object=cell.tree.root.sec(0.5), param='_ref_i_CadepK',
               description='Soma CadepK_i')
"""