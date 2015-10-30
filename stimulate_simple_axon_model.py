__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import sys
"""

Sweeps through hyperpolarizing and depolarizing current injections and saves to file.
Once files are generated from various mechanism dictionaries, the output can be compared and distance-dependent
propagation of the fast and slow components can be quantified and visualized.

"""

#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'


def zero_na():
    """

    """
    for sec_type in ['axon_hill', 'ais', 'soma']:
        cell.modify_mech_param(sec_type, 'nax', 'gbar', 0.)
    for sec_type in ['axon', 'basal', 'trunk', 'apical', 'tuft']:
        cell.reinitialize_subset_mechanisms(sec_type, 'nax')


def update_na_ka(soma_na_gbar, x):
    """

    :param x: array [soma.sh_nax, soma.gkabar, soma.gkdrbar, trunk.ka factor]
    """
    cell.modify_mech_param('soma', 'nax', 'gbar', soma_na_gbar)
    cell.modify_mech_param('soma', 'kdr', 'gkdrbar', x[2])
    cell.modify_mech_param('soma', 'kap', 'gkabar', x[1])
    cell.modify_mech_param('soma', 'nax', 'sh', x[0])
    slope = (x[3] - 1.) * x[1] / 300.
    for sec_type in ['basal', 'apical', 'trunk', 'tuft']:
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', min_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', max_loc=75., slope=slope, replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', max_loc=75., value=0.)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=75., max_loc=300., slope=slope,
                               value=x[1]+slope*75., replace=False)
        cell.modify_mech_param(sec_type, 'kad', 'gkabar', origin='soma', min_loc=300., value=x[1]+slope*300.,
                               replace=False)
        cell.modify_mech_param(sec_type, 'nax', 'gbar', origin='soma')
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='soma')
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param('axon_hill', 'nax', 'gbar', soma_na_gbar * 2.)
    cell.modify_mech_param('ais', 'nax', 'gbar', soma_na_gbar * 2 * 5.)
    cell.modify_mech_param('axon', 'nax', 'gbar', origin='axon_hill')
    for sec_type in ['axon_hill', 'ais', 'axon']:
        cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
        cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma')
        cell.modify_mech_param(sec_type, 'nax', 'sh', origin='soma')


def offset_vm():
    """

    :param sec_type: str
    """
    sim.modify_stim(0, amp=0.)
    sim.tstop = equilibrate
    i_holding = 0.
    t = np.arange(0., equilibrate, dt)
    offset = True
    while offset:
        sim.modify_stim(1, node=cell.tree.root, loc=0., amp=i_holding)
        sim.run(v_init)
        rec = sim.rec_list[0]['vec']
        vm = np.interp(t, sim.tvec, rec)
        v_rest = np.mean(vm[int((equilibrate - 3.)/dt):int((equilibrate - 1.)/dt)])
        if v_rest < v_init - 0.25:
            i_holding += 0.005
            if sim.verbose:
                print 'increasing i_holding to %.3f' % (i_holding)
        elif v_rest > v_init + 0.25:
            i_holding -= 0.005
            if sim.verbose:
                print 'decreasing i_holding to %.3f' % (i_holding)
        else:
            offset = False
    sim.tstop = duration


def get_plateau(vm):
    """

    :param t:
    :param vm:
    :return:
    """
    left = int((equilibrate-3.) / dt)
    right = int((equilibrate-1.) / dt)
    start = int((equilibrate+stim_dur-11.) / dt)
    end = int((equilibrate+stim_dur-1.) / dt)
    baseline = np.mean(vm[left:right])
    plateau = np.min(vm[start:end]) - baseline
    return plateau

def stim_sweep(f, vm_amp_targets=[-5., 15., 35.], step_sizes=[0.01, 0.01, 0.01]):
    """

    :param f: :class:'h5py.File'
    """
    sim.tstop = duration
    simiter = 0
    amp = 0.
    t = np.arange(0., duration, dt)
    last_step_size = step_sizes[0]
    for i, target in enumerate(vm_amp_targets):
        matched = False
        if not step_sizes[i] == last_step_size:
            last_step_size = step_sizes[i]
            amp = step_sizes[i]
        else:
            amp += step_sizes[i]
        direction = np.sign(target)
        while not matched:
            start_time = time.time()
            sim.modify_stim(0, amp=amp * direction)
            sim.run(v_init)
            vm = np.interp(t, sim.tvec, sim.rec_list[0]['vec'])
            plateau = get_plateau(vm)
            if plateau*direction >= target*direction:
                sim.parameters['vm_amp_target'] = target
                sim.parameters['amp'] = amp * direction
                sim.parameters['plateau'] = plateau
                sim.export_to_file(f, simiter)
                print 'Simulation took %i s with amp %.2f' % (time.time()-start_time, amp*direction)
                matched = True
            else:
                amp += step_sizes[i]
                print 'Changing amp to %.2f' % (amp*direction)
        simiter += 1


equilibrate = 250.  # time to steady-state
stim_dur = 100.
duration = equilibrate + stim_dur + 100.
dt = 0.02
th_dvdt = 20.
v_init = -65.

#mech_filename = '101815 simple_axon_model_uniform_km2'
#mech_filename = '102615 simple_axon_model_tuning_thinner_axon'
#mech_filename = '102615 simple_axon_model_tuning_thinner_axon_no_na'
#mech_filename = '102615 simple_axon_model_tuning_thinner_axon_no_na_reduced_kap_kdr_km2'

#mech_filename = '102715 simple_axon_model_no_na_reduced_k'
#mech_filename = '102715 simple_axon_model_no_na'
mech_filename = '102715 simple_axon_model'

#mech_filename = '101815 simple_axon_model_uniform_km2_no_na'
#mech_filename = '101815 simple_axon_model_uniform_km2_no_na_reduced_kap_kdr_km2'

if len(sys.argv) > 1:
    mech_filename = str(sys.argv[1])

rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+mech_filename

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=False)

#axon_seg_locs = [seg.x for seg in cell.axon[2].sec]
axon_seg_locs = np.arange(30., 500., 30.) / 500.

sim = QuickSim(duration)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=equilibrate, dur=stim_dur)
sim.append_stim(cell, cell.tree.root, loc=0., amp=0., delay=0., dur=duration)
sim.append_rec(cell, cell.tree.root, loc=0., description='soma')
sim.append_rec(cell, cell.axon[1], loc=1., description='ais')
for loc in axon_seg_locs:
    sim.append_rec(cell, cell.axon[2], loc=loc, description='axon')
sim.parameters['duration'] = duration
sim.parameters['stim_dur'] = stim_dur
sim.parameters['equilibrate'] = equilibrate

offset_vm()

#[soma.sh_nax, soma.gkabar, soma.gkdrbar, trunk.ka factor]
#x0['na_ka'] = [2.7, 4.367E-02, 1.282E-02, 1.5]

#rec_filename = '102715 test simple_axon_model_no_na_reduced_k'
#rec_filename = '102715 test simple_axon_model_no_na'
rec_filename = '102715 test simple_axon_model'
#rec_filename = '102815 test simple_axon_model_spike_height'


with h5py.File(data_dir+rec_filename+'.hdf5', 'w') as f:
    stim_sweep(f)