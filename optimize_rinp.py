__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
"""
Optimizes input resistance by injecting hyperpolarizing current at soma, trunk, and tuft. Modifies g_pas to fit target
values from Magee 1995, Bittner et al., 2012 in the absence of Ih.
"""
morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
#mech_filename = '032315 kap_kad_ampar_scale nmda kd pas no_ih no_na.pkl'
mech_filename = '032315 kap_kad_ih_ampar_scale nmda kd pas no_na.pkl'
rec_filename = 'calibrate_rinp'


def plot_rinp_sec_types():
    """

    :param x: array
    """
    f = h5py.File(data_dir+rec_filename+'.hdf5', 'w')
    rinps = []
    for i, (sec_type, node) in enumerate(zip(['soma', 'trunk', 'tuft'], [cell.tree.root, trunk, tuft])):
        sim.modify_rec(0, node=node, description=sec_type)
        sim.modify_stim(0, node=node)
        sim.run(v_init)
        sim.export_to_file(f, i)
        rinps.append(get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[2])
    f.close()
    print('R_Inp: soma: %i, trunk: %i, tuft: %i' % (rinps[0], rinps[1], rinps[2]))
    plot_superimpose_conditions(rec_filename)


def rinp_error(x, plot=0):
    """
    :param x: array
    :param plot: int
    :return: float
    """
    start_time = time.time()
    cell.modify_mech_param('soma', 'pas', 'g', x[0])
    for sec_type in ['axon', 'ais', 'spine_neck', 'spine_head']:
        cell.reinitialize_subset_mechanisms(sec_type, 'pas')
    cell.modify_mech_param('basal', 'pas', 'g', x[1])
    cell.modify_mech_param('trunk', 'pas', 'g', x[1])
    cell.reinitialize_subset_mechanisms('apical', 'pas')
    cell.modify_mech_param('tuft', 'pas', 'g', x[2])
    result = {}
    if plot:
        f = h5py.File(data_dir+rec_filename+'.hdf5', 'w')
        i = 0
    for sec_type, node in zip(['soma', 'trunk', 'tuft'], [cell.tree.root, trunk, tuft]):
        sim.modify_rec(0, node=node, description=sec_type)
        sim.modify_stim(0, node=node)
        sim.run(v_init)
        if plot:
            sim.export_to_file(f, i)
            i += 1
        result[sec_type] = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[2]
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print('Current x: [%.4E, %.4E, %.4E], took %i s' % (x[0], x[1], x[2], time.time()-start_time))
    print('Error: %.4E, R_Inp: soma: %i, trunk: %i, tuft: %i' % (Err, result['soma'], result['trunk'], result['tuft']))
    if plot:
        f.close()
        plot_superimpose_conditions(rec_filename)
    else:
        return Err


equilibrate = 200.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
amp = -0.1
v_init = -65.

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

sim = QuickSim(duration, verbose=False)
sim.parameters['description'] = 'RInp (Steady-state)'
sim.append_stim(cell, cell.tree.root, 1., amp, equilibrate, stim_dur)
sim.append_rec(cell, cell.tree.root, loc=1.)
# this selection of specific branches requires manual selection for different morphologies
trunk = cell.tree.get_node_with_index(46)  # ~50 um from SLM along thicker of two trunk branches (EB2)
tuft = cell.tree.get_node_with_index(52)  # ~50 um into SLM along thick primary tuft branch (EB2)

#the target values and acceptable ranges
target_val = {'soma': 110., 'trunk': 100., 'tuft': 120.}
target_range = {'soma': 5., 'trunk': 5., 'tuft': 5.}

#the initial guess and bounds
# x = [soma.g_pas, trunk.g_pas slope, tuft.g_pas slope]
#x0 = [1e-6, 5e-5, 1e-4]
xmin = [1e-7, 1e-6, 1e-5]  # first-pass bounds
xmax = [1e-4, 1e-4, 1e-3]
x0 = [8.95E-06, 3.04E-05, 3.85E-04]  # following first pass basinhopping with steps taken in log space

# rewrite the bounds in the way required by optimize.minimize
xbounds = [(low, high) for low, high in zip(xmin, xmax)]

blocksize = 0.5  # defines the fraction of the xrange that will be explored at each step
                 #  basinhopping starts with this value and reduces it by 10% every 'interval' iterations

mytakestep = MyTakeStep(blocksize, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)
"""
result = optimize.basinhopping(rinp_error, x0, niter= 720, niter_success=100, disp=True, interval=20,
                                                    minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
new_x = result.x
print('g_pas: soma, trunk, tuft: [%.2E, %.2E, %.2E]' % (new_x[0], new_x[1], new_x[2]))
"""
#rinp_error(new_x, plot=1)
rinp_error(x0, plot=1)