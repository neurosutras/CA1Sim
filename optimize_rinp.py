__author__ = 'milsteina'
from specify_cells import *
from plot_results import *
import scipy.optimize as optimize
"""
Optimizes input resistance by injecting hyperpolarizing current at soma, trunk, and tuft. Modifies g_pas to fit target
values from Magee 1995, Bittner et al., 2012 in the absence of Ih.
"""
morph_filename = 'EB022715-stitched-proofread.swc'  # EB2
mech_filename = '031815 calibrate nmda gmax.pkl'
rec_filename = 'calibrate_rinp'


def null_minimizer(fun, x0, args, **options):
    """
    Rather than allow basinhopping to pass each local mimimum to a gradient descent algorithm for polishing, this method
    just catches and passes all local minima so basinhopping can proceed.
    """
    return optimize.OptimizeResult(x=x0, fun=fun(x0, *args), success=True, nfev=1)


class MyTakeStep(object):
    """
    Converts basinhopping absolute stepsize into different stepsizes for each parameter such that the stepsizes are
    some fraction of the ranges specified by xmin and xmax. Also enforces bounds for x.
    """
    def __init__(self, blocksize, xmin, xmax, stepsize=0.5):
        self.stepsize = stepsize
        self.blocksize = blocksize
        self.xmin = xmin
        self.xmax = xmax
        self.xrange = []
        for i in range(len(self.xmin)):
            self.xrange.append(self.xmax[i] - self.xmin[i])

    def __call__(self, x):
        for i in range(len(x)):
            if x[i] < self.xmin[i]:
                x[i] = self.xmin[i]
            if x[i] > self.xmax[i]:
                x[i] = self.xmax[i]
            snew = self.stepsize / 0.5 * self.blocksize * self.xrange[i] / 2.
            sinc = min(self.xmax[i] - x[i], snew)
            sdec = min(x[i]-self.xmin[i], snew)
            x[i] += np.random.uniform(-sdec, sinc)
        return x


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
        rinps.append(get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[1])
    f.close()
    print('R_Inp: soma: %i, trunk: %i, tuft: %i' % (rinps[0], rinps[1], rinps[2]))
    plot_superimpose_conditions(rec_filename)



def rinp_error(x, sec_type):
    """
    :param x: array
    :param sec_type: str
    :return: float
    """
    start_time = time.time()
    if sec_type == 'soma':
        node = cell.tree.root
        cell.modify_mech_param('soma', 'pas', 'g', x[0])
        cell.modify_mech_param('basal', 'pas', 'g', origin='soma')
        cell.modify_mech_param('trunk', 'pas', 'g', origin='soma')
        cell.modify_mech_param('apical', 'pas', 'g', origin='trunk')
        cell.modify_mech_param('tuft', 'pas', 'g', origin='trunk')
    elif sec_type == 'trunk':
        node = trunk
        cell.modify_mech_param('trunk', 'pas', 'g', origin='soma', slope=x[0])
        cell.modify_mech_param('apical', 'pas', 'g', origin='trunk')
        cell.modify_mech_param('tuft', 'pas', 'g', origin='trunk')
    elif sec_type == 'tuft':
        node = tuft
        cell.modify_mech_param('tuft', 'pas', 'g', origin='trunk', slope=x[0])
    sim.modify_rec(0, node=node, description=sec_type)
    sim.modify_stim(0, node=node)
    sim.run(v_init)
    rinp = get_Rinp(np.array(sim.tvec), np.array(sim.rec_list[0]['vec']), equilibrate, duration, amp)[1]
    result = {sec_type: rinp}
    Err = 0.
    for target in result:
        Err += ((target_val[target] - result[target])/target_range[target])**2.
    print('Error: %.4f, sec_type: %s , R_Inp: %i, current x: %.4E, took %i s' % (Err, sec_type, rinp, x[0],
                                                                                 time.time()-start_time))
    return Err


equilibrate = 200.  # time to steady-state
stim_dur = 500.
duration = equilibrate + stim_dur
amp = -0.1
v_init = -65.

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

sim = QuickSim(duration, verbose=False)
sim.parameters['description'] = 'RInp (Steady-state)'
sim.append_stim(cell, cell.tree.root, 0.5, amp, equilibrate, stim_dur)
sim.append_rec(cell, cell.tree.root)
trunk = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type == 'trunk' and
                                               trunk.children[1].type == 'trunk'][0]  # trunk bifurcation
# 50 um into tuft
tuft = [tuft for tuft in cell.tuft if cell.get_distance_to_node(cell.get_dendrite_origin(tuft), tuft, 0.5) >= 50.][0]

#the target values and acceptable ranges
target_val = {'soma': 110., 'trunk': 100., 'tuft': 120.}
target_range = {'soma': 5., 'trunk': 5., 'tuft': 5.}

#the initial guess and bounds
# x = [soma.g_pas, trunk.g_pas slope, tuft.g_pas slope]
x0 = [1.5e-5, 1e-7, 1e-7]
xmin = [1e-7, 1e-9, 1e-9]  # first-pass bounds
xmax = [1e-3, 1e-4, 1e-4]

# rewrite the bounds in the way required by optimize.minimize
xbounds = [(low, high) for low, high in zip(xmin, xmax)]

blocksize = 0.5  # defines the fraction of the xrange that will be explored at each step
                 #  basinhopping starts with this value and reduces it by 10% every 'interval' iterations

mytakestep = MyTakeStep(blocksize, xmin, xmax)

minimizer_kwargs = dict(method=null_minimizer)
"""
result = optimize.basinhopping(rinp_error, x0, niter= 720, niter_success=100, disp=True, interval=20,
                                                                minimizer_kwargs=minimizer_kwargs, take_step=mytakestep)
"""
result = []
result.append(optimize.minimize(rinp_error, [x0[0]], args=['soma'], method='L-BFGS-B', bounds=[xbounds[0]],
                                    options={'ftol': 1e-4}))
cell.modify_mech_param('soma', 'pas', 'g', result[0].x[0])
cell.modify_mech_param('basal', 'pas', 'g', origin='soma')
result.append(optimize.minimize(rinp_error, [x0[1]], args=['trunk'], method='L-BFGS-B', bounds=[xbounds[1]],
                                    options={'ftol': 1e-4}))
cell.modify_mech_param('trunk', 'pas', 'g', origin='soma', slope=result[1].x[0])
cell.modify_mech_param('apical', 'pas', 'g', origin='trunk')
result.append(optimize.minimize(rinp_error, [x0[2]], args=['tuft'], method='L-BFGS-B', bounds=[xbounds[2]],
                                    options={'ftol': 1e-4}))
cell.modify_mech_param('trunk', 'pas', 'g', origin='trunk', slope=result[2].x[0])
new_x = [result[0].x[0], result[1].x[0], result[2].x[0]]
print('soma.g_pas: %.2E, trunk slope: %.2E, tuft slope: %.2E' % (new_x[0], new_x[1], new_x[2]))
plot_rinp_sec_types(new_x)