__author__ = 'milsteina'
import math
import pickle
import os.path
import datetime
import copy
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mm
import h5py
import scipy.optimize as optimize
import scipy.signal as signal
from neuron import h  # must be found in system $PYTHONPATH


#---------------------------------------Some global variables and functions------------------------------

data_dir = 'data/'
morph_dir = 'morphologies/'

freq = 100      # Hz, frequency at which AC length constant will be computed
d_lambda = 0.1  # no segment will be longer than this fraction of the AC length constant

"""
Structure of Mechanism Dictionary: dict of dicts

keys:               description:
'mechanism name':   Value is dictionary specifying how to set parameters at the mechanism level.
'cable':            Value is dictionary specifying how to set basic cable parameters at the section level. Includes
                        'Ra', 'cm', and the special parameter 'spatial_res', which scales the number of segments per
                        section for the specified sec_type by a factor of an exponent of 3.
'ions':             Value is dictionary specifying how to set parameters for ions at the section or segment level.
                    These parameters must be specified **after** all other mechanisms have been inserted.
values:
None:               Use default values for all parameters within this mechanism.
dict:
    keys:
    'parameter name':
    values:     dict:
                        keys:        value:
                        'origin':   'self':     Use 'value' as a baseline value.
                                    sec_type:   Inherit value from last seg of the closest node with sec of
                                                sec_type along the path to root.
                        'value':    float:      If 'origin' is 'self', contains the baseline value.
                        'slope':    float:      If exists, contains slope in units per um. If not, use
                                                constant 'value' for the full length of sec.
                        'max':      float:      If 'slope' exists, 'max' is an upper limit for the value
                        'min':      float:      If 'slope' exists, min is a lower limit for the value

"""

default_mech_dict = {'ais': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                             'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'apical': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                                'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'axon': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'basal': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                               'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'soma': {'cable': {'Ra': {'value': 150.}, 'cm': {'value': 1.}},
                              'pas': {'e': {'value': -67.}, 'g': {'value': 2.5e-05}}},
                     'trunk': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                               'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'tuft': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'spine_neck': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'spine_head': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}}}


def lambda_f(sec, f=freq):
    """
    Calculates the AC length constant for the given section at the frequency f
    Used to determine the number of segments per hoc section to achieve the desired spatial and temporal resolution
    :param sec : :class:'h.Section'
    :param f : int
    :return : int
    """
    diam = sec(0.5).diam
    Ra = sec.Ra
    cm = sec.cm
    return 1e5*math.sqrt(diam/(4*math.pi*f*Ra*cm))


def d_lambda_nseg(sec, lam=d_lambda, f=freq):
    """
    The AC length constant for this section and the user-defined fraction is used to determine the maximum size of each
    segment to achieve the desired spatial and temporal resolution. This method returns the number of segments to set
    the nseg parameter for this section. For tapered cylindrical sections, the diam parameter will need to be
    reinitialized after nseg changes.
    :param sec : :class:'h.Section'
    :param lam : int
    :param f : int
    :return : int
    """
    L = sec.L
    return int((L/(lam*lambda_f(sec, f))+0.9)/2)*2+1


def scaleSWC(filenameBase, mag=100, scope='neurolucida'):
    # this function rescales the SWC file with the real distances.
    f = open(morph_dir+filenameBase+'.swc')
    lines = f.readlines()
    f.close()
    Points = []
    if mag == 100:
        if scope == 'neurolucida':
            xyDist = 0.036909375  # 0.07381875
            zDist = 1.0
        else:
            xyDist = 0.065
            zDist = 0.05
    else:
        raise Exception('Calibration for {}X objective unknown.'.format(mag))
    for line in lines:
        ll = line.split(' ')
        nn = int(float(ll[0]))    # label of the point
        tp = int(float(ll[1]))  # point type
        py = float(ll[2])    # note the inversion of x, y.
        px = float(ll[3])
        z = float(ll[4])    # z
        r = float(ll[5])    # radius of the sphere.
        np = int(float(ll[6]))    # parent point id.
        # get the length in micron
        py *= xyDist; px *= xyDist; r = r*xyDist; z *= zDist
        Points.append([nn,tp,py,px,z,r,np])

    print 'Saving SWC to file '+filenameBase+'-scaled.swc'
    f = open(morph_dir+filenameBase+'-scaled.swc', 'w')
    for [nn,tp,py,px,z,r,np] in Points:
        ll = str(int(nn))+' '+str(int(tp))+' '+str(py)+' '+str(px)+' '+str(z)+' '+str(r)+' '+str(int(np))+'\n'
        f.write(ll)
    f.close()


def investigateSWC(filenameBase):
    # this function reports the min and max values for y, x, z, and radius from an SWC file.
    f = open(morph_dir+filenameBase+'.swc')
    lines = f.readlines()
    f.close()
    xvals = []
    yvals = []
    zvals = []
    rvals = []
    for line in lines:
        ll = line.split(' ')
        yvals.append(float(ll[2]))    # note the inversion of x, y.
        xvals.append(float(ll[3]))
        zvals.append(float(ll[4]))    # z
        rvals.append(float(ll[5]))    # radius of the sphere.
    print 'x - ',min(xvals),':',max(xvals)
    print 'y - ',min(yvals),':',max(yvals)
    print 'z - ',min(zvals),':',max(zvals)
    print 'r - ',min(rvals),':',max(rvals)


def translateSWCs():
    """
    Eric Bloss has produced high resolution .swc files that each contain a volume 10 um deep in z. This method
    determines from the filename the z offset of each file and translates the z coordinates of the .swc files to
    facilitate stitching them together into a single volume. Also changes the sec_type of any node that is not a root
    and has no children within a file to 7 to indicate a leaf that potentially needs to be connected to a nearby root.
    Also attempts to connect unconnected nodes that are within 2 um of each other across consecutive slices, and labels
    them with sec_type = 8. This doesn't work particularly well and files must be extensively proofread in NeuTuMac.
    """
    num_nodes = 0
    outputname = 'combined-offset-connected.swc'
    out_f = open(outputname, 'w')
    # out_test = open('combined-offset-connected.swc', 'w')
    prev_nodes = {}
    filenames = []
    z_offsets = []
    for filename in os.listdir('.'):
        if '.swc' in filename and not '-offset' in filename:
            filenames.append(filename)
            z_offsets.append(float(filename.split('z=')[1].split(' ')[0])/10.0)
    indexes = range(len(z_offsets))
    indexes.sort(key=z_offsets.__getitem__)
    for i in indexes:
        f = open(filenames[i])
        lines = f.readlines()
        f.close()
        num_nodes += len(prev_nodes)
        nodes = {}
        leaves = []
        for line in [line.split(' ') for line in lines if not line.split(' ')[0] in ['#', '\r\n']]:
            index = int(float(line[0])) + num_nodes  # node index
            nodes[index] = {}
            nodes[index]['type'] = int(float(line[1]))  # sec_type
            nodes[index]['y'] = float(line[2])  # note the inversion of x, y.
            nodes[index]['x'] = float(line[3])
            nodes[index]['z'] = float(line[4]) + z_offsets[i]
            nodes[index]['r'] = float(line[5])  # radius of the sphere.
            nodes[index]['parent'] = int(float(line[6]))  # index of parent node
            if not nodes[index]['parent'] == -1:
                nodes[index]['parent'] += num_nodes
                leaves.append(index)
        for index in nodes:  # keep nodes with no children
            parent = nodes[index]['parent']
            if parent in leaves:
                leaves.remove(parent)
        for index in leaves:
            nodes[index]['type'] = 7
        print 'Saving '+filenames[i]+' to '+outputname
        if prev_nodes:
            leaves = [index for index in nodes if (nodes[index]['type'] == 7 or nodes[index]['parent'] == -1)]
            for prev_index in [index for index in prev_nodes if (prev_nodes[index]['type'] == 7 or
                                                                prev_nodes[index]['parent'] == -1)]:
                for index in leaves:
                    distance = math.sqrt((prev_nodes[prev_index]['x']-nodes[index]['x'])**2 +
                                         (prev_nodes[prev_index]['y']-nodes[index]['y'])**2 +
                                         (prev_nodes[prev_index]['z']-nodes[index]['z'])**2)
                    # print prev_index, index, distance
                    if distance < 2.:
                        prev_nodes[prev_index]['type'] = 8
                        nodes[index]['type'] = 8
                        nodes[index]['parent'] = prev_index
                        leaves.remove(index)
                        break
        for index in prev_nodes:
            line = str(index)+' '+str(prev_nodes[index]['type'])+' '+str(prev_nodes[index]['y'])+' '+\
                   str(prev_nodes[index]['x'])+' '+str(prev_nodes[index]['z'])+' '+str(prev_nodes[index]['r'])+' '+\
                   str(prev_nodes[index]['parent'])+'\n'
            out_f.write(line)
        prev_nodes = copy.deepcopy(nodes)
    for index in prev_nodes:
        line = str(index)+' '+str(prev_nodes[index]['type'])+' '+str(prev_nodes[index]['y'])+' '+\
               str(prev_nodes[index]['x'])+' '+str(prev_nodes[index]['z'])+' '+str(prev_nodes[index]['r'])+' '+\
               str(prev_nodes[index]['parent'])+'\n'
        out_f.write(line)
    out_f.close()


def write_to_pkl(fname, data):
    """
    HocCell objects maintain a nested dictionary specifying membrane mechanism parameters for each subcellular
    compartment. This method is used to save that dictionary to a .pkl file that can be read in during model
    specification or after parameter optimization.
    :param fname: str
    :param data: picklable object
    """
    output = open(fname, 'wb')
    pickle.dump(data, output, 2)
    output.close()


def read_from_pkl(fname):
    """
    HocCell objects maintain a nested dictionary specifying membrane mechanism parameters for each subcellular
    compartment. This method is used to load that dictionary from a .pkl file during model specification.
    :param fname: str
    :return: unpickled object
    """
    if os.path.isfile(fname):
        pkl_file = open(fname, 'rb')
        data = pickle.load(pkl_file)
#        pprint.pprint(data)
        pkl_file.close()
        return data
    else:
        raise Exception('File: {} does not exist.'.format(fname))


class QuickSim(object):
    """
    This method is used to run a quick simulation with a set of current injections and a set of recording sites.
    Can save detailed information about the simulation to an HDF5 file after each run. Once defined, IClamp objects
    persist when using an interactive console, but not when executing standalone scripts. Therefore, the best practice
    is simply to set amp to zero to turn off current injections, or move individual IClamp processes to different
    locations rather then adding and deleting them.
    class params:
    self.stim_list:
    self.rec_list:
    """
    def __init__(self, tstop=400., cvode=1, dt=None, verbose=1):
        self.rec_list = []  # list of dicts with keys for 'cell', 'node', 'loc' and 'vec': pointer to hoc Vector object.
                            # Also contains keys for 'ylabel' and 'units' for recording parameters other than Vm.
        self.stim_list = []  # list of dicts with keys for 'cell', 'node', 'stim': pointer to hoc IClamp object, and
                             # 'vec': recording of actual stimulus for plotting later
        self.tstop = tstop
        h.load_file('stdrun.hoc')
        h.celsius = 35.0
        if cvode:
            self.cvode = h.CVode()
            self.cvode.active(1)
            self.cvode.atol(0.001)  # 0.0001
        else:
            self.cvode = None
        if dt is None:
            self.dt = h.dt
        else:
            self.dt = dt
        self.verbose = verbose
        self.tvec = h.Vector()
        self.tvec.record(h._ref_t)
        self.parameters = {}

    def run(self, v_init=-65.):
        start_time = time.time()
        h.tstop = self.tstop
        if self.cvode is None:
            h.dt = self.dt
        h.v_init = v_init
        h.init()
        h.run()
        if self.verbose:
            print 'Simulation runtime: ', time.time()-start_time, ' sec'

    def append_rec(self, cell, node, loc=None, param='_ref_v', object=None, ylabel='Vm', units='mV', description=None):
        rec_dict = {'cell': cell, 'node': node, 'ylabel': ylabel, 'units': units}
        if not description is None:
            rec_dict['description'] = description
        rec_dict['vec'] = h.Vector()
        if object is None:
            if loc is None:
                loc = 0.5
            rec_dict['vec'].record(getattr(node.sec(loc), param))
        else:
            if loc is None:
                try:
                    loc = object.get_segment().x  # this should not push the section to the hoc stack or risk overflow
                except:
                    loc = 0.5  # if the object doesn't have a .get_loc() method, default to 0.5
            if param is None:
                rec_dict['vec'].record(object)
            else:
                rec_dict['vec'].record(getattr(object, param))
        rec_dict['loc'] = loc
        self.rec_list.append(rec_dict)

    def append_stim(self, cell, node, loc, amp, delay, dur, description='IClamp'):
        stim_dict = {'cell': cell, 'node': node, 'description': description}
        stim_dict['stim'] = h.IClamp(node.sec(loc))
        stim_dict['stim'].amp = amp
        stim_dict['stim'].delay = delay
        stim_dict['stim'].dur = dur
        stim_dict['vec'] = h.Vector()
        stim_dict['vec'].record(stim_dict['stim']._ref_i)
        self.stim_list.append(stim_dict)

    def modify_stim(self, index=0, node=None, loc=None, amp=None, delay=None, dur=None, description=None):
        stim_dict = self.stim_list[index]
        if not (node is None and loc is None):
            if not node is None:
                stim_dict['node'] = node
            if loc is None:
                loc = stim_dict['stim'].get_segment().x
            stim_dict['stim'].loc(stim_dict['node'].sec(loc))
        if not amp is None:
            stim_dict['stim'].amp = amp
        if not delay is None:
            stim_dict['stim'].delay = delay
        if not dur is None:
            stim_dict['stim'].dur = dur
        if not description is None:
            stim_dict['description'] = description

    def modify_rec(self, index=0, node=None, loc=None, object=None, param='_ref_v', ylabel=None, units=None,
                                                                                        description=None):
        rec_dict = self.rec_list[index]
        if not ylabel is None:
            rec_dict['ylabel'] = ylabel
        if not units is None:
            rec_dict['units'] = units
        if not node is None:
            rec_dict['node'] = node
        if not loc is None:
            rec_dict['loc'] = loc
        if object is None:
            rec_dict['vec'].record(getattr(rec_dict['node'].sec(rec_dict['loc']), param))
        elif param is None:
            rec_dict['vec'].record(object)
        else:
            rec_dict['vec'].record(getattr(object, param))
        if not description is None:
            rec_dict['description'] = description

    def plot(self):
        for rec_dict in self.rec_list:
            if 'description' in rec_dict:
                description = str(rec_dict['description'])
            else:
                description = ''
            plt.plot(self.tvec, rec_dict['vec'], label=rec_dict['node'].name+'('+str(rec_dict['loc'])+') - '+
                                                                                description)
            plt.xlabel("Time (ms)")
            plt.ylabel(rec_dict['ylabel']+' ('+rec_dict['units']+')')
        plt.legend(loc='upper right')
        if 'description' in self.parameters:
            plt.title(self.parameters['description'])
        plt.show()
        plt.close()

    def export_to_file(self, f, simiter=0):
        """
        Extracts important parameters from the lists of stimulation and recording sites, and exports to an HDF5
        database. Arrays are saved as datasets and metadata is saved as attributes.
        :param f: :class:'h5py.File'
        :param simiter: int
        """
        start_time = time.time()
        f.create_group(str(simiter))
        f[str(simiter)].create_dataset('time', compression='gzip', compression_opts=9, data=self.tvec)
        f[str(simiter)]['time'].attrs['dt'] = self.dt
        for parameter in self.parameters:
            f[str(simiter)].attrs[parameter] = self.parameters[parameter]
        if self.stim_list:
            f[str(simiter)].create_group('stim')
            for index, stim in enumerate(self.stim_list):
                stim_out = f[str(simiter)]['stim'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                                  data=stim['vec'])
                cell = stim['cell']
                stim_out.attrs['cell'] = cell.gid
                node = stim['node']
                stim_out.attrs['index'] = node.index
                stim_out.attrs['type'] = node.type
                loc = stim['stim'].get_segment().x
                stim_out.attrs['loc'] = loc
                distance = cell.get_distance_to_node(cell.tree.root, node, loc)
                stim_out.attrs['soma_distance'] = distance
                distance = cell.get_distance_to_node(cell.get_dendrite_origin(node), node, loc)
                stim_out.attrs['branch_distance'] = distance
                stim_out.attrs['amp'] = stim['stim'].amp
                stim_out.attrs['delay'] = stim['stim'].delay
                stim_out.attrs['dur'] = stim['stim'].dur
                stim_out.attrs['description'] = stim['description']
        f[str(simiter)].create_group('rec')
        for index, rec in enumerate(self.rec_list):
            rec_out = f[str(simiter)]['rec'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                            data=rec['vec'])
            cell = rec['cell']
            rec_out.attrs['cell'] = cell.gid
            node = rec['node']
            rec_out.attrs['index'] = node.index
            rec_out.attrs['type'] = node.type
            rec_out.attrs['loc'] = rec['loc']
            distance = cell.get_distance_to_node(cell.tree.root, node, rec['loc'])
            rec_out.attrs['soma_distance'] = distance
            distance = cell.get_distance_to_node(cell.get_dendrite_origin(node), node, rec['loc'])
            rec_out.attrs['branch_distance'] = distance
            rec_out.attrs['ylabel'] = rec['ylabel']
            rec_out.attrs['units'] = rec['units']
            if 'description' in rec:
                rec_out.attrs['description'] = rec['description']
        if self.verbose:
            print 'Simulation ', simiter, ': exporting took: ', time.time()-start_time, ' s'


class Normalized_Step(object):
    """
    For use with scipy.optimize packages like basinhopping that allow a customized step-taking method.
    Converts basinhopping absolute stepsize into different stepsizes for each parameter such that the stepsizes are
    some fraction of the ranges specified by xmin and xmax. Also enforces bounds for x, and explores the range in
    log10 space when the range is greater than 2 orders of magnitude.
    xmin and xmax are delivered as raw, not relative values. Can handle negative values and ranges that cross zero. If
    xmin and xmax are not provided, or contain None as values, the default is 0.1 and 10. * x0.
    """
    def __init__(self, x0, xmin=None, xmax=None, stepsize=0.5):
        self.stepsize = stepsize
        if xmin is None:
            xmin = [None for i in range(len(x0))]
        for i in range(len(xmin)):
            if xmin[i] is None:
                if x0[i] > 0.:
                    xmin[i] = 0.1 * x0[i]
                else:
                    xmin[i] = 10. * x0[i]
        if xmax is None:
            xmax = [None for i in range(len(x0))]
        for i in range(len(xmax)):
            if xmax[i] is None:
                if x0[i] > 0.:
                    xmax[i] = 10. * x0[i]
                else:
                    xmax[i] = 0.1 * x0[i]
        self.x0 = x0
        self.x_range = np.subtract(xmax, xmin)
        self.order_mag = np.abs(np.log10(np.abs(np.divide(xmax, xmin))))
        self.log10_range = np.log10(np.add(1., self.x_range))
        self.x_offset = np.subtract(1., xmin)

    def __call__(self, current_x):
        x = np.add(current_x, self.x_offset)
        x = np.maximum(x, 1.)
        x = np.minimum(x, np.add(1., self.x_range))
        for i in range(len(x)):
            if self.order_mag[i] >= 2.:
                x[i] = self.log10_step(i, x[i])
            else:
                x[i] = self.linear_step(i, x[i])
        new_x = np.subtract(x, self.x_offset)
        return new_x

    def linear_step(self, i, xi):
        step = self.stepsize * self.x_range[i] / 2.
        new_xi = np.random.uniform(max(1., xi-step), min(xi+step, 1.+self.x_range[i]))
        return new_xi

    def log10_step(self, i, xi):
        step = self.stepsize * self.log10_range[i] / 2.
        xi = np.log10(xi)
        new_xi = np.random.uniform(max(0., xi-step), min(xi+step, self.log10_range[i]))
        new_xi = np.power(10., new_xi)
        return new_xi


def combine_output_files(rec_file_list, new_rec_filename=None):
    """
    List contains names of files generated by "embarassingly parallel" execution of simulations on separate cores.
    This function combines the contents of the files into one .hdf5 file.
    :param rec_file_list: list
    :param new_rec_filename: str or None
    """
    if new_rec_filename is None:
        new_rec_filename = 'combined_output_'+datetime.datetime.today().strftime('%m%d%Y%H%M')
    new_f = h5py.File(data_dir+new_rec_filename+'.hdf5', 'w')
    simiter = 0
    for rec_filename in rec_file_list:
        old_f = h5py.File(data_dir+rec_filename+'.hdf5', 'r')
        for old_group in old_f.itervalues():
            new_f.copy(old_group, new_f, name=str(simiter))
            simiter += 1
        old_f.close()
    new_f.close()
    print 'Combined data in list of files and exported to: '+new_rec_filename+'.hdf5'
    return new_rec_filename


def time2index(tvec, start, stop):
    """
    When using adaptive time step (cvode), indices corresponding to specific time points cannot be calculated from a
    fixed dt. This method returns the indices closest to the duration bounded by the specified time points.
    :param tvec: :class:'numpy.array'
    :param start: float
    :param stop: float
    :return: tuple of int
    """
    left = np.where(tvec >= start)[0]
    if np.any(left):  # at least one value was found
        left = left[0]
    else:
        right = len(tvec) - 1  # just take the last two indices
        left = right - 1
        return left, right
    if tvec[left] >= stop:
        right = left
        left -= 1
        return left, right
    right = np.where(tvec <= stop)[0][-1]
    if right == left:
        left -= 1
    return left, right


def get_Rinp(tvec, vec, start, stop, amp):
    """
    Calculate peak and steady-state input resistance from a step current injection. For waveform current injections, the
    peak but not the steady-state will have meaning.
    :param tvec:
    :param vec:
    :param start:
    :param stop:
    :param amp:
    :return: tuple of float
    """

    interp_t = np.arange(0, stop, 0.01)
    interp_vm = np.interp(interp_t, tvec, vec)
    left, right = time2index(interp_t, start-3., start-1.)
    baseline = np.average(interp_vm[left:right])
    temp_vec = np.abs(interp_vm - baseline)
    peak = np.max(temp_vec[right:])
    left, right = time2index(interp_t, stop-3., stop-1.)
    plateau = np.average(temp_vec[left:right])
    return baseline, peak/abs(amp), plateau/abs(amp)


def model_exp_rise_decay(t, tau_rise, tau_decay):
    shape = np.exp(-t/tau_decay)-np.exp(-t/tau_rise)
    return shape/np.max(shape)


def model_exp_rise(t, tau):
    return 1-np.exp(-t/tau)


def model_exp_decay(t, tau):
    return np.exp(-t/tau)


def model_scaled_exp(t, A, tau, A0=0):
    return A*np.exp(t/tau)+A0


def null_minimizer(fun, x0, args, **options):
    """
    Rather than allow basinhopping to pass each local mimimum to a gradient descent algorithm for polishing, this method
    just catches and passes all local minima so basinhopping can proceed.
    """
    return optimize.OptimizeResult(x=x0, fun=fun(x0, *args), success=True, nfev=1)


class MyTakeStep(object):
    """
    For use with scipy.optimize packages like basinhopping that allow a customized step-taking method.
    Converts basinhopping absolute stepsize into different stepsizes for each parameter such that the stepsizes are
    some fraction of the ranges specified by xmin and xmax. Also enforces bounds for x, and explores the range in
    log space when the range is greater than 3 orders of magnitude.
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
            #  chooses new value in log space to allow fair sampling across orders of magnitude
            if np.log10(self.xmax[i]) - np.log10(self.xmin[i]) >= 3.:
                x[i] = np.power(10, np.random.uniform(np.log10(x[i]-sdec), np.log10(x[i]+sinc)))
            else:
                x[i] += np.random.uniform(-sdec, sinc)
        return x


def get_expected_spine_index_map(sim_file):
    """
    There is a bug with HDF5 when reading from a file too often within a session. Instead of constantly reading from the
    HDF5 file directly and searching for content by spine_index or path_index, the number of calls to the sim_file can
    be reduced by creating a mapping from spine_index or path_index to HDF5 group key. It is possible for a spine to
    have more than one entry in an expected_file, with branch recordings in different locations and therefore different
    expected EPSP waveforms, so it is necessary to also distinguish those entries by path_index.

    :param sim_file:
    :return: dict
    """
    index_map = {}
    for key, sim in sim_file.iteritems():
        path_index = sim.attrs['path_index']
        spine_index = sim.attrs['spine_index']
        if path_index not in index_map:
            index_map[path_index] = {}
        index_map[path_index][spine_index] = key
    return index_map


def get_expected_EPSP(sim_file, group_index, equilibrate, duration):
    """
    Given an output file generated by parallel_expected_EPSP, this method returns a dict of numpy arrays, each
    containing the depolarization-rectified expected EPSP for each recording site from stimulating a single spine.
    :param sim_file: :class:'h5py.File'
    :param group_index: int
    :param equilibrate: float
    :param duration: float
    :return: dict of :class:'numpy.array'
    """
    sim = sim_file[str(group_index)]
    t = sim['time'][:]
    interp_t = np.arange(0, duration, 0.01)
    left, right = time2index(interp_t, equilibrate-3., equilibrate-1.)
    start, stop = time2index(interp_t, equilibrate-2., duration)
    trace_dict = {}
    for rec in sim['rec'].itervalues():
        location = rec.attrs['description']
        vm = rec[:]
        interp_vm = np.interp(interp_t, t, vm)
        baseline = np.average(interp_vm[left:right])
        interp_vm -= baseline
        interp_vm = interp_vm[start:stop]
        rectified = np.zeros(len(interp_vm))
        rectified[np.where(interp_vm>0)[0]] += interp_vm[np.where(interp_vm>0)[0]]
        trace_dict[location] = rectified
    interp_t = interp_t[start:stop]
    interp_t -= interp_t[0] + 2.
    trace_dict['time'] = interp_t
    return trace_dict


def get_expected_vs_actual(expected_sim_file, actual_sim_file, expected_index_map, sorted_actual_sim_keys,
                           interval=0.3):
    """
    Given an output file generated by parallel_expected_EPSP, and an output file generated by
    parallel_branch_cooperativity, this method returns a dict of lists, each containing an input-output function
    relating expected to actual peak depolarization for each recording site from stimulating a group of spines on a
    single branch or path. The variable expected_index_map contains a dictionary that converts an integer spine_index to
    a string group_index to locate the expected EPSP for a given spine in the expected_sim_file. The variable
    sorted_actual_sim_keys contains the indexes of the simulations in the actual_sim_file corresponding to the branch or
    path, ordered by number of stimulated spines. These variables must be pre-computed.
    :param expected_sim_file: :class:'h5py.File'
    :param actual_sim_file: :class:'h5py.File'
    :param expected_index_map: dict
    :param sorted_actual_sim_keys: list of str
    :param interval: float
    :return: dict of list
    """
    equilibrate = actual_sim_file[sorted_actual_sim_keys[0]].attrs['equilibrate']
    duration = actual_sim_file[sorted_actual_sim_keys[0]].attrs['duration']
    actual = {}
    for sim in [actual_sim_file[key] for key in sorted_actual_sim_keys]:
        t = sim['time'][:]
        interp_t = np.arange(0, duration, 0.01)
        left, right = time2index(interp_t, equilibrate-3., equilibrate-1.)
        start, stop = time2index(interp_t, equilibrate-2., duration)
        for rec in sim['rec'].itervalues():
            location = rec.attrs['description']
            if not location in actual:
                actual[location] = []
            vm = rec[:]
            interp_vm = np.interp(interp_t, t, vm)
            baseline = np.average(interp_vm[left:right])
            interp_vm -= baseline
            interp_vm = interp_vm[start:stop]
            actual[location].append(np.max(interp_vm))
    spine_list = sim.attrs['syn_indexes']
    interp_t = interp_t[start:stop]
    interp_t -= interp_t[0] + 2.
    expected = {}
    summed_traces = {}
    equilibrate = expected_sim_file.values()[0].attrs['equilibrate']
    duration = expected_sim_file.values()[0].attrs['duration']
    for i, spine_index in enumerate(spine_list):
        group_index = expected_index_map[spine_index]
        trace_dict = get_expected_EPSP(expected_sim_file, group_index, equilibrate, duration)
        t = trace_dict['time']
        left, right = time2index(interp_t, -2.+i*interval, interp_t[-1])
        right = min(right, left+len(t))
        for location in [location for location in trace_dict if not location == 'time']:
            trace = trace_dict[location]
            if not location in expected:
                expected[location] = []
                summed_traces[location] = np.zeros(len(interp_t))
            summed_traces[location][left:right] += trace[:right-left]
            expected[location].append(np.max(summed_traces[location]))
    return expected, actual


def export_nmdar_cooperativity(expected_filename, actual_filename, description="", output_filename=None):
    """
    Expects expected and actual files to be generated by parallel_clustered_ or
    parallel_distributed_branch_cooperativity. Files contain simultaneous voltage recordings from 3 locations (soma,
    trunk, dendrite origin) during synchronous stimulation of branches, and an NMDAR conductance recording from a single
    spine in each group. Spines are distributed across 4 dendritic sec_types (basal, trunk, apical, tuft).
    Generates a processed output file containing expected vs. actual data and metadata for each group of spines.
    Can be used to generate plots of supralinearity, NMDAR conductance, or average across conditions, etc.
    :param expected_filename: str
    :param actual_filename: str
    :param description: str
    :param output_filename: str
    """
    sim_key_dict = {}
    with h5py.File(data_dir+actual_filename+'.hdf5', 'r') as actual_file:
        for key, sim in actual_file.iteritems():
            path_index = sim.attrs['path_index']
            if path_index not in sim_key_dict:
                sim_key_dict[path_index] = []
            sim_key_dict[path_index].append(key)
        with h5py.File(data_dir+expected_filename+'.hdf5', 'r') as expected_file:
            expected_index_map = get_expected_spine_index_map(expected_file)
            with h5py.File(data_dir+output_filename+'.hdf5', 'w') as output_file:
                output_file.attrs['description'] = description
                for path_index in sim_key_dict:
                    path_group = output_file.create_group(str(path_index))
                    sim_keys = sim_key_dict[path_index]
                    sim_keys.sort(key=lambda x: len(actual_file[x].attrs['syn_indexes']))
                    sim = actual_file[sim_keys[0]]
                    path_type = sim.attrs['path_type']
                    path_category = sim.attrs['path_category']
                    soma_distance = sim['rec']['4'].attrs['soma_distance']
                    branch_distance = sim['rec']['4'].attrs['branch_distance']
                    origin_distance = soma_distance - branch_distance
                    path_group.attrs['path_type'] = path_type
                    path_group.attrs['path_category'] = path_category
                    path_group.attrs['soma_distance'] = soma_distance
                    path_group.attrs['branch_distance'] = branch_distance
                    path_group.attrs['origin_distance'] = origin_distance
                    expected_dict, actual_dict = get_expected_vs_actual(expected_file, actual_file,
                                                                        expected_index_map[path_index], sim_keys)
                    for rec in sim['rec'].itervalues():
                        location = rec.attrs['description']
                        rec_group = path_group.create_group(location)
                        rec_group.create_dataset('expected', compression='gzip', compression_opts=9,
                                                           data=expected_dict[location])
                        rec_group.create_dataset('actual', compression='gzip', compression_opts=9,
                                                           data=actual_dict[location])


def get_binned_firing_rate(spike_times, t, bin_dur=10.):
    """

    :param spike_times: list of lists
    :param t: array
    :param bin_dur: float (ms)
    :return: array (Hz)
    """
    bin_centers = np.arange(t[0]+bin_dur/2., t[-1], bin_dur)
    count = np.zeros(len(bin_centers))
    for train in spike_times:
        for spike_time in train:
            if t[0] <= spike_time <= bin_centers[-1] + bin_dur / 2.:
                i = np.where(bin_centers + bin_dur / 2. >= spike_time)[0][0]
                count[i] += 1
    rate = count / bin_dur * 1000.
    rate = np.interp(t, bin_centers, rate)
    return rate


def get_smoothed_firing_rate(spike_times, t, bin_dur=50., bin_step=10., dt=0.1):
    """

    :param spike_times: list of lists
    :param t: array
    :param bin_dur: float (ms)
    :param dt: float (ms)
    :return: array (Hz)
    """
    count = np.zeros(len(t))
    bin_centers = np.arange(t[0]+bin_dur/2., t[-1]-bin_dur/2., bin_step)
    rate = np.zeros(len(bin_centers))
    for train in spike_times:
        for spike_time in train:
            if spike_time >= t[0] and spike_time <= t[-1]:
                i = np.where(t >= spike_time)[0][0]
                count[i] += 1
    left = 0
    right = int(bin_dur / dt)
    interval = int(bin_step / dt)
    for i, bin_t in enumerate(bin_centers):
        rate[i] = np.sum(count[left:right]) / bin_dur * 1000.
        left += interval
        right += interval
    rate = np.interp(t, bin_centers, rate)
    return rate


def get_removed_spikes(rec_filename, before=1.6, after=6., dt=0.02, th=20., plot=1):
    """

    :param rec_filename: str
    :param before: float : time to remove before spike
    :param after: float : time to remove after spike in case where trough or voltage recovery cannot be used
    :param dt: float : temporal resolution for interpolation and dvdt
    :param th: float : slope threshold
    :param plot: int
    :return: list of :class:'numpy.array'
    """
    removed = []
    with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
        sim = f.values()[0]
        equilibrate = sim.attrs['equilibrate']
        duration = sim.attrs['duration']
        track_equilibrate = sim.attrs['track_equilibrate']
        for rec in f.values():
            t = np.arange(0., duration, dt)
            vm = np.interp(t, rec['time'], rec['rec']['0'])
            start = int((equilibrate + track_equilibrate) / dt)
            t = np.subtract(t[start:], equilibrate + track_equilibrate)
            vm = vm[start:]
            dvdt = np.gradient(vm, [dt])
            crossings = np.where(dvdt >= th)[0]
            if not np.any(crossings):
                removed.append(vm)
            else:
                start = 0.
                i = 0
                left = max(0, crossings[i] - int(before/dt))
                previous_vm = vm[left]
                while i < len(crossings) and start < len(vm):
                    start = crossings[i]
                    left = max(0, crossings[i] - int(before/dt))
                    if not np.isnan(vm[left]):
                        previous_vm = vm[left]
                    recovers = np.where(vm[crossings[i]:] < previous_vm)[0]
                    if np.any(recovers):
                        recovers = crossings[i] + recovers[0]
                    falling = np.where(dvdt[crossings[i]:] < 0.)[0]
                    if np.any(falling):
                        falling = crossings[i] + falling[0]
                        rising = np.where(dvdt[falling:] >= 0.)[0]
                        if np.any(rising):
                            rising = falling + rising[0]
                    else:
                        rising = []
                    if np.any(recovers):
                        if np.any(rising):
                            right = min(recovers, rising)
                        else:
                            right = recovers
                    elif np.any(rising):
                        right = rising
                    else:
                        right = min(crossings[i] + int(after/dt), len(vm)-1)
                    # added to remove majority of complex spike:
                    if vm[right] >= -45. and np.any(recovers):
                        right = recovers
                    for j in range(left, right):
                        vm[j] = np.nan
                    i += 1
                    while i < len(crossings) and crossings[i] < right:
                        i += 1
                not_blank = np.where(~np.isnan(vm))[0]
                vm = np.interp(t, t[not_blank], vm[not_blank])
                removed.append(vm)
            temp_t = np.arange(0., duration, dt)
            temp_vm = np.interp(temp_t, rec['time'], rec['rec']['0'])
            start = int((equilibrate + track_equilibrate) / dt)
            temp_t = np.subtract(temp_t[start:], equilibrate + track_equilibrate)
            temp_vm = temp_vm[start:]
            if plot:
                plt.plot(temp_t, temp_vm)
                plt.plot(t, vm)
                plt.show()
                plt.close()
    return removed


def get_theta_filtered_traces(rec_filename, dt=0.02):
    """

    :param rec_file_name: str
    # remember .attrs['phase_offset'] could be inside ['train'] for old files
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
        sim = f.values()[0]
        equilibrate = sim.attrs['equilibrate']
        track_equilibrate = sim.attrs['track_equilibrate']
        track_length = sim.attrs['track_length']
        input_field_duration = sim.attrs['input_field_duration']
        duration = sim.attrs['duration']
        stim_dt = sim.attrs['stim_dt']
        track_duration = duration - equilibrate - track_equilibrate
        stim_t = np.arange(-track_equilibrate, track_duration, stim_dt)
        exc_input = [get_binned_firing_rate(sim['train'].values(), stim_t) for sim in f.values()]
        inh_input = [get_binned_firing_rate(sim['inh_train'].values(), stim_t) for sim in f.values()]
        phase_offsets = []
        for sim in f.values():
            if 'phase_offset' in sim.attrs:
                phase_offsets.append(sim.attrs['phase_offset'])
            elif 'phase_offset' in sim['train'].attrs:
                phase_offsets.append(sim['train'].attrs['phase_offset'])
    rec_t = np.arange(0., track_duration, dt)
    spikes_removed = get_removed_spikes(rec_filename, plot=0)
    # down_sample traces to 2 kHz after clipping spikes for theta and ramp filtering
    down_dt = 0.5
    down_stim_t = np.arange(-track_equilibrate, track_duration, down_dt)
    down_rec_t = np.arange(0., track_duration, down_dt)
    # 2000 ms Hamming window, ~3 Hz low-pass for ramp, ~5 - 10 Hz bandpass for theta
    window_len = min(int(2000./down_dt), len(down_rec_t) - 1)
    theta_filter = signal.firwin(window_len, [5., 10.], nyq=1000./2./down_dt, pass_zero=False)
    pop_exc_theta = []
    pop_inh_theta = []
    intra_theta = []
    for pop in exc_input:
        down_sampled = np.interp(down_stim_t, stim_t, pop)
        filtered = signal.filtfilt(theta_filter, [1.], down_sampled, padtype='even', padlen=window_len)
        up_sampled = np.interp(stim_t, down_stim_t, filtered)
        pop_exc_theta.append(up_sampled)
    for pop in inh_input:
        down_sampled = np.interp(down_stim_t, stim_t, pop)
        filtered = signal.filtfilt(theta_filter, [1.], down_sampled, padtype='even', padlen=window_len)
        up_sampled = np.interp(stim_t, down_stim_t, filtered)
        pop_inh_theta.append(up_sampled)
    for trace in spikes_removed:
        down_sampled = np.interp(down_rec_t, rec_t, trace)
        filtered = signal.filtfilt(theta_filter, [1.], down_sampled, padtype='even', padlen=window_len)
        up_sampled = np.interp(rec_t, down_rec_t, filtered)
        intra_theta.append(up_sampled)
    return stim_t, pop_exc_theta, pop_inh_theta, rec_t, intra_theta, phase_offsets


def get_phase_precession(rec_filename, start_loc=None, end_loc=None, dt=0.02):
    """

    :param rec_file_name: str
    # remember .attrs['phase_offset'] could be inside ['train'] for old files
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
        sim = f.values()[0]
        equilibrate = sim.attrs['equilibrate']
        track_equilibrate = sim.attrs['track_equilibrate']
        track_length = sim.attrs['track_length']
        input_field_duration = sim.attrs['input_field_duration']
        duration = sim.attrs['duration']
        stim_dt = sim.attrs['stim_dt']
        track_duration = duration - equilibrate - track_equilibrate
        stim_t = np.arange(-track_equilibrate, track_duration, stim_dt)
        if 'global_theta_cycle_duration' in sim.attrs:
            theta_duration = sim.attrs['global_theta_cycle_duration']
        else:
            theta_duration = 150.
        if start_loc is None:
            start_loc = 0.
        if end_loc is None:
            end_loc = track_duration
        phase_offsets = []
        for sim in f.values():
            if 'phase_offset' in sim.attrs:
                phase_offsets.append(sim.attrs['phase_offset'])
            elif 'train' in sim and 'phase_offset' in sim['train'].attrs:
                phase_offsets.append(sim['train'].attrs['phase_offset'])
            else:
                phase_offsets.append(0.)
        output_trains = [np.array(sim['output']) for sim in f.values() if 'output' in sim]
    spike_phase_array = []
    spike_time_array = []
    for i, train in enumerate(output_trains):
        time_offset = phase_offsets[i]
        on_track = np.where((train >= start_loc) & (train <= end_loc))[0]
        if not np.any(on_track):
            spike_phase_array.append([])
            spike_time_array.append([])
        else:
            spike_times = train[on_track]
            spike_time_array.append(spike_times)
            spike_times = np.subtract(spike_times, time_offset)
            spike_phases = np.mod(spike_times, theta_duration)
            spike_phases /= theta_duration
            spike_phases *= 360.
            spike_phase_array.append(spike_phases)
    rec_t = np.arange(0., track_duration, dt)
    spikes_removed = get_removed_spikes(rec_filename, plot=0)
    # down_sample traces to 2 kHz after clipping spikes for theta and ramp filtering
    down_dt = 0.5
    down_rec_t = np.arange(0., track_duration, down_dt)
    # 2000 ms Hamming window, ~3 Hz low-pass for ramp, ~5 - 10 Hz bandpass for theta
    window_len = min(len(down_rec_t) - 1, int(2000./down_dt))
    theta_filter = signal.firwin(window_len, [5., 10.], nyq=1000./2./down_dt, pass_zero=False)
    intra_theta = []
    for trace in spikes_removed:
        down_sampled = np.interp(down_rec_t, rec_t, trace)
        filtered = signal.filtfilt(theta_filter, [1.], down_sampled, padtype='even', padlen=window_len)
        up_sampled = np.interp(rec_t, down_rec_t, filtered)
        intra_theta.append(up_sampled)
    intra_peak_array = []
    intra_phase_array = []
    for i, trial in enumerate(intra_theta):
        time_offset = phase_offsets[i]
        peak_locs = signal.argrelmax(trial)[0]
        peak_times = rec_t[peak_locs]
        intra_peak_array.append(peak_times)
        peak_times = np.subtract(peak_times, time_offset)
        peak_phases = np.mod(peak_times, theta_duration)
        peak_phases /= theta_duration
        peak_phases *= 360.
        intra_phase_array.append(peak_phases)
    return spike_time_array, spike_phase_array, intra_peak_array, intra_phase_array


def get_subset_downsampled_recordings(rec_filename, description, dt=0.1):
    """

    :param rec_file_name: str
    # remember .attrs['phase_offset'] could be inside ['train'] for old files
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
        sim = f.values()[0]
        equilibrate = sim.attrs['equilibrate']
        track_equilibrate = sim.attrs['track_equilibrate']
        duration = sim.attrs['duration']
        rec_t = np.arange(0., duration, dt)
        sim_list = []
        for sim in f.values():
            rec_list = []
            index_list = []
            for rec in [rec for rec in sim['rec'].values() if 'description' in rec.attrs and rec.attrs['description']
                                                                                                        == description]:
                down_sampled = np.interp(rec_t, sim['time'], rec)
                rec_list.append(down_sampled[int((equilibrate + track_equilibrate) / dt):])
                index_list.append(rec.attrs['index'])
            sim_list.append({'index_list': index_list, 'rec_list': rec_list})
    rec_t = np.arange(0., duration - track_equilibrate - equilibrate, dt)
    return rec_t, sim_list