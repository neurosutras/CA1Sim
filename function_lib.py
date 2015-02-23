__author__ = 'milsteina'
import math
import pickle
import os.path
import datetime
import time
import numpy as np
import matplotlib.pyplot as plt
import h5py
from neuron import h


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


default_mech_dict = {'ais': {'cable': {'Ra': {'origin': 'soma'}},
                             'kdr': {'gkdrbar': {'origin': 'soma'}},
                             'hh2': {'gnabar': {'origin': 'self', 'value': 0.25},
                                     'gkbar': {'origin': 'soma'}},
                             'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'apical': {'cable': {'Ra': {'origin': 'soma'}},
                                'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'axon': {'cable': {'Ra': {'origin': 'soma'}},
                              'kdr': {'gkdrbar': {'origin': 'soma'}},
                              'hh2': {'gnabar': {'origin': 'soma'},
                                      'gkbar': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'basal': {'cable': {'Ra': {'origin': 'soma'}},
                               'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'soma': {'cable': {'Ra': {'origin': 'self', 'value': 150}},
                              'kdr': None,
                              'hh2': {'gnabar': {'origin': 'self', 'value': 0.05}},
                              'pas': {'e': {'origin': 'self', 'value': -65},
                                      'g': {'origin': 'self', 'value': 2.5e-05}}},
                     'trunk': {'cable': {'Ra': {'origin': 'soma'}},
                               'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}},
                     'tuft': {'cable': {'Ra': {'origin': 'soma'}},
                              'pas': {'e': {'origin': 'soma'}, 'g': {'origin': 'soma'}}}}
"""
default_mech_dict = {}


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


def scaleSWC(filenameBase, mag=100):
    # this function rescales the SWC file with the real distances.
    f = open(filenameBase+'.swc')
    lines = f.readlines()
    f.close()
    Points = []
    if mag == 100:
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
    f = open(filenameBase+'-scaled.swc','w')
    for [nn,tp,py,px,z,r,np] in Points:
        ll = str(int(nn))+' '+str(int(tp))+' '+str(py)+' '+str(px)+' '+str(z)+' '+str(r)+' '+str(int(np))+'\n'
        f.write(ll)
    f.close()


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
    def __init__(self, tstop=400, cvode=1, dt=None):
        self.rec_list = []  # list of dicts with keys for 'cell', 'node', 'loc' and 'vec': pointer to hoc Vector object.
                            # Also contains keys for 'ylabel' and 'units' for recording parameters other than Vm.
        self.stim_list = []  # list of dicts with keys for 'cell', 'node', 'stim': pointer to hoc IClamp object, and
                             # 'vec': recording of actual stimulus for plotting later
        self.tstop = tstop
        h.load_file('stdrun.hoc')
        h.celsius = 35.0
        if cvode:
            h.cvode_active(1)
            h.cvode.atol(0.0001)
        if dt is None:
            self.dt = h.dt
        else:
            self.dt = dt
        self.tvec = h.Vector()
        self.tvec.record(h._ref_t)
        self.parameters = {}

    def run(self, v_init=-65.):
        start_time = time.time()
        h.tstop = self.tstop
        h.dt = self.dt
        h.v_init = v_init
        h.init()
        h.run()
        print 'Simulation runtime: ', time.time()-start_time, ' sec'

    def append_rec(self, cell, node, loc=0.5, param='_ref_v', object=None, ylabel='Vm', units='mV', description=None):
        rec_dict = {'cell': cell, 'node': node, 'ylabel': ylabel, 'units': units}
        if not description is None:
            rec_dict['description'] = description
        rec_dict['vec'] = h.Vector()
        if object is None:
            rec_dict['vec'].record(getattr(node.sec(loc), param))
        else:
            try:
                loc = object.get_loc()
            except:
                loc = 0.5  # if the object doesn't have a .get_loc() method, default to 0.5
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
        if (node is None) | (loc is None):
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

    def modify_rec(self, index=0, node=None, loc=None, object=None, param=None, ylabel=None, units=None,
                                                                                        description=None):
        rec_dict = self.rec_list[index]
        if param is None:
            param = '_ref_v'
        else:
            if not ylabel is None:
                rec_dict['ylabel'] = ylabel
            if not units is None:
                rec_dict['units'] = units
        if object is None:
            if not node is None:
                rec_dict['node'] = node
                rec_dict['vec'].record(getattr(node.sec(rec_dict['loc']), param))
            if not loc is None:
                rec_dict['loc'] = loc
                rec_dict['vec'].record(getattr(rec_dict['node'].sec(loc), param))
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
        print 'Simulation ',simiter,': exporting took: ', time.time()-start_time, ' s'


def time2index(tvec, start, stop):
    """
    When using adaptive time step (cvode), indices corresponding to specific time points cannot be calculated from a
    fixed dt. This method returns the indices closest to the duration bounded by the specified time points.
    :param tvec: :class:'numpy.array'
    :param start: float
    :param stop: float
    :return: tuple of int
    """
    left = np.where(tvec >= start)[0][0]
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
    :return:
    """
    left, right = time2index(tvec, start-3., start-1.)
    baseline = np.average(vec[left:right])
    temp_vec = np.abs(vec - baseline)
    peak = np.max(temp_vec)
    left, right = time2index(tvec, stop-3., stop-1.)
    plateau = np.average(temp_vec[left:right])
    return peak/abs(amp), plateau/abs(amp)
