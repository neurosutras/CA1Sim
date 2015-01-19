__author__ = 'Aaron D. Milstein'
import math
import pickle
import os.path
import datetime
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
"""

default_mech_dict = {'soma': {'pas': {'g': {'origin': 'self', 'value': 0.0001}},
                              'hh2': {'gnabar': {'origin': 'self', 'value': 0.05}}},
                     'ais': {'hh2': {'gnabar': {'origin': 'self', 'value': 0.25}}}}

#default_mech_dict = {}


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
    The AC length constant for this section and the user-defined fraction is used to determine the maximum size
    of each segment to achieve the desired spatial and temporal resolution. This method returns the number of
    segments to set the nseg parameter for this section. For tapered cylindrical sections, the diam parameter
    will need to be reinitialized after nseg changes.
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
    Can iterate through simulations, each with some modifications of the stimulation parameters, and save all the
    recorded vectors to an HDF5 file.
    Once defined, IClamp objects persist when using an interactive console, but not when executing standalone scripts.
    Therefore, the best practice is simply to set amp to zero to turn off current injections, or move individual IClamp
    processes to different locations rather then adding and deleting them.
    class params:
    self.stim_list: list of dicts with keys for 'cell', 'node', and 'stim': pointer to hoc IClamp object
    self.rec_list: list of dicts with keys for 'cell', 'node', 'loc' and 'vec': pointer to hoc Vector object. Also
    contains keys for 'name' and 'units' for recording parameters other than Vm.
    """
    def __init__(self, tstop=400, dt=None):
        self.rec_list = []
        self.stim_list = []
        self.tstop = tstop
        if dt is None:
            self.dt = h.dt
        else:
            self.dt = dt
        self.tvec = h.Vector()
        self.tvec.record(h._ref_t)

    def run(self, v_init=-65): # these are going to have to be lists of lists so there can be more than one stim per run
        h.tstop = self.tstop
        h.dt = self.dt
        h.finitialize(v_init)
        h.run()
#        self.plot()

    def append_rec(self, cell, node, loc, param='_ref_v', ylabel='Vm', units='mV'):
        rec_dict = {'cell': cell, 'node': node, 'loc': loc, 'ylabel': ylabel, 'units': units}
        rec_dict['vec'] = h.Vector()
        rec_dict['vec'].record(getattr(node.sec(loc), param))
        self.rec_list.append(rec_dict)

    def append_stim(self, cell, node, loc, amp, delay, dur):
        stim_dict = {'cell': cell, 'node': node}
        stim_dict['stim'] = h.IClamp(node.sec(loc))
        stim_dict['stim'].amp = amp
        stim_dict['stim'].delay = delay
        stim_dict['stim'].dur = dur
        stim_dict['vec'] = h.Vector()
        stim_dict['vec'].record(stim_dict['stim']._ref_i)
        self.stim_list.append(stim_dict)

    def modify_stim(self, index=0, node=None, loc=None, amp=None, delay=None, dur=None):
        stim_dict = self.stim_list[index]
        if (node is None) | (loc is None):
            if not node is None:
                stim_dict['node'] = node
            if loc is None:
                loc = stim_dict['stim'].get_loc()
            stim_dict['stim'].loc(stim_dict['node'].sec(loc))
        if not amp is None:
            stim_dict['stim'].amp = amp
        if not delay is None:
            stim_dict['stim'].delay = delay
        if not dur is None:
            stim_dict['stim'].dur = dur

    def modify_rec(self, index=0, node=None, loc=None):
        rec_dict = self.rec_list[index]
        if not node is None:
            rec_dict['node'] = node
            rec_dict['vec'].record(node.sec(rec_dict['loc'])._ref_v)
        if not loc is None:
            rec_dict['loc'] = loc
            rec_dict['vec'].record(rec_dict['node'].sec(loc)._ref_v)

    def plot(self):
        for rec_dict in self.rec_list:
            plt.plot(self.tvec, rec_dict['vec'], label=rec_dict['node'].name+'('+str(rec_dict['loc'])+')')
            plt.xlabel("Time (ms)")
            plt.ylabel(rec_dict['ylabel']+' ('+rec_dict['units']+')')
        legend = plt.legend(loc='upper right')
        plt.show()
        plt.close()

    def export_to_file(self, f, simiter=0):
        if simiter == 0:
            f.create_dataset('time', compression='gzip', compression_opts=9, data=self.tvec)
        f.create_group(str(simiter))
        f[str(simiter)].create_group('stim')
        for index, stim in enumerate(self.stim_list):
            stim_out = f[str(simiter)]['stim'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                              data=stim['vec'])
            stim_out.attrs['cell'] = stim['cell'].gid
            stim_out.attrs['node'] = stim['node'].name
            stim_out.attrs['loc'] = stim['stim'].get_loc()
            stim_out.attrs['amp'] = stim['stim'].amp
            stim_out.attrs['delay'] = stim['stim'].delay
            stim_out.attrs['dur'] = stim['stim'].dur
        f[str(simiter)].create_group('rec')
        for index, rec in enumerate(self.rec_list):
            rec_out = f[str(simiter)]['rec'].create_dataset(str(index), compression='gzip', compression_opts=9,
                                                            data=rec['vec'])
            rec_out.attrs['cell'] = rec['cell'].gid
            rec_out.attrs['node'] = rec['node'].name
            rec_out.attrs['loc'] = rec['loc']
            rec_out.attrs['ylabel'] = rec['ylabel']
            rec_out.attrs['units'] = rec['units']


def sim_multiple(sim, stim_mods=None, rec_filename=None):
    """
    Iterate through modifications to stimuli and export all recorded vectors to unique hdf5 file.
    :param sim: :class:'QuickSim'
    :param stim_mods: list of dicts with {index: dict of kwargs} to modify stim parameters
        example: [{0: {'amp': 0.3, 'loc': 1.0}}, {0: {'amp': 0.4}}]
        would change the 1st IClamp object in the sim object to a new location and iterate for two simulations with
        different current amplitudes.
    """
    if stim_mods is None:
        stim_mods = [{}]
    if rec_filename is None:
        rec_filename = 'sim_output_'+str(datetime.datetime.today().strftime('%m%d%Y%H%M'))
    f = h5py.File(data_dir+rec_filename+'.hdf5')
    for simiter, stim_mod in enumerate(stim_mods):
        if stim_mod:  # don't change any stim parameters for empty dictionaries in list
            for index in range(len(sim.stim_list)):  # iterate through IClamp objects in sim.stim_list
                if index in stim_mod:  # apply any modifications to IClamp object if specified by dictionary
                    sim.modify_stim(index, **stim_mod[index])
        sim.run()
        sim.export_to_file(f, simiter)
    f.close()
    print 'Simulation output saved to: {}{}.hdf5'.format(data_dir, rec_filename)
    global exported_filename
    exported_filename = rec_filename


def sweep_and_plot(sim, param, min, max, inc, index=0):
    pass


def plot_from_file(rec_filename=None):
    if rec_filename is None:
        global exported_filename
        rec_filename = exported_filename
    f = h5py.File(data_dir+rec_filename+'.hdf5')
    time = f['time']
    fig, axes = plt.subplots(2, len(f)-1)
    for column in range(len(f) - 1):
        simiter = str(column)
        for rec in f[simiter]['rec'].itervalues():
            axes[0][column].plot(time, rec, label=rec.attrs['node']+'('+str(rec.attrs['loc'])+')')
        axes[0][column].set_xlabel("Time (ms)")
        axes[0][column].set_ylabel(rec.attrs['ylabel']+' ('+rec.attrs['units']+')')
        axes[0][column].legend(loc='upper right')
        for stim in f[simiter]['stim'].itervalues():
            axes[1][column].plot(time, stim, label=stim.attrs['node']+'('+str(stim.attrs['loc']))
        axes[1][column].set_xlabel("Time (ms)")
        axes[1][column].set_ylabel('Current (nA)')
        axes[1][column].legend(loc='upper right')
    plt.show()
    plt.close()
    f.close()


def plot_from_hdf(f):
    time = f['time']
    fig, axes = plt.subplots(2, len(f)-1)
    for column in range(len(f) - 1):
        simiter = str(column)
        for rec in f[simiter]['rec'].itervalues():
            axes[0][column].plot(time, rec, label=rec.attrs['node']+'('+str(rec.attrs['loc'])+')')
        axes[0][column].set_xlabel("Time (ms)")
        axes[0][column].set_ylabel(rec.attrs['ylabel']+' ('+rec.attrs['units']+')')
        axes[0][column].legend(loc='upper right')
        for stim in f[simiter]['stim'].itervalues():
            axes[1][column].plot(time, stim, label=stim.attrs['node']+'('+str(stim.attrs['loc']))
        axes[1][column].set_xlabel("Time (ms)")
        axes[1][column].set_ylabel('Current (nA)')
        axes[1][column].legend(loc='upper right')
    plt.show()
    plt.close()