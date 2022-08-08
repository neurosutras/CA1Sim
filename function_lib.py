__author__ = 'Aaron D. Milstein'
from mpi4py import MPI
from collections import Iterable
import h5py
import math
import datetime
import copy
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import random
import pprint
import sys
import os
import yaml


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
                     'axon_hill': {'cable': {'Ra': {'origin': 'soma'}, 'cm': {'origin': 'soma'}},
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
    diam = np.mean([seg.diam for seg in sec])
    Ra = sec.Ra
    cm = np.mean([seg.cm for seg in sec])
    return 1e5*math.sqrt(diam/(4.*math.pi*f*Ra*cm))


def d_lambda_nseg(sec, lam=d_lambda, f=freq):
    """
    The AC length constant for this section and the user-defined fraction is used to determine the maximum size of each
    segment to achieve the d esired spatial and temporal resolution. This method returns the number of segments to set
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

    print('Saving SWC to file '+filenameBase+'-scaled.swc')
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
    print('x - ',min(xvals),':',max(xvals))
    print('y - ',min(yvals),':',max(yvals))
    print('z - ',min(zvals),':',max(zvals))
    print('r - ',min(rvals),':',max(rvals))


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
    indexes = list(range(len(z_offsets)))
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
        print('Saving '+filenames[i]+' to '+outputname)
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


class ExplicitDumper(yaml.SafeDumper):
    """
    YAML dumper that will never emit aliases.
    """

    def ignore_aliases(self, data):
        return True


def nested_convert_scalars(data):
    """
    Crawls a nested dictionary, and converts any scalar objects from numpy types to python types.
    :param data: dict
    :return: dict
    """
    if isinstance(data, dict):
        for key in data:
            data[key] = nested_convert_scalars(data[key])
    elif isinstance(data, Iterable) and not isinstance(data, (str, tuple)):
        data = list(data)
        for i in range(len(data)):
            data[i] = nested_convert_scalars(data[i])
    elif hasattr(data, 'item'):
        data = data.item()
    return data


def write_to_yaml(file_path, data, convert_scalars=True):
    """

    :param file_path: str (should end in '.yaml')
    :param data: dict
    :param convert_scalars: bool
    :return:
    """
    with open(file_path, 'w') as outfile:
        if convert_scalars:
            data = nested_convert_scalars(data)
        yaml.dump(data, outfile, default_flow_style=False, Dumper=ExplicitDumper)


def read_from_yaml(file_path, include_loader=None):
    """

    :param file_path: str (should end in '.yaml')
    :return:
    """
    if os.path.isfile(file_path):
        with open(file_path, 'r') as stream:
            if include_loader is None:
                Loader = yaml.FullLoader
            else:
                Loader = include_loader
            data = yaml.load(stream, Loader=Loader)
        return data
    else:
        raise IOError('read_from_yaml: invalid file_path: %s' % file_path)


def combine_hdf5_file_paths(file_path_list, new_file_path=None):
    """
    List contains names of files generated by "embarassingly parallel" execution of simulations on separate cores.
    This function combines the contents of the files into one .hdf5 file.
    :param file_path_list: list of str (paths)
    :param new_file_path: str (path)
    """
    if new_file_path is None:
        raise ValueError('combine_output_file_paths: invalid file path provided: %s' % new_file_path)
    new_f = h5py.File(new_file_path, 'w')
    iter = 0
    for old_file_path in file_path_list:
        old_f = h5py.File(old_file_path, 'r')
        for old_group in old_f.values():
            new_f.copy(old_group, new_f, name=str(iter))
            iter += 1
        old_f.close()
    new_f.close()
    print('combine_output_file_paths: exported to file path: %s' % new_file_path)


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

def interpolate_tvec_vec(tvec, vec, duration, dt=0.02):
    """
    Interpolates the array for tvec from t=0 to t=duration according to the dt time step, and interpolates the
    vec array to correlate with the interpolated tvec array.
    :param tvec: vector of times
    :param vec: vector of voltage recordings
    :param duration: length of time of voltage trace
    :param dt:
    :return:
    """
    interp_t = np.arange(0., duration, dt)
    interp_vm = np.interp(interp_t, tvec, vec)
    return interp_t, interp_vm



def get_instantaneous_spike_probability(rate, dt=0.02, generator=None):
    """
    Given an instantaneous spike rate in Hz and a time bin in ms, calculate whether a spike has occurred in that time
    bin, assuming an exponential distribution of inter-spike intervals.
    :param rate: float (Hz)
    :param dt: float (ms)
    :param generator: :class:'random.Random'
    :return: bool
    """
    if generator is None:
        generator = random
    x = generator.uniform(0., 1.)
    rate /= 1000.
    p = 1. - np.exp(-rate * dt)
    return bool(x < p)


def get_inhom_poisson_spike_times_by_thinning(rate, t, dt=0.02, refractory=3., generator=None):
    """
    Given a time series of instantaneous spike rates in Hz, produce a spike train consistent with an inhomogeneous
    Poisson process with a refractory period after each spike.
    :param rate: instantaneous rates in time (Hz)
    :param t: corresponding time values (ms)
    :param dt: temporal resolution for spike times (ms)
    :param refractory: absolute deadtime following a spike (ms)
    :param generator: :class:'random.Random()'
    :return: list of m spike times (ms)
    """
    if generator is None:
        generator = random
    interp_t = np.arange(t[0], t[-1] + dt, dt)
    interp_rate = np.interp(interp_t, t, rate)
    interp_rate /= 1000.
    non_zero = np.where(interp_rate > 0.)[0]
    interp_rate[non_zero] = 1. / (1. / interp_rate[non_zero] - refractory)
    spike_times = []
    max_rate = np.max(interp_rate)
    i = 0
    ISI_memory = 0.
    while i < len(interp_t):
        x = generator.random()
        if x > 0.:
            ISI = -np.log(x) / max_rate
            i += int(ISI / dt)
            ISI_memory += ISI
            if (i < len(interp_t)) and (generator.random() <= interp_rate[i] / max_rate) and ISI_memory >= 0.:
                spike_times.append(interp_t[i])
                ISI_memory = -refractory
    return spike_times


def get_patterned_input_component_traces(rec_filename, dt=0.02):
    """

    :param rec_file_name: str
    # remember .attrs['phase_offset'] could be inside ['train'] for old files
    """
    with h5py.File(data_dir+rec_filename+'.hdf5', 'r') as f:
        sim = next(iter(f.values()))
        equilibrate = sim.attrs['equilibrate']
        track_equilibrate = sim.attrs['track_equilibrate']
        duration = sim.attrs['duration']
        track_duration = duration - equilibrate - track_equilibrate
        start = int((equilibrate + track_equilibrate)/dt)
        vm_array = []
        for sim in f.values():
            t = np.arange(0., duration, dt)
            vm = np.interp(t, sim['time'], sim['rec']['0'])
            vm = vm[start:]
            vm_array.append(vm)
    rec_t = np.arange(0., track_duration, dt)
    spikes_removed = get_removed_spikes(rec_filename, plot=0, dt=dt)
    # down_sample traces to 2 kHz after clipping spikes for theta and ramp filtering
    down_dt = 0.5
    down_t = np.arange(0., track_duration, down_dt)
    # 2000 ms Hamming window, ~2 Hz low-pass for ramp, ~5 - 10 Hz bandpass for theta
    window_len = int(2000./down_dt)
    pad_len = int(window_len/2.)
    theta_filter = signal.firwin(window_len, [5., 10.], nyq=1000./2./down_dt, pass_zero=False)
    ramp_filter = signal.firwin(window_len, 2., nyq=1000./2./down_dt)
    theta_traces = []
    ramp_traces = []
    for trace in spikes_removed:
        down_sampled = np.interp(down_t, rec_t, trace)
        padded_trace = np.zeros(len(down_sampled)+window_len)
        padded_trace[pad_len:-pad_len] = down_sampled
        padded_trace[:pad_len] = down_sampled[::-1][-pad_len:]
        padded_trace[-pad_len:] = down_sampled[::-1][:pad_len]
        filtered = signal.filtfilt(theta_filter, [1.], padded_trace, padlen=pad_len)
        filtered = filtered[pad_len:-pad_len]
        up_sampled = np.interp(rec_t, down_t, filtered)
        theta_traces.append(up_sampled)
        filtered = signal.filtfilt(ramp_filter, [1.], padded_trace, padlen=pad_len)
        filtered = filtered[pad_len:-pad_len]
        up_sampled = np.interp(rec_t, down_t, filtered)
        ramp_traces.append(up_sampled)
    return rec_t, vm_array, theta_traces, ramp_traces


def clean_axes(axes):
    """
    Remove top and right axes from pyplot axes object.
    :param axes:
    """
    if not type(axes) in [np.ndarray, list]:
        axes = [axes]
    elif type(axes) == np.ndarray:
        axes = axes.flatten()
    for axis in axes:
        axis.tick_params(direction='out')
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
        axis.get_xaxis().tick_bottom()
        axis.get_yaxis().tick_left()

