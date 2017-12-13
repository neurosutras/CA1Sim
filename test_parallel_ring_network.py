# Source: https://neuron.yale.edu/neuron/static/docs/neuronpython/ballandstick5.html
# Also reference https://www.neuron.yale.edu/neuron/static/papers/jnm/parallelizing_models_jnm2008.pdf (published article)
#    and https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=96444&file=/HinesCarnevaleJNM2007/ring/ringpar.hoc#tabs-2
#        (hoc code)
# Parallel context documentation: https://www.neuron.yale.edu/neuron/static/docs/help/neuron/neuron/classes/parcon.html

# To run: mpiexec -n 1 python test_parallel_ring_network.py
# Must also intall neuronpy: pip install neuronpy

import os
from mpi4py import MPI
from neuron import h
import numpy
from itertools import izip
from neuronpy.util import spiketrain
from math import sin, cos, pi
from neuronpy.graphics import spikeplot

h.load_file('stdrun.hoc')

class Cell(object):
    """Generic cell template."""
    def __init__(self):
        self.x, self.y, self.z = 0, 0, 0
        self.synlist = [] #### NEW CONSTRUCT IN THIS WORKSHEET
        self.create_sections()
        self.build_topology()
        self.build_subsets()
        self.define_geometry()
        self.define_biophysics()
        self.create_synapses()
    #
    def create_sections(self):
        """Create the sections of the cell. Remember to do this
        in the form::
            h.Section(name='soma', cell=self)
        """
        raise NotImplementedError("create_sections() is not implemented.")
    #
    def build_topology(self):
        """Connect the sections of the cell to build a tree."""
        raise NotImplementedError("build_topology() is not implemented.")
    #
    def define_geometry(self):
        """Set the 3D geometry of the cell."""
        raise NotImplementedError("define_geometry() is not implemented.")
    #
    def define_biophysics(self):
        """Assign the membrane properties across the cell."""
        raise NotImplementedError("define_biophysics() is not implemented.")
    #
    def create_synapses(self):
        """Subclasses should create synapses (such as ExpSyn) at various
        segments and add them to self.synlist."""
        pass # Ignore if child does not implement.
    #
    def build_subsets(self):
        """Build subset lists. This defines 'all', but subclasses may
        want to define others. If overridden, call super() to include 'all'."""
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)
    #
    def connect2target(self, target, thresh=10):
        """Make a new NetCon with this cell's membrane
        potential at the soma as the source (i.e. the spike detector)
        onto the target passed in (i.e. a synapse on a cell).
        Subclasses may override with other spike detectors."""
        nc = h.NetCon(self.soma(1)._ref_v, target, sec = self.soma)
        nc.threshold = thresh
        return nc
    #
    def is_art(self):
        """Flag to check if we are an integrate-and-fire artificial cell."""
        return 0
    #
    def set_position(self, x, y, z):
        """
        Set the base location in 3D and move all other
        parts of the cell relative to that location.
        """
        for sec in self.all:
            for i in range(int(h.n3d())):
                h.pt3dchange(i,
                        x - self.x + h.x3d(i),
                        y - self.y + h.y3d(i),
                        z - self.z + h.z3d(i),
                        h.diam3d(i))
        self.x, self.y, self.z = x, y, z
    #
    def rotateZ(self, theta):
        """Rotate the cell about the Z axis."""
        for sec in self.all:
            for i in range(int(h.n3d(sec=sec))):
                x = h.x3d(i, sec=sec)
                y = h.y3d(i, sec=sec)
                c = cos(theta)
                s = sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                h.pt3dchange(i, xprime, yprime, h.z3d(i, sec=sec), h.diam3d(i, sec=sec), sec=sec)


class BallAndStick(Cell):  #### Inherits from Cell
    """Two-section cell: A soma with active channels and
    a dendrite with passive properties."""
    #### __init__ is gone and handled in Cell.
    #### We can override __init__ completely, or do some of
    #### our own initialization first, and then let Cell do its
    #### thing, and then do a bit more ourselves with "super".
    ####
    #### def __init__(self):
    ####     # Do some stuff
    ####     super(Cell, self).__init__()
    ####     # Do some more stuff
    #
    def create_sections(self):
        """Create the sections of the cell."""
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
    #
    def build_topology(self):
        """Connect the sections of the cell to build a tree."""
        self.dend.connect(self.soma(1))
    #
    def define_geometry(self):
        """Set the 3D geometry of the cell."""
        self.soma.L = self.soma.diam = 12.6157 # microns
        self.dend.L = 200                      # microns
        self.dend.diam = 1                     # microns
        self.dend.nseg = 5
        self.shape_3D()
    #
    def define_biophysics(self):
        """Assign the membrane properties across the cell."""
        for sec in self.all: # 'all' exists in parent object.
            sec.Ra = 100    # Axial resistance in Ohm * cm
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
        # Insert active Hodgkin-Huxley current in the soma
        self.soma.insert('hh')
        self.soma.gnabar_hh = 0.12  # Sodium conductance in S/cm2
        self.soma.gkbar_hh = 0.036  # Potassium conductance in S/cm2
        self.soma.gl_hh = 0.0003    # Leak conductance in S/cm2
        self.soma.el_hh = -54.3     # Reversal potential in mV
        # Insert passive current in the dendrite
        self.dend.insert('pas')
        self.dend.g_pas = 0.001  # Passive conductance in S/cm2
        self.dend.e_pas = -65    # Leak reversal potential mV
    #
    def shape_3D(self):
        """
        Set the default shape of the cell in 3D coordinates.
        Set soma(0) to the origin (0,0,0) and dend extending along
        the X-axis.
        """
        len1 = self.soma.L
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(0, 0, 0, self.soma.diam, sec=self.soma)
        h.pt3dadd(len1, 0, 0, self.soma.diam, sec=self.soma)
        len2 = self.dend.L
        h.pt3dclear(sec=self.dend)
        h.pt3dadd(len1, 0, 0, self.dend.diam, sec=self.dend)
        h.pt3dadd(len1 + len2, 0, 0, self.dend.diam, sec=self.dend)
    #
    #### build_subsets, rotateZ, and set_location are gone. ####
    #
    #### NEW STUFF ####
    #
    def create_synapses(self):
        """Add an exponentially decaying synapse in the middle
        of the dendrite. Set its tau to 2ms, and append this
        synapse to the synlist of the cell."""
        syn = h.ExpSyn(self.dend(0.5))
        syn.tau = 2
        self.synlist.append(syn) # synlist is defined in Cell


class Ring:
    """A network of *N* ball-and-stick cells where cell n makes an
    excitatory synapse onto cell n + 1 and the last, Nth cell in the
    network projects to the first cell.
    """
    def __init__(self, N=5, stim_w=0.004, stim_spike_num=1, syn_w=0.01,
            syn_delay=5):
        """
        :param N: Number of cells.
        :param stim_w: Weight of the stimulus
        :param stim_spike_num: Number of spikes generated in the stimulus
        :param syn_w: Synaptic weight
        :param syn_delay: Delay of the synapse
        """
        self._N = N              # Total number of cells in the net
        self.cells = []          # Cells on this host
        self.nclist = []         # NetCon list on this host
        self.gidlist = []        # List of global identifiers on this host
        self.stim = None         # Stimulator
        self.stim_w = stim_w     # Weight of stim
        self.stim_spike_num = stim_spike_num  # Number of stim spikes
        self.syn_w = syn_w       # Synaptic weight
        self.syn_delay = syn_delay  # Synaptic delay
        self.t_vec = h.Vector()   # Spike time of all cells on this host
        self.id_vec = h.Vector()  # Ids of spike times on this host

        #### Make a new ParallelContext object
        self.pc = h.ParallelContext()

        self.set_numcells(N)  # Actually build the net -- at least the portion
                              # of cells on this host.

    def set_numcells(self, N, radius=50):
        """Create, layout, and connect N cells."""
        self._N = N
        self.set_gids() #### Used when creating and connecting cells
        self.create_cells()
        self.connect_cells()
        self.connect_stim()

    def set_gids(self):
        """Set the gidlist on this host."""
        self.gidlist = []
        #### Round-robin counting.
        #### Each host as an id from 0 to pc.nhost() - 1.
        for i in range(int(self.pc.id()), self._N, int(self.pc.nhost())):
            self.gidlist.append(i)

    def create_cells(self):
        """Create cell objects on this host and set their location."""
        self.cells = []
        N = self._N
        r = 50 # Radius of cell locations from origin (0,0,0) in microns

        for i in self.gidlist:
            cell = BallAndStick()
            # When cells are created, the soma location is at (0,0,0) and
            # the dendrite extends along the X-axis.
            # First, at the origin, rotate about Z.
            cell.rotateZ(i*2*numpy.pi/N)

            # Then reposition
            x_loc = float(numpy.cos(i*2*numpy.pi/N))*r
            y_loc = float(numpy.sin(i*2*numpy.pi/N))*r
            cell.set_position(x_loc, y_loc, 0)

            self.cells.append(cell)

            #### Tell this host it has this gid
            #### gids can be any integer, they just need to be unique.
            #### In this simple case, we set the gid to i.
            self.pc.set_gid2node(i, int(self.pc.id()))

            #### Means to tell the ParallelContext that this cell is
            #### a source for all other hosts. NetCon is temporary.
            nc = cell.connect2target(None)
            self.pc.cell(i, nc) # Associate the cell with this host and gid

            #### Record spikes of this cell
            self.pc.spike_record(i, self.t_vec, self.id_vec)

# OLD WAY
#    def connect_cells(self):
#        self.nclist = []
#        N = self._N
#        for i in range(N):
#            src = self.cells[i]
#            tgt_syn = self.cells[(i+1)%N].synlist[0]
#            nc = src.connect2target(tgt_syn)
#            nc.weight[0] = self.syn_w
#            nc.delay = self.syn_delay
#            nc.record(self.t_vec, self.id_vec, i)
#            self.nclist.append(nc)

    def connect_cells(self):
        """Connect cell n to cell n + 1."""
        self.nclist = []
        N = self._N
        for i in self.gidlist:
            src_gid = (i-1+N) % N
            tgt_gid = i
            if self.pc.gid_exists(tgt_gid):
                target = self.pc.gid2cell(tgt_gid)
                syn = target.synlist[0]
                nc = self.pc.gid_connect(src_gid, syn)
                nc.weight[0] = self.syn_w
                nc.delay = self.syn_delay
                self.nclist.append(nc)

    def connect_stim(self):
        """Connect a spike generator on the first cell in the network."""
        #### If the first cell is not on this host, return
        if not self.pc.gid_exists(0):
            return
        self.stim = h.NetStim()
        self.stim.number = self.stim_spike_num
        self.stim.start = 9
        self.ncstim = h.NetCon(self.stim, self.cells[0].synlist[0])
        self.ncstim.delay = 1
        self.ncstim.weight[0] = self.stim_w # NetCon weight is a vector.

    def get_spikes(self):
        """Get the spikes as a list of lists."""
        return spiketrain.netconvecs_to_listoflists(self.t_vec, self.id_vec)

    def write_spikes(self, file_name='out.spk'):
        """Append the spike output file with spikes on this host. The output
        format is the timestamp followed by a tab then the gid of the source
        followed by a newline.

        :param file_name: is the full or relative path to a spike output file.

        .. note::

            When parallelized, each process will write to the same file so it
            is opened in append mode. The order in which the processes write is
            arbitrary so while the spikes within the process may be ordered by
            time, the output file will be unsorted. A quick way to sort a file
            is with the bash command sort, which can be called after all
            processes have written the file with the following format::

                exec_cmd = 'sort -k 1n,1n -k 2n,2n ' + file_name + \
                        ' > ' + 'sorted_' + file_name
                os.system(exec_cmd)
        """
        for i in range(int(self.pc.nhost())):
            self.pc.barrier() # Sync all processes at this point
            if i == int(self.pc.id()):
                if i == 0:
                    mode = 'w' # write
                else:
                    mode = 'a' # append
                with open(file_name, mode) as spk_file: # Append
                    for (t, id) in izip(self.t_vec, self.id_vec):
                        spk_file.write('%.3f\t%d\n' %(t, id)) # timestamp, id
        self.pc.barrier()


def prun():
    pc = h.ParallelContext()
    my_ring = Ring()
    pc.set_maxstep(10)
    h.stdinit()
    h.dt = 0.025 # Fixed dt
    print 'world %d of %d,  bbs %d of %d,  net %d of %d' %(int(pc.id_world()), int(pc.nhost_world()), int(pc.id_bbs()),
                                                         int(pc.nhost_bbs()), int(pc.id()), int(pc.nhost()))
    pc.psolve(100)
    my_ring.write_spikes('spike_out_file')
    pc.runworker()
    exec_cmd = 'sort -k 1n,1n -k 2n,2n spike_out_file > sorted_spike_out_file'
    os.system(exec_cmd)
    """
    #This is possible if script is run on one rank, so that there is only one ring object.
    spikes = my_ring.get_spikes()
    """
    f = open('sorted_spike_out_file', 'r')
    lines = f.readlines()
    spike_dict = {}
    for line in lines:
        time, id = line.split('\t')
        time = float(time)
        id = int(id)
        if id not in spike_dict:
            spike_dict[id] = []
        spike_dict[id].append(time)
    spikes = spike_dict.values()
    sp = spikeplot.SpikePlot()
    sp.plot_spikes(spikes)



prun()
