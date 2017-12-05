from neuron import h
h.load_file('nrngui.hoc')
pc = h.ParallelContext()
rank = int(pc.id())
nhost = int(pc.nhost())
from cell import BallStick

class Ring(object):

  def __init__(self, ncell, delay):
    #print "construct ", self
    self.delay = delay
    self.ncell = int(ncell)
    self.mkring(self.ncell)
    self.mkstim()
    self.spike_record()

  def __del__(self):
    pc.gid_clear()
    #print "delete ", self

  def mkring(self, ncell):
    self.mkcells(ncell)
    self.connectcells(ncell)

  def mkcells(self, ncell):
    global rank, nhost
    self.cells = []
    for i in range(rank, ncell, nhost):
      cell = BallStick()
      self.cells.append(cell)
      pc.set_gid2node(i, rank)
      nc = cell.connect2target(None)
      pc.cell(i, nc)

  def connectcells(self, ncell):
    global rank, nhost
    self.nclist = []
    # not efficient but demonstrates use of pc.gid_exists
    for i in range(ncell):
      targid = (i+1)%ncell
      if pc.gid_exists(targid):
        target = pc.gid2cell(targid)
        syn = target.synlist[0]
        nc = pc.gid_connect(i, syn)
        self.nclist.append(nc)
        nc.delay = self.delay
        nc.weight[0] = 0.01


  #Instrumentation - stimulation and recording
  def mkstim(self):
    if not pc.gid_exists(0):
      return
    self.stim = h.NetStim()
    self.stim.number = 1
    self.stim.start = 0
    self.ncstim = h.NetCon(self.stim, pc.gid2cell(0).synlist[0])
    self.ncstim.delay = 0
    self.ncstim.weight[0] = 0.01

  def spike_record(self):
    self.tvec = h.Vector()
    self.idvec = h.Vector()
    for i in range(len(self.cells)):
      nc = self.cells[i].connect2target(None)
      pc.spike_record(nc.srcgid(), self.tvec, self.idvec)


def runring(ncell=5, delay=1, tstop=100):
  ring = Ring(ncell, delay)
  pc.set_maxstep(10)
  h.stdinit()
  pc.psolve(tstop)
  spkcnt = pc.allreduce(len(ring.tvec), 1)
  tmax = pc.allreduce(ring.tvec.x[-1], 2)
  tt = h.Vector()
  idv = h.Vector()
  pc.allgather(ring.tvec.x[-1], tt)
  pc.allgather(ring.idvec.x[-1], idv)
  idmax = int(idv.x[int(tt.max_ind())])
  return (int(spkcnt), tmax, idmax, (ncell, delay, tstop, (pc.id_world(), pc.nhost_world())))

"""
from neuron import h
h.load_file('nrngui.hoc')
#pc = h.ParallelContext()
from cell import BallStick

class Ring(object):

  def __init__(self, ncell, delay, pc):
    #print "construct ", self
    self.pc = pc
    self.rank = int(pc.id())
    self.nhost = int(pc.nhost())
    self.delay = delay
    self.ncell = int(ncell)
    self.mkring(self.ncell, self.rank, self.nhost)
    self.mkstim()
    self.spike_record()


  def __del__(self):
    self.pc.gid_clear()
    #print "delete ", self

  def mkring(self, ncell, rank, nhost):
    self.mkcells(ncell, rank, nhost)
    self.connectcells(ncell)

  def mkcells(self, ncell, rank, nhost):
    self.cells = []
    for i in range(rank, ncell, nhost):
      cell = BallStick()
      self.cells.append(cell)
      self.pc.set_gid2node(i, rank)
      nc = cell.connect2target(None)
      self.pc.cell(i, nc)

  def connectcells(self, ncell):
    self.nclist = []
    # not efficient but demonstrates use of pc.gid_exists
    for i in range(ncell):
      targid = (i+1)%ncell
      if self.pc.gid_exists(targid):
        target = self.pc.gid2cell(targid)
        syn = target.synlist[0]
        nc = self.pc.gid_connect(i, syn)
        self.nclist.append(nc)
        nc.delay = self.delay
        nc.weight[0] = 0.01


  #Instrumentation - stimulation and recording
  def mkstim(self):
    if not self.pc.gid_exists(0):
      return
    self.stim = h.NetStim()
    self.stim.number = 1
    self.stim.start = 0
    self.ncstim = h.NetCon(self.stim, self.pc.gid2cell(0).synlist[0])
    self.ncstim.delay = 0
    self.ncstim.weight[0] = 0.01

  def spike_record(self):
    self.tvec = h.Vector()
    self.idvec = h.Vector()
    for i in range(len(self.cells)):
      nc = self.cells[i].connect2target(None)
      self.pc.spike_record(nc.srcgid(), self.tvec, self.idvec)


def runring(ncell=5, delay=1, tstop=100):
  ring = Ring(ncell, delay)
  ring.pc.set_maxstep(10)
  h.stdinit()
  ring.pc.psolve(tstop)
  spkcnt = ring.pc.allreduce(len(ring.tvec), 1)
  tmax = ring.pc.allreduce(ring.tvec.x[-1], 2)
  tt = h.Vector()
  idv = h.Vector()
  ring.pc.allgather(ring.tvec.x[-1], tt)
  ring.pc.allgather(ring.idvec.x[-1], idv)
  idmax = int(idv.x[int(tt.max_ind())])
  return (int(spkcnt), tmax, idmax, (ncell, delay, tstop, (ring.pc.id_world(), ring.pc.nhost_world())))
"""
