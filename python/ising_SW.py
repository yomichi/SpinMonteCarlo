import math
import numpy as np
import numpy.random as random

execfile('observable.py')
execfile('unionfind.py')
execfile('square.py')

class SW:
  def __init__(self, L, T):
    ## parameters
    self.L = L
    self.lattice = Square(L)
    self.T = T
    self.beta = 1.0/T

    ## configurations
    self.spins = np.ones(self.num_site(), dtype=np.int)
    self.mag = self.num_site()
    self.ene = -self.num_bond()

    ## algorithm parameter
    self.prob = 1.0-math.exp(-2*self.beta)
    self.UF = UnionFind(self.num_site())

  ## delegate to lattice
  def num_site(self) : return self.lattice.num_site()
  def num_bond(self) : return self.lattice.num_bond()
  def neighbors(self, site) : return self.lattice.neighbors(site)
  def site2coord(self, site) : return self.lattice.site2coord(site)
  def coord2site(self, x,y) : return self.lattice.coord2site(x,y)

  ## delegate to UnionFind
  def unify(self, lsite, rsite) : return self.UF.unify(lsite, rsite)
  def cluster_id(self, site) : return self.UF.cluster_id(site)

  def update(self):
    self.UF.reset(self.num_site())

    num_activated = 0

    for site in xrange(self.num_site()):
      x,y = self.site2coord(site)
      spin = self.spins[site]
      for neighbor in self.neighbors(site)[0:2]:
        if self.spins[neighbor] == spin and random.rand() < self.prob :
          self.unify(site, neighbor)
          num_activated += 1
    num_clusters = self.UF.clusterize()
    cl_spins = random.randint(2, size=num_clusters)*2-1
    cl_sizes = np.zeros(num_clusters)

    for site in xrange(self.num_site()):
      cid = self.cluster_id(site)
      self.spins[site] = cl_spins[cid]
      cl_sizes[cid] += 1.0

    am = np.sum(self.spins)
    am *= 1.0/self.num_site()
    am = abs(am)
    cl_sizes *= 1.0/self.num_site()
    cl_sizes *= cl_sizes
    m2 = np.sum(cl_sizes)
    cs = np.cumsum(cl_sizes)
    m4 = 12*np.dot(cl_sizes[1:], cs[:-1]) + np.dot(cl_sizes,cl_sizes)

    ene = self.num_bond() - (1.0+1.0/math.tanh(self.beta))*num_activated
    ene /= self.num_bond()

    return am, m2, m4, ene

  def run(self, MCS=8192, Thermalization=1024):
    for mcs in xrange(Thermalization):
      self.update()
    absmag = Observable()
    mag2 = Observable()
    mag4 = Observable()
    sus = Observable()
    ene = Observable()
    for mcs in xrange(MCS):
      am,m2, m4, en = self.update()
      chi = m2 * self.beta/self.num_site()
      absmag.add(am)
      mag2.add(m2)
      mag4.add(m4)
      sus.add(chi)
      ene.add(en)
    return absmag, mag2, mag4, sus, ene
