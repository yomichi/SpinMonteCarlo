import numpy as np
import numpy.random as random

execfile('observable.py')
execfile('square.py')


class Wolff:
  def __init__(self, L, T):
  ## parameters
  self.lattice = Square(L)
  self.L = L
  self.N = L*L
  self.T = T
  self.beta = 1/T

  ## configurations
  self.spins = np.ones(self.N, dtype=np.int)
  self.mag = self.num_site()
  self.ene = -self.num_bond()

  def num_site(self): return self.lattice.num_site()
  def num_bond(self): return self.lattice.num_bond()
  def neighbors(self, site): return self.lattice(site)

  def flip(self, site): self.spins[site] *= -1
  
  def update(self):
  stack = []
  site = random.randint(self.N)
  centerspin = self.spins[site]
  self.flip(site)
  stack.append(site)
  csize = 1
  while len(stack) > 0:
    site = stack.pop()
    for neighbor in self.neighbors(site):
    if self.spins[neighbor] == centerspin and random.rand() < self.prob:
      self.flip(neighbor)
      stack.append(neighbor)
      csize += 1
  improved_mag = self.mag - centerspin*csize
  self.mag -= 2*centerspin*csize
  return improved_mag

  def run(self, MCS = 8192, Thermalize=int(MCS/8)):
  for mcs in range(Thermalize):
    self.update()
  mag = Observable()
  for mcs in range(MCS):
    m = self.update()
    mag.add(m)
