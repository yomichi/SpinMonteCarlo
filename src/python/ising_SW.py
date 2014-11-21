import math
import numpy
import numpy.random as random

from observable import *
from unionfind import *
from square import *

class Cluster:
    def __init__(self):
        self.spin = 2*random.randint(2)-1
        self.mag = 0

class SW:
    def __init__(self, L, T):
        ## parameters
        self.L = L
        self.lattice = Square(L)
        self.T = T
        self.beta = 1.0/T

        ## configurations
        self.spins = numpy.ones(self.num_site(), dtype=numpy.int)
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
    def union(self, lsite, rsite) : return self.UF.union(lsite, rsite)
    def cluster_id(self, site) : return self.UF.cluster_id(site)

    def update(self):
        self.UF.reset(self.num_site())

        num_activated = 0
        for site in range(self.num_site()):
            x,y = self.site2coord(site)
            spin = self.spins[site]
            for neighbor in [self.coord2site((x+1)%self.L,y), self.coord2site(x,(y+1)%self.L)]:
                if self.spins[neighbor] == spin and random.rand() < self.prob :
                    self.union(site, neighbor)
                    num_activated += 1
        self.numcluster = self.UF.clusterize()
        self.clusters = [ Cluster() for ic in range(self.numcluster)]

        for site in range(self.num_site()):
            cid = self.cluster_id(site)
            self.spins[site] = self.clusters[cid].spin
            self.clusters[cid].mag += 1

        m2 = 0.0
        for cluster in self.clusters:
            m2 += cluster.mag**2

        ene = self.num_bond() - (1.0+1.0/math.tanh(self.beta))*num_activated

        return m2, ene

    def run(self, MCS=8192, Thermalization=1024):
        for mcs in range(Thermalization):
            self.update()
        mag2 = Observable()
        ene = Observable()
        for mcs in range(MCS):
            m2, en = self.update()
            m2 *= self.beta/self.num_site()
            mag2.add(m2)
            ene.add(en)
        return mag2, ene
