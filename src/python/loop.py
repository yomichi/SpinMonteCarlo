import random
import numpy

from observable import *
from unionfind  import *

class GraphElement:
    def __init__(self, bond, time, bottom=0, top=0):
        self.bond = bond
        self.time = time
        self.bottom = bottom  # index of bottom node
        self.top = top        # index of top node
        self.isdiagonal = True

    def flip(self): self.isdiagonal ^= True

class Looper:
    """
    Loop algorithm for S=1/2 antiferromagnetic Heisenberg periodic chain with length L at temperature T
    """
    def __init__(self, L, T, J=1.0):

        ## parameter
        self.L = L
        self.beta = 1.0/T
        self.J = J  ## coupling constant

        ## configuration

        # 1 (-1) denotes up (down) spin
        self.spins = numpy.ones(L, dtype=numpy.int)

        self.operators = []
        self.UF = UnionFind(L)
        self.numcluster = 0
        self.bj2 = 0.5*L*self.beta*self.J

    def leftsite(self, bond):
        return bond

    def rightsite(self, bond):
        return (bond + 1)%self.L

    # delegete to UnionFind
    def add(self, index): return self.UF.add(index)
    def union(self, index0, index1): return self.UF.union(index0, index1)
    def clusterize(self): return self.UF.clusterize()
    def cluster_id(self, index): return self.UF.cluster_id(index)

    def update(self):
        nodes = range(self.L)  # union find node for each site
        self.UF.reset(self.L)

        # diagonal update -- scattering graph elements which connect anti-parallel spins
        #
        # |  |       |  |
        # |  |  <->  ----
        # |  |       ----
        # |  |       |  |
        #
        # weight of this graph is
        # 
        # w = \beta J / 2  for anti-parallel spins
        #     0            for parallel spins
        #

        oi = 0
        t = random.expovariate(self.bj2)
        while t < 1.0 or oi < len(self.operators) :
            if oi == len(self.operators) or t < self.operators[oi].time :
                # insert new graph element at t

                b = numpy.random.randint(self.L)
                lsite, rsite = self.leftsite(b), self.rightsite(b)
                if self.spins[lsite] != self.spins[rsite] :
                    newop = GraphElement(bond = b, time = t)
                    self.operators.insert(oi, newop)
                    t += random.expovariate(self.bj2)
                else :
                    t += random.expovariate(self.bj2)
                    continue
            else:
                # remove passed graph element

                if self.operators[oi].isdiagonal :
                    self.operators.remove(self.operators[oi])
                    continue
            b = self.operators[oi].bond
            lsite, rsite = self.leftsite(b), self.rightsite(b)
            bottom = self.UF.union(nodes[lsite], nodes[rsite])
            top = self.UF.add()
            self.operators[oi].bottom = bottom
            self.operators[oi].top = nodes[lsite] = nodes[rsite] = top
            if not(self.operators[oi].isdiagonal) :
                self.spins[lsite] *= -1
                self.spins[rsite] *= -1
            oi += 1
        # end of diagonal update
        
        # boundary condition along time is periodic
        for site in range(self.L):
            self.UF.union(site, nodes[site])
        
        self.numcluster = self.clusterize()
        flips = numpy.random.randint(2, size = self.numcluster)
        for op in self.operators:
            bID, tID = self.UF.cluster_id(op.bottom), self.UF.cluster_id(op.top)
            # print(len(flips), bID, tID)
            if flips[bID] ^ flips[tID] :
                op.flip()
        for site in range(self.L):
            if flips[self.UF.cluster_id(self.UF.root(site))]:
                self.spins[site] *= -1

    def measure(self) :
        """
        measure uniform susceptibility, staggered susceptibility, and energy with improved estimator
        """
        mags = numpy.zeros(self.numcluster)  ## magnetization of loops => uni. suscep.
        lens = numpy.zeros(self.numcluster)  ## length of loops => stag. suscep.
        nodes = range(self.L)
        for op in self.operators:
            bottom, top = self.UF.cluster_id(op.bottom), self.UF.cluster_id(op.top)
            lens[bottom] += 2 * op.time
            lens[top] -= 2 * op.time
        for site in range(self.L):
            sid = self.UF.cluster_id(site)
            lens[sid] += 1.0
            mags[sid] += 0.5*self.spins[site]
        m2 = mags.dot(mags)
        l2 = lens.dot(lens)
        m2 *= self.beta/self.L
        l2 *= 0.25*self.beta/self.L
        ene = 0.25*self.L*self.J - len(self.operators)/self.beta
        return m2, l2

    def run(self, MCS=8192, Thermalization=1024):
        for mcs in range(Thermalization):
            self.update()
        usus = Observable()
        ssus = Observable()
        ene = Observable()
        for mcs in range(MCS):
            self.update()
            m2,l2,e = self.measure()
            usus.add(m2)
            ssus.add(l2)
            ene.add(e)

        return usus, ssus, ene
