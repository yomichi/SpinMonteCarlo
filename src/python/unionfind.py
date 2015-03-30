class UnionFind:
    def __init__(self, num=0):
        self.table = range(num)
        self.cluster_IDs = None

    def reset(self, num=0):
        self.table = range(num)
        self.cluster_IDs = None
        return self

    def num(self):
        return len(self.table)
    
    def isroot(self, index):
        return self.table[index] == index

    def root(self, index):
        while not(self.isroot(index)):
            index = self.table[index]
        return index

    def add(self):
        self.cluster_IDs = None
        N = self.num()
        self.table.append(N)
        return N

    def union(self, index0, index1):
        if index0 == index1 :
            return index0
        self.cluster_IDs = None
        level0, level1 = 0,0
        rt0,rt1 = index0, index1
        while not(self.isroot(rt0)):
            level0 += 1
            rt0 = self.table[rt0]
        while not(self.isroot(rt1)):
            level1 += 1
            rt1 = self.table[rt1]
        if level0 < level1 : 
            self.table[rt0] = self.table[index0] = self.table[index1] = rt1
            return rt1
        else:
            self.table[rt1] = self.table[index0] = self.table[index1] = rt0
            return rt0

    def clusterize(self):
        N = 0
        self.cluster_IDs = range(self.num())
        for index in xrange(self.num()):
            if self.isroot(index):
                self.cluster_IDs[index] = N
                N += 1
        for index in xrange(self.num()):
            rt = self.root(index)
            self.cluster_IDs[index] = self.cluster_IDs[rt]
        return N

    def cluster_id(self, index):
        if self.cluster_IDs == None:
            self.clusterize()
        return self.cluster_IDs[index]
