class UnionFind:
  def __init__(self, num=0):
    self.parents = range(num)
    self.weights = [1]*num
    self.cluster_IDs = None

  def reset(self, num=0):
    self.parents = range(num)
    self.weights = [1] * num
    self.cluster_IDs = None
    return self

  def num(self):
    return len(self.parents)
  
  def isroot(self, index):
    return self.parents[index] == index

  def root(self, index):
    while not(self.isroot(index)):
      index = self.parents[index]
    return index

  def add(self):
    self.cluster_IDs = None
    N = self.num()
    self.parents.append(N)
    self.weights.append(1)
    return N

  def unify(self, index0, index1):
    if index0 == index1 :
      return index0
    self.cluster_IDs = None
    rt0,rt1 = self.root(index0), self.root(index1)
    w0,w1 = self.weights[rt0], self.weights[rt1]
    if w0 < w1 :
      rt0,rt1 = rt1,rt0

    self.parents[rt1] = rt0
    if w0 == w1:
      self.weights[rt0] += 1
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
