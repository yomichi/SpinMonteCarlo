import math

class Observable:
  """
  Simple Monte Calro mean and error estimator

  [NOTICE] No error check will be performed for simplicity :p
  """
  def __init__(self):
    self.num = 0
    self.sum = 0
    self.sum2 = 0

  def reset(self):
    self.num = 0
    self.sum = self.sum2 = 0.0

  def add(self, value):
    self.num += 1
    self.sum += value
    self.sum2 += value*value

  def mean(self):
    return self.sum / self.num

  def variance(self):
    v = (self.sum2 - self.sum*self.sum/self.num)/(self.num-1)

    # When all added samples are (almost) the same,
    # result `v` is sometimes negative due to a floating point error.
    # This `max` function call avoids this undesirable behavior.
    return max(v, 0.0) 

  def error(self):
    return math.sqrt(self.variance()/self.num)

