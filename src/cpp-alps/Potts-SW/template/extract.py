import numpy as np
import pyalps
import pyalps.plot as plot

#         Observable name : output_prefix
names = { 'Time' : 'time',
          'Speed' : 'speed',
          'Binder Ratio' : 'binder',
          'Susceptibility' : 'sus',
          'Energy' : 'ene',
          'Magnetization^2' : 'mag2',
          'Magnetization^4' : 'mag4',
          'Specific Heat' : 'spec',
          'Largest Cluster Size' : 'largest',
          'Number of Clusters' : 'cluster',
          }

xnames   = [ 'L', ]
foreachs = [ ['T'], ]
fe_types = [ [np.float], ]

def extract(data, xname, names, foreach, fe_types):
  if np.isscalar(foreach):
    foreach = [foreach]
  if np.isscalar(fe_types):
    fe_types = [fetypes]
  for name in names:
    for obs in pyalps.collectXY(data, xname, name, foreach=foreach):
      vals = [ typ(obs.props[sym]) for sym, typ in zip(foreach, fe_types) ]
      filename = names[name]
      for sym, val in zip(foreach, vals):
        filename += '-{}{}'.format(sym,val)
      filename += '.dat'
      with open(filename, 'w') as f:
        f.write(plot.convertToText([obs]).replace(' +/- ', ' '))


result_files = pyalps.getResultFiles(prefix='params')

data = pyalps.loadMeasurements(result_files, names.keys())

for xname, fe, fet in zip(xnames, foreachs, fe_types):
  extract(data, xname, names, fe, fet)

