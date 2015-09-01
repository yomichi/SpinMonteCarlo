import pyalps
import numpy as np

params = []

MCS = 65536

for L in [16, 24, 32]:
  for T in np.linspace(1.1, 1.15, 51):
    params.append({
      'ALGORITHM' : 'cluster',
      'MODEL' : 'Potts'
      'q' : 2,
      'LATTICE' : 'square lattice',
      'L' : L,
      'T' : T,
      'J' : 1.0,
      'SWEEPS' : MCS,
      'THERMALIZATION' : MCS >> 3,
      })

pyalps.writeInputFiles('params', params)

