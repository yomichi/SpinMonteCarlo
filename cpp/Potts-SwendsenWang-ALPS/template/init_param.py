import pyalps
import numpy as np

params = []

MCS = 65536

for q in [2,3,4]:
  bc = np.log(1.0 + np.sqrt(q))
  for L in [16, 24, 32]:
    for b in np.linspace(0.9, 1.1, 21):
      beta = bc * b
      params.append({
        'ALGORITHM' : 'cluster',
        'MODEL' : 'Potts',
        'q' : 2,
        'LATTICE' : 'square lattice',
        'L' : L,
        'T' : 1.0/beta,
        'beta' : beta,
        'J' : 1.0,
        'SWEEPS' : MCS,
        'THERMALIZATION' : MCS >> 3,
        })

pyalps.writeInputFiles('params', params)

