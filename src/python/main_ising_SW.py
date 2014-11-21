from  ising_SW import *

L = 16
T = 1.0

sw = SW(L, T)

m2, ene = sw.run()

output = str(T)
output += " " + str(m2.mean()) + " " + str(m2.error())
output += " " + str(ene.mean()) + " " + str(ene.error())

print(output)
