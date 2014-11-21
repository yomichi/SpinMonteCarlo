from loop import *

L = 8
J = 1.0
Ts = numpy.arange(0.1, 1.01, 0.1)

usus_output = open('ssus.dat', 'w')
ssus_output = open('ssus.dat', 'w')
ene_output = open('energy.dat', 'w')

for T in Ts:
    looper = Looper(L,T,J)
    usus, ssus, ene = looper.run()
    usus_output.write(str(T) + " " + str(usus.mean()) + " " + str(usus.error()))
    ssus_output.write(str(T) + " " + str(ssus.mean()) + " " + str(ssus.error()))
    ene_output.write(str(T) + " " + str(ene.mean()) + " " + str(ene.error()))
