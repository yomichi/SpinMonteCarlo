execfile('loop.py')

L = 8
J = 1.0
Ts = np.arange(0.1, 1.01, 0.1)

usus_output = open('usus.dat', 'w')
ssus_output = open('ssus.dat', 'w')
ene_output = open('energy.dat', 'w')

for T in Ts:
    looper = Looper(L,T,J)
    usus, ssus, ene = looper.run()
    usus_output.write(str(T) + " " + str(usus.mean()) + " " + str(usus.error()) + '\n')
    ssus_output.write(str(T) + " " + str(ssus.mean()) + " " + str(ssus.error()) + '\n')
    ene_output.write(str(T) + " " + str(ene.mean()) + " " + str(ene.error()) + '\n')
