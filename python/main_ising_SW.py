execfile('ising_SW.py')


out_am = open('abs_mag.dat', 'w')
out_m2 = open('mag2.dat', 'w')
out_m4 = open('mag4.dat', 'w')
out_sus = open('suscep.dat', 'w')
out_ene = open('ene.dat', 'w')
out_bin = open('binder.dat', 'w')

Ls = [8, 16, 24]
Ts = [1.9, 2.1, 2.3, 2.5]
for L in Ls:
  for T in Ts:
    sw = SW(L, T)
    m2, m4, sus, ene = sw.run()

    m2m = m2.mean()
    m4m = m4.mean()
    binder = m4m/(3*m2m*m2m)

    out_am.write(str(T) + " " + str(am.mean()) + " " + str(am.error()) + "\n")
    out_m2.write(str(T) + " " + str(m2.mean()) + " " + str(m2.error()) + "\n")
    out_m4.write( str(T) + " " + str(m4.mean()) + " " + str(m4.error()) + "\n")
    out_sus.write( str(T) + " " + str(sus.mean()) + " " + str(sus.error()) + "\n")
    out_ene.write( str(T) + " " + str(ene.mean()) + " " + str(ene.error()) + "\n")
    out_bin.write( str(T) + " " + str(binder) + "\n")
    print 'L = ', L, ', T = ', T, ' done.'
  out_am.write('\n\n')
  out_m2.write('\n\n')
  out_m4.write('\n\n')
  out_sus.write('\n\n')
  out_ene.write('\n\n')
  out_bin.write('\n\n')

