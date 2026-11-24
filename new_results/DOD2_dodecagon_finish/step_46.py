import mpmath as mp, time
mp.mp.dps = 1350
t0=time.time()
r3 = mp.sqrt(3)
X = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
T = mp.quad(lambda tt: mp.log(1-X[2]*tt)*mp.log(1-X[5]*tt)/(1-X[1]*tt), [0,1])
print("one T3a quad at 1350 dps: %ds" % (time.time()-t0))
C = mp.quad(lambda t: mp.log(abs(1-mp.mpf('1.268')*t))*mp.log(abs(1-mp.mpf('0.732')*t))/t, [0,1])
print("one Creal quad at 1350 dps: %ds total" % (time.time()-t0))
L = mp.polylog(3, mp.mpf('0.423'))
print("polylog: %ds total" % (time.time()-t0))
