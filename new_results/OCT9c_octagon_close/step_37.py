import pickle, time
import mpmath as mp
mp.mp.dps = 350
t0=time.time()
v = mp.polylog(3, mp.mpf(-4374))
print('Li3(-4374):', time.time()-t0, 's')
t0=time.time()
v2 = mp.re(mp.polylog(3, mp.mpf(4375)))
print('Li3(4375) re:', time.time()-t0, 's')
t0=time.time()
v3 = mp.polylog(3, mp.mpf('-1/4374'))
print('Li3(-1/4374):', time.time()-t0, 's')
# check: for real arg>1, mp.polylog returns complex with imaginary part; real part convention
print(mp.polylog(3, mp.mpf(10)))
