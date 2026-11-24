
import os, time
print("loadavg:", os.getloadavg())
print("cpus:", os.cpu_count())
import mpmath as mp
mp.mp.dps = 30
t = time.time()
s = mp.mpf(0)
for i in range(1, 200000):
    s += mp.mpf(i)**2
print("busywork ok", time.time()-t)
