import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps = 40
import corner_series as cs

# Test M_exact at complex u on circle r=0.9
a = mp.pi/2
r = 0.9
for th in [0, mp.pi/2, mp.pi, 3*mp.pi/2]:
    u = r*mp.e**(1j*th)
    t=time.time()
    v = cs.M_exact(u, a)
    print(f"theta={mp.nstr(th,4)}, u={mp.nstr(u,6)}, M={mp.nstr(v,12)}, time={round(time.time()-t,2)}s")
