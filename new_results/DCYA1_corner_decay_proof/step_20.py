import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps = 40
import corner_series as cs
for N in [30, 60]:
    t=time.time()
    coeffs = cs.corner_coeffs(mp.pi/2, N)
    print(f"N={N}: {round(time.time()-t,2)}s, M_{N}={mp.nstr(coeffs[N],12)}")
