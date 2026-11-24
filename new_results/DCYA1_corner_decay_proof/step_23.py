import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
import corner_series as cs

# Worst point u=-0.9, a=pi/2. Compute M at increasing dps to see quadrature convergence.
a=mp.pi/2
for dps in [50, 80, 120]:
    mp.mp.dps=dps
    t=time.time()
    v = cs.M_exact(mp.mpf('-0.9'), a, pts=[0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1])
    print(f"dps={dps}: M(-0.9)={mp.nstr(v,30)} (t={round(time.time()-t,1)}s)")
