from mpmath import mp, mpf, pi, quad
mp.dps = 60
import sys, time; sys.path.insert(0,'.')
from squircle import H1u, H2u, H3u, L1, L2, L3
for phi in [0.5, 1.0, 1.4, 1.55, 1.5707]:
    t0=time.time(); v=L2(phi); dt=time.time()-t0
    print(f"L2({phi}): {dt:.2f}s  {mp.nstr(v,20)}")
