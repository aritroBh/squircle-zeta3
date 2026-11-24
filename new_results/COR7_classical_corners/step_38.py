import sys, time
sys.path.insert(0, '<project>/research/10_frontier/G1_zeta3_any_angle')
import mpmath as mp
import corner_series as cs

c = mp.mpf(1)/3
for dps in (100, 200):
    mp.mp.dps = dps
    alpha = 2*mp.acos(c)
    t0 = time.time()
    M3q = cs.corner_coeffs(alpha, 3)[3]
    print(f"dps {dps}: {time.time()-t0:.1f}s, M3 = {mp.nstr(M3q, 30)}")
