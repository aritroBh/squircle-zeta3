import sys, time
sys.path.insert(0, "<project>/research/10_frontier/G1_zeta3_any_angle")
import mpmath as mp
import corner_series as csmod
mp.mp.dps = 700
t0 = time.time()
Ms700 = csmod.corner_coeffs(2*mp.pi/3, 3)
t1 = time.time()
print(f"dps 700: {t1-t0:.1f}s")
print("M3(2pi/3) 700d:", mp.nstr(Ms700[3], 30))
open("M3_2pi3_700.txt", "w").write(mp.nstr(Ms700[3], 695))
