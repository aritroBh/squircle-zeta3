import sys, time
sys.path.insert(0, "<project>/research/10_frontier/G1_zeta3_any_angle")
import mpmath as mp
import corner_series as csmod
mp.mp.dps = 1100
t0 = time.time()
Ms1100 = csmod.corner_coeffs(2*mp.pi/3, 3)
t1 = time.time()
print(f"dps 1100: {t1-t0:.1f}s")
M3 = Ms1100[3]
print("M3(2pi/3) first 60 digits:", mp.nstr(M3, 60))
open("M3_2pi3_1100.txt", "w").write(mp.nstr(M3, 1095))
# cross-check against the 700-digit value
v700 = mp.mpf(open("M3_2pi3_700.txt").read())
print("|diff 1100 vs 700| =", mp.nstr(abs(M3 - v700), 3))
