import sys
sys.path.insert(0, '<project>/research/10_frontier/G1_zeta3_any_angle')
import mpmath as mp
import corner_series as cs

c = mp.cos(2*mp.acos(mp.mpf(1)/3))
for gtest in [mp.mpf(1)/2, mp.mpf('1e-8'), mp.mpf('1e-50'), mp.mpf('1e-150'), mp.mpf('1e-300')]:
    mp.mp.dps = 200
    v = cs.integrand_series(gtest, c, 3)[3]
    mp.mp.dps = 450
    vref = cs.integrand_series(gtest, c, 3)[3]
    mp.mp.dps = 200
    print(f"g = {mp.nstr(gtest,3)}: err(200 vs 450) = {mp.nstr(abs(v-vref),3)}  (|v|={mp.nstr(abs(vref),3)})")
