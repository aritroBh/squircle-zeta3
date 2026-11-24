import time
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from cor9_lib import L7
mp.mp.dps = 100
t0 = time.time()
rt2 = mp.sqrt(2); b = rt2 - 1; c = b*b
alpha = 2*mp.acos(c)
M3q = L7.corner_coeffs_hp(alpha, 3)[3]
dt = time.time()-t0
M3h = mp.mpf(open(WD+"b2sq_M3_hpl_1150.txt").read())
res = abs(M3q - M3h)
print(f"quadrature {dt:.0f}s at dps 100")
print("M3(quad) =", mp.nstr(M3q, 40))
print("M3(HPL)  =", mp.nstr(M3h, 40))
print("|quad - HPL| =", mp.nstr(res, 3))
open(WD+'b2sq_M3_quad_100.txt','w').write(mp.nstr(M3q, 95))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_26: b2sq quadrature dps 100 ({dt:.0f}s), |quad-HPL| = {mp.nstr(res,2)}.\n")
