import time
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from cor9_lib import L7
mp.mp.dps = 100
r5 = mp.sqrt(5); phi = (1+r5)/2
for tag, c in [("gol1", 1/phi), ("gol3", r5-2)]:
    t0 = time.time()
    M3q = L7.corner_coeffs_hp(2*mp.acos(c), 3)[3]
    M3h = mp.mpf(open(WD+f"{tag}_M3_hpl_1150.txt").read())
    res = abs(M3q - M3h)
    print(f"{tag}: quad {time.time()-t0:.0f}s  M3q = {mp.nstr(M3q, 30)}  |quad-HPL| = {mp.nstr(res, 2)}")
    open(WD+f"{tag}_M3_quad_100.txt", "w").write(mp.nstr(M3q, 95))
    open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_34: {tag} quadrature dps 100, |quad-HPL| = {mp.nstr(res,2)}.\n")
