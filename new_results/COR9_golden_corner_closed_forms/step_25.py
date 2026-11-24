import time
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
mp.mp.dps = 1150
from cor9_lib import L7
t0 = time.time()
rt2 = mp.sqrt(2); b = rt2 - 1; c = b*b   # 3 - 2 sqrt2
print("c0 =", mp.nstr(c, 30))
print("alpha = 2 arccos c0 =", mp.nstr(2*mp.acos(c), 30))
M3, W4, W3, W2 = L7.M3_weights(c)
dt = time.time()-t0
print(f"HPL series done, {dt:.0f}s")
print("M3 =", mp.nstr(M3, 40))
print("W4 =", mp.nstr(W4, 30))
print("W3 =", mp.nstr(W3, 30))
print("W2 =", mp.nstr(W2, 30))
for nm, v in [("M3", M3), ("W4", W4), ("W3", W3), ("W2", W2)]:
    open(WD+f"b2sq_{nm}_hpl_1150.txt", "w").write(mp.nstr(v, 1145))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_25: b2sq M3/W4/W3/W2 at dps 1150 via HPL ({dt:.0f}s).\n")
print("saved")
