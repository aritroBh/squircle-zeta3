import time
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
mp.mp.dps = 1150
from cor9_lib import L7
r5 = mp.sqrt(5); phi = (1+r5)/2
for tag, c in [("gol1", 1/phi), ("gol3", r5-2)]:
    t0 = time.time()
    M3, W4, W3, W2 = L7.M3_weights(c)
    dt = time.time()-t0
    print(f"=== {tag}: c0 = {mp.nstr(c, 25)}  alpha = {mp.nstr(2*mp.acos(c), 25)}  ({dt:.0f}s)")
    print("  M3 =", mp.nstr(M3, 40))
    print("  W4 =", mp.nstr(W4, 25))
    print("  W3 =", mp.nstr(W3, 25))
    print("  W2 =", mp.nstr(W2, 25))
    for nm, v in [("M3", M3), ("W4", W4), ("W3", W3), ("W2", W2)]:
        open(WD+f"{tag}_{nm}_hpl_1150.txt", "w").write(mp.nstr(v, 1145))
    open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_33: {tag} M3/W4/W3/W2 dps 1150 via HPL ({dt:.0f}s).\n")
print("saved")
