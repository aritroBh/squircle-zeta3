
import time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 95
exec(open('av2core.py').read())
print("=== method B (polar) convergence ===")
for p in [1, 2]:
    print("p=", p)
    prev = None
    for (h, n) in [(0.25, 40), (0.2, 50), (0.16, 60), (0.125, 80), (0.1, 100), (0.08, 120)]:
        t0=time.time(); v = B3_tsB(p, h, n); t1=time.time()
        d = fabs(v-prev) if prev is not None else mpf(0)
        print("  h=%.3f n=%3d  %s  |d|=%s (%.1fs)" % (h, n, mp.nstr(v, 40), mp.nstr(d, 3), t1-t0))
        prev = v
print("=== method A at p=2, h=0.1 ===")
print(mp.nstr(B3_tsA(2, 0.1, 100), 40))
print("=== method B at p=2, h=0.1 ===")
print(mp.nstr(B3_tsB(2, 0.1, 100), 40))
