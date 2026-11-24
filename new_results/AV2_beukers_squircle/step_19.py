
import time
from mpmath import mp, mpf, fabs
mp.dps = 110
exec(open('av2core.py').read())
print("=== B3_vform3 (v-outer, split inner) p=64, 1024 ===")
for p in [64, 1024]:
    prev = None
    for (h, n1, n2) in [('0.12', 100, 90), ('0.08', 140, 110), ('0.06', 180, 130)]:
        t0=time.time(); v = B3_vform3(p, h, n1, h, n2); t1=time.time()
        d = fabs(v-prev) if prev is not None else mpf(0)
        print("p=%5s h=%s: |d|=%s (%.0fs)  %s" % (p, h, mp.nstr(d,3), t1-t0, mp.nstr(v,42)))
        prev = v
