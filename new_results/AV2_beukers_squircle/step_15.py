
import time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 105
exec(open('av2core.py').read())
print("=== B3_cform3 (split inner) ===")
for p in [2, 4, 64, 1024]:
    prev = None
    for (h, n) in [('0.09', 110), ('0.06', 160), ('0.045', 210)]:
        t0=time.time(); v = B3_cform3(p, h, n, h, n); t1=time.time()
        d = fabs(v-prev) if prev is not None else mpf(0)
        print("p=%5s h=%s: |d|=%s (%.0fs)  %s" % (p, h, mp.nstr(d,3), t1-t0, mp.nstr(v,40)))
        prev = v
