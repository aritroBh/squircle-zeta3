
import time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 105
exec(open('av2core.py').read())
print("=== cform self-convergence p=2 and p=4 ===")
for p in [2, 4]:
    prev = None
    for (h, n) in [('0.09', 110), ('0.07', 140), ('0.05', 200), ('0.04', 240)]:
        t0=time.time(); v = B3_cform(p, h, n, h, n); t1=time.time()
        d = fabs(v-prev) if prev is not None else mpf(0)
        print("p=%s h=%s n=%d: |d|=%s (%.0fs)  %s" % (p, h, n, mp.nstr(d,3), t1-t0, mp.nstr(v,45)))
        prev = v
