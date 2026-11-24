
import time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 105
exec(open('av2core.py').read())
print("=== B3_cform4 cheap levels, p=64,1024: does outer split fix p-dependence? ===")
for p in [64, 1024]:
    prev = None
    for (h, n1, n2) in [('0.09', 110, 80), ('0.06', 150, 100)]:
        t0=time.time(); v = B3_cform4(p, h, n1, h, n2); t1=time.time()
        d = fabs(v-prev) if prev is not None else mpf(0)
        print("p=%5s h=%s: |d|=%s (%.0fs)  %s" % (p, h, mp.nstr(d,3), t1-t0, mp.nstr(v,42)))
        prev = v
print("=== compare cform3 (no outer split) same params p=1024 ===")
v3 = B3_cform3(1024, '0.06', 150, '0.06', 100)
print("cform3:", mp.nstr(v3, 42))
