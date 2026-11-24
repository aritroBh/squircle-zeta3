
import time
from mpmath import mp, mpf, fabs
mp.dps = 110
exec(open('av2core.py').read())
# reference F values from F_direct earlier
refs = {('0.3', 1024): mpf('1.25925552088194598763469778447'),
        ('0.9', 1024): mpf('0.856840633762388243601622360654')}
for (cs, p) in [('0.3', 1024), ('0.9', 1024)]:
    ref = refs[(cs, p)]
    for (h, n) in [('0.09', 100), ('0.06', 140), ('0.045', 180)]:
        v = F_cpU(mpf(cs), p, h, n)
        print("c=%s p=%d h=%s: F_cpU err vs ref=%s" % (cs, p, h, mp.nstr(fabs(v-ref),4)))
print("=== B3_cformU p=2, 64, 1024 ===")
for p in [2, 64, 1024]:
    prev = None
    for (h, n1, n2) in [('0.09', 110, 90), ('0.06', 150, 110), ('0.045', 190, 140)]:
        t0=time.time(); v = B3_cformU(p, h, n1, h, n2); t1=time.time()
        d = fabs(v-prev) if prev is not None else mpf(0)
        print("p=%5s h=%s: |d|=%s (%.0fs)  %s" % (p, h, mp.nstr(d,3), t1-t0, mp.nstr(v,42)))
        prev = v
