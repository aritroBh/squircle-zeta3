
import time
from mpmath import mp, mpf, fabs
mp.dps = 120
exec(open('av2core.py').read())
FP = FPrecomp(1024, '0.003', 1530)
print("=== B3_fast2 (peak subtracted) p=1024 ===")
prev = None
for (h1, n1) in [('0.06', 80), ('0.04', 115), ('0.03', 152), ('0.024', 190)]:
    t0=time.time(); v = B3_fast2(1024, h1, n1, FP); t1=time.time()
    d = fabs(v-prev) if prev is not None else mpf(0)
    print("h1=%s: |d|=%s (%.0fs) %s" % (h1, mp.nstr(d,3), t1-t0, mp.nstr(v,45)))
    prev = v
