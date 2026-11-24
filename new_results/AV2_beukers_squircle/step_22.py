
import time
from mpmath import mp, mpf, fabs
mp.dps = 120
exec(open('av2core.py').read())
print("=== inner rate: F(0.3, 1024) ===")
# best reference: h2=0.002
FPref = FPrecomp(1024, '0.002', 2300)
ref = FPref.F(mpf('0.3'))
print("ref(h2=0.002):", mp.nstr(ref, 55))
for (h2, n2) in [('0.008', 570), ('0.006', 760), ('0.005', 910), ('0.004', 1150), ('0.003', 1530)]:
    FP = FPrecomp(1024, h2, n2)
    v = FP.F(mpf('0.3'))
    print("h2=%s: err=%s" % (h2, mp.nstr(fabs(v-ref),4)))
print("=== outer rate at p=1024 (inner h2=0.003) ===")
FP = FPrecomp(1024, '0.003', 1530)
prev = None
for (h1, n1) in [('0.04', 115), ('0.03', 152), ('0.024', 190), ('0.018', 250)]:
    t0=time.time(); v = B3_fast(p=1024, h1=h1, n1=n1, FP=FP); t1=time.time()
    d = fabs(v-prev) if prev is not None else mpf(0)
    print("h1=%s: |d|=%s (%.0fs) %s" % (h1, mp.nstr(d,3), t1-t0, mp.nstr(v,45)))
    prev = v
