
import time
from mpmath import mp, mpf, fabs
mp.dps = 115
exec(open('av2core.py').read())
print("=== inner F accuracy with precomputed grid, p=1024 ===")
t0=time.time(); FP = FPrecomp(1024, '0.004', 1150); t1=time.time()
print("grid build: %.1fs, %d nodes" % (t1-t0, len(FP.grid)))
ref = F_cp(mpf('0.3'), 1024, '0.018', 500)
t0=time.time(); v = FP.F(mpf('0.3')); t1=time.time()
print("F(0.3) err=%s  (%.3fs per F call)" % (mp.nstr(fabs(v-ref),4), t1-t0))
print("=== outer convergence p=1024, 64 ===")
for p in [1024, 64]:
    FP = FPrecomp(p, '0.004', 1150)
    prev = None
    for (h1, n1) in [('0.06', 100), ('0.04', 140), ('0.03', 180)]:
        t0=time.time(); v = B3_fast(p, h1, n1, FP); t1=time.time()
        d = fabs(v-prev) if prev is not None else mpf(0)
        print("p=%5s h1=%s: |d|=%s (%.0fs)  %s" % (p, h1, mp.nstr(d,3), t1-t0, mp.nstr(v,42)))
        prev = v
