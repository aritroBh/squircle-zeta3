
import time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 95
exec(open('av2core.py').read())
# reference: B3(1) from mpmath quad run at dps=70
ref1 = mpf('0.246617302384923740565168068839592055997186102')
print("=== convergence of custom TS, p=1, h=0.25 ===")
for n in [40, 60, 80, 100]:
    t0=time.time(); v = B3_tsA(1, 0.25, n); t1=time.time()
    print("n=%3d  err vs ref = %s   (%.1fs)" % (n, mp.nstr(fabs(v-ref1),4), t1-t0))
print("=== p=2, h=0.25 ===")
for n in [60, 80, 100]:
    v = B3_tsA(2, 0.25, n)
    print("n=%3d  A=%s" % (n, mp.nstr(v, 50)))
print("=== cross-method at p=2, n=100 ===")
print("B3_tsB:", mp.nstr(B3_tsB(2, 0.25, 100), 50))
