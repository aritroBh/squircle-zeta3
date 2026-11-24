
import time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 95
exec(open('av2core.py').read())
ref1 = mpf('0.246617302384923740565168068839592055997186102')
print("=== custom TS convergence, p=1, method A, h=0.25 ===")
for n in [40, 60, 80, 100, 120]:
    t0=time.time(); v = B3_tsA(1, 0.25, n); t1=time.time()
    print("n=%3d  err vs ref = %s   (%.1fs)" % (n, mp.nstr(fabs(v-ref1),4), t1-t0))
print("=== h=0.2 ===")
for n in [80, 100, 120]:
    t0=time.time(); v = B3_tsA(1, 0.2, n); t1=time.time()
    print("n=%3d  err vs ref = %s   (%.1fs)" % (n, mp.nstr(fabs(v-ref1),4), t1-t0))
