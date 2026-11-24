
import time
from mpmath import mp, mpf, zeta, fabs
mp.dps = 70
exec(open('av2core.py').read())
print("=== B3 timing + agreement test, dps=70 (stable) ===")
for p in [1, 2, 4, 16, 64]:
    t0=time.time(); A = B3_quadA(p); t1=time.time()
    B = B3_quadB(p); t2=time.time()
    print("p=%3d  A=%s (%.1fs)" % (p, mp.nstr(A,40), t1-t0))
    print("       B=%s (%.1fs)  |A-B|=%s" % (mp.nstr(B,40), t2-t1, mp.nstr(fabs(A-B),4)))
