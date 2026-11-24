
import time
from mpmath import mp, mpf, fabs
mp.dps = 110
exec(open('av2core.py').read())
# validate F_cpQ against F_direct values from previous run
import mpmath
def F_direct(c, p):
    ep = 1/p
    def f(r):
        lQ = ep*log1p(-r**p); Q = exp(lQ)
        one_Q = -mpmath.expm1(lQ)
        val = one_Q + (r*r*c)*Q
        return r*(-log(val))/(1-r*r*c)
    return mpmath.quad(f, [0, 1])
for c in [mpf('0.3'), mpf('0.9')]:
    ref = F_direct(c, 1024)
    for (h,n) in [('0.09',100),('0.06',140),('0.045',180)]:
        v = F_cpQ(c, 1024, h, n)
        print("c=%s p=1024 h=%s: F_cpQ err=%s" % (c, h, mp.nstr(fabs(v-ref),4)))
print("=== B3_cformQ convergence p=2, 64, 1024 ===")
for p in [2, 64, 1024]:
    prev = None
    for (h, n) in [('0.09', 110), ('0.06', 150), ('0.045', 190)]:
        t0=time.time(); v = B3_cformQ(p, h, n, h, n); t1=time.time()
        d = fabs(v-prev) if prev is not None else mpf(0)
        print("p=%5s h=%s: |d|=%s (%.0fs)  %s" % (p, h, mp.nstr(d,3), t1-t0, mp.nstr(v,40)))
        prev = v
