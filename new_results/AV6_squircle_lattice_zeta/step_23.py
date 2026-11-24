code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, quad, exp, log, power, fsum
setdps(70)
p = mpf(50); s = mpf(4); a = s/p
mpows = [power(m, p) for m in range(1, 4000)]
def theta(t):
    f = mpf(0); cutoff = mpf(10)**(-85); m = 0
    while True:
        m += 1
        e = exp(-t*mpows[m-1])
        f += e
        if e < cutoff: break
    return 1+2*f
def F(t):
    th = theta(t); return th*th-1
delta = mpf('2.46e-150'); ld = float(log(delta))
import math
def chunked(chunk):
    pts = [ld]
    while pts[-1] < -1e-9:
        pts.append(min(0.0, pts[-1]+chunk))
    pts[-1] = 0.0
    tot = mpf(0)
    for i in range(len(pts)-1):
        lo, hi = pts[i], pts[i+1]
        tot += quad(lambda x: exp(a*x)*F(exp(x)), [lo, hi])
    return tot
for ch in [50, 25, 12]:
    t0=__import__('time').time()
    v = chunked(ch)
    print(f"chunk={ch}: {nstr(v, 45)}  ({__import__('time').time()-t0:.0f}s)")
'''
exec(code)
