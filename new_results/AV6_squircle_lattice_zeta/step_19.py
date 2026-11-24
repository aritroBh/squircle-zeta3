src = open("squircle.py").read()
src = src.replace(
"""    for k in range(N+2):
        if k > pint >= 0:
            bk = mpf(0)
        else:
            bk = exp(loggamma(p+1)-loggamma(k+1)-loggamma(p-k+1))
        v[k] = bk * yp / y**k""",
"""    bk = mpf(1)
    for k in range(N+2):
        if k > pint >= 0:
            bk = mpf(0)
        else:
            if k > 0:
                bk *= (p-(k-1))/k
            v[k] = bk * yp / y**k
        if k > pint >= 0:
            v[k] = mpf(0)""")
open("squircle.py","w").write(src)
test = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs
setdps(60)
p = mpf('0.5'); s = mpf(4); a = s/p
ck = taylor_h(5, 250, p, a, 6)
print([type(c).__name__ for c in ck])
# verify against exact binomial series value
g = lambda t: (5**p + (250+t)**p)**(-a)
from mpmath import diff, fac
for k in range(5):
    print(k, nstr(ck[k], 12), nstr(diff(g,0,n=k)/fac(k), 12))
import time
t0=time.time(); print("ZA(0.5,4):", nstr(ZA(0.5,4), 30), f"{time.time()-t0:.1f}s")
'''
exec(test)
