src = open("squircle.py").read()
src = src.replace(
"""    # large-t quadrature
    def theta(t):
        f = mpf(0)
        m = 1
        while True:
            e = exp(-t*power(m, p))
            f += e
            if e < mpf(10)**(-mp.dps-10): break
            m += 1
        return 1+2*f
    def integrand(t):
        th = theta(t)
        return t**(a-1)*(th*th-1)
    large = quad(integrand, [delta, 1, inf])
    Z = (small + large)/ga
    if info: return Z, dict(small=small, large=large, delta=delta, est1=est1, jstar=jstar)
    return Z""",
"""    # large-t quadrature
    mpows_cache = {}
    def theta(t):
        f = mpf(0)
        m = 1
        cutoff = mpf(10)**(-mp.dps-10)
        while True:
            mpw = mpows_cache.get(m)
            if mpw is None:
                mpw = power(m, p); mpows_cache[m] = mpw
            e = exp(-t*mpw)
            f += e
            if e < cutoff: break
            m += 1
        return 1+2*f
    def F(t):
        th = theta(t)
        return th*th-1
    if delta < 1:
        ld = log(delta)
        I1 = quad(lambda x: exp(a*x)*F(exp(x)), [ld, 0])
        I2 = quad(lambda t: t**(a-1)*F(t), [1, inf])
        large = I1 + I2
    else:
        large = quad(lambda t: t**(a-1)*F(t), [delta, inf])
    Z = (small + large)/ga
    if info: return Z, dict(small=small, large=large, delta=delta, est1=est1, jstar=jstar)
    return Z""")
open("squircle.py","w").write(src)
test = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time
setdps(60)
for (p_,s_) in [(0.5,4),(25,4),(50,4)]:
    t0=time.time(); ZA_=ZA(p_,s_); t1=time.time()
    ZB_,iB=ZB(p_,s_,info=True); t2=time.time()
    print(f"p={p_} s={s_}: ZA {t1-t0:.1f}s ZB {t2-t1:.1f}s agree={nstr(fabs(ZA_-ZB_)/fabs(ZA_),4)} est1={nstr(iB['est1'],4)}")
'''
exec(test)
with open("progress.log","a") as f: f.write("step_18: taylor_h binomial fix; quad ln-transform\n")
