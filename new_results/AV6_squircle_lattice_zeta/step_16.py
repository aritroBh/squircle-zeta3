src = open("squircle.py").read()

# 1) Jint returns real part
src = src.replace(
"""    if x <= 1:
        return Aint(p,a) - x*hyper([a, 1/p],[1+1/p], -x**p)
    else:
        return x**(1-a*p)/(a*p-1)*hyper([a, a-1/p],[a-1/p+1], -x**(-p))""",
"""    if x <= 1:
        return (Aint(p,a) - x*hyper([a, 1/p],[1+1/p], -x**p)).real
    else:
        return (x**(1-a*p)/(a*p-1)*hyper([a, a-1/p],[a-1/p+1], -x**(-p))).real""")

# 2) adaptive delta: only use kappa formula for p>=1.5; p<=1 uses delta=4 with more series terms
src = src.replace(
"""    if delta is None:
        if p <= 1:
            delta = mpf('1.0')
        else:
            d1 = (2*pi/(1.3*D))**p
            kap = 2*pi*(1-1/p)*(2*pi/p)**(1/(p-1))*sin(pi/(2*(p-1)))
            d2 = (kap/(1.3*D))**(p-1)
            delta = min(mpf('0.02'), d1, d2)""",
"""    if delta is None:
        if p <= 1:
            delta = mpf('4.0')
        elif p < 1.5:
            delta = min(mpf('0.02'), (2*pi/(1.3*D))**p)
        else:
            d1 = (2*pi/(1.3*D))**p
            kap = 2*pi*(1-1/p)*(2*pi/p)**(1/(p-1))*sin(pi/(2*(p-1)))
            d2 = (abs(kap)/(1.3*D))**(p-1)
            delta = min(mpf('0.02'), d1, d2)""")

# 3) Jmax default larger for p<=1
src = src.replace("def ZB(p, s, delta=None, Jmax=200, info=False):",
"""def ZB(p, s, delta=None, Jmax=None, info=False):
    if Jmax is None:
        Jmax = 160 if mpf(p) <= 1 else 200""")

open("squircle.py","w").write(src)
test = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time
setdps(60)
for (p_,s_) in [(0.5,4),(1.05,4),(25,4)]:
    t0=time.time(); ZA_=ZA(p_,s_); t1=time.time()
    ZB_,iB=ZB(p_,s_,info=True); t2=time.time()
    print(f"p={p_} s={s_}: ZA {t1-t0:.1f}s ZB {t2-t1:.1f}s agree={nstr(fabs(ZA_-ZB_)/fabs(ZA_),4)}")
    print("   ZA:", nstr(ZA_,30)); print("   ZB:", nstr(ZB_,30))
'''
exec(test)
