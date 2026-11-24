code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, power, gamma
import time
setdps(70)
p_, s_ = 50, 4
ZA_ = ZA(p_, s_, N=40, M=300, K=250, R=16)
print("ZA      :", nstr(ZA_, 45))
t0=time.time()
for dexp in [150, 200, 250, 300]:
    dd = mpf(10)**(-dexp)
    ZB_, iB = ZB(p_, s_, delta=dd, info=True)
    print(f"delta=1e-{dexp}: {nstr(ZB_, 45)}  (small={nstr(iB['small'],8)}) {time.time()-t0:.0f}s")
'''
exec(code)
