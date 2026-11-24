code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, quad, exp, log, power
import time
setdps(70)
p_, s_ = 50, 4
t0=time.time()
ZA1 = ZA(p_, s_, N=40, M=300, K=250, R=16, Icorr=16)
ZA3 = ZA(p_, s_, N=200, M=500, K=300, R=18, Icorr=20)
print("ZA1:", nstr(ZA1, 45), time.time()-t0)
print("ZA3:", nstr(ZA3, 45))
# ZB with much larger delta (still valid for even p since r-series=0; error e^{-kappa delta^{-1/(p-1)}})
kap = 2*pi*(1-1/p_)*(2*pi/p_)**(1/(p_-1))*sin(pi/(2*(p_-1)))
print("kappa:", nstr(kap, 6))
for dexp in [120, 80, 60]:
    dd = mpf(10)**(-dexp)
    print(f"delta=1e-{dexp}: R~e^{nstr(-kap*dd**(-1/mpf(p_-1)),4)}")
    ZBd = ZB(p_, s_, delta=dd)
    print("   ZB:", nstr(ZBd, 45))
'''
exec(code)
