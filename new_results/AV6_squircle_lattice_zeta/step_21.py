code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time
setdps(70)
p_, s_ = 50, 4
t0=time.time()
ZA1 = ZA(p_, s_, N=40, M=300, K=250, R=16, Icorr=16)
ZA2 = ZA(p_, s_, N=80, M=800, K=400, R=20, Icorr=24)
print("ZA base :", nstr(ZA1, 40), f"{time.time()-t0:.0f}s")
print("ZA big  :", nstr(ZA2, 40))
ZB1, iB1 = ZB(p_, s_, info=True)
print("ZB base :", nstr(ZB1, 40), " delta:", nstr(iB1['delta'],3))
ZB2 = ZB(p_, s_, delta=iB1['delta']*30)
print("ZB d30  :", nstr(ZB2, 40))
ZB3 = ZB(p_, s_, delta=iB1['delta']/30)
print("ZB d/30 :", nstr(ZB3, 40))
'''
exec(code)
