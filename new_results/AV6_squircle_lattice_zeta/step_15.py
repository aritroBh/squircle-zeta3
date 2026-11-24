code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time
setdps(60)
for (p_,s_) in [(0.5,4),(1.0,4),(1.05,4),(25,4),(50,4)]:
    t0=time.time(); ZA_=ZA(p_,s_); t1=time.time()
    print(f"p={p_} s={s_}: ZA {t1-t0:.1f}s -> {nstr(ZA_,20)}")
    t2=time.time(); ZB_,iB=ZB(p_,s_,info=True); t3=time.time()
    print(f"        ZB {t3-t2:.1f}s -> {nstr(ZB_,20)}  delta={nstr(iB['delta'],3)}")
'''
exec(code)
