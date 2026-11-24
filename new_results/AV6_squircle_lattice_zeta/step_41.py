code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time
setdps(60)
def rel(a,b):
    d = fabs(a-b)
    return d/fabs(b) if d else mpf('1e-999')
t0=time.time()
for p_ in [4, 6, 8, 10, 12]:
    ZA_ = ZA(p_, 4)
    ZB_ = ZB(p_, 4)
    print(f"p={p_}: ZA-ZB rel = {nstr(rel(ZA_,ZB_),4)}  ({time.time()-t0:.0f}s)")
'''
exec(code)
with open("progress.log","a") as f: f.write("step_30: ZB degradation map p=4..12\n")
