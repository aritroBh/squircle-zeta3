code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time, json
setdps(65)
vals = {}
t0=time.time()
ps = [24, 32, 40, 48, 64, 96, 128]
for p_ in ps:
    K = max(500, 6*p_)
    for s_ in [5, 6]:
        ZC_ = ZC(p_, s_, K=K, R=18)
        vals[(p_,s_)] = str(ZC_)
        print(f"p={p_} s={s_}: {nstr(ZC_, 30)}  ({time.time()-t0:.0f}s)")
json.dump({f"{a},{b}":v for (a,b),v in vals.items()}, open("epsgrid56.json","w"))
print("total", time.time()-t0)
'''
exec(code)
with open("progress.log","a") as f: f.write("step_37: eps grid s=5,6 computed\n")
