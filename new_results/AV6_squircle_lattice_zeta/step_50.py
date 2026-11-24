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
    for s_ in [3, 4]:
        ta=time.time()
        ZC_ = ZC(p_, s_, K=K, R=18)
        vals[(p_,s_)] = str(ZC_)
        print(f"p={p_} s={s_}: {nstr(ZC_, 30)}  ({time.time()-ta:.0f}s)")
json.dump({f"{a},{b}":v for (a,b),v in vals.items()}, open("epsgrid34.json","w"))
print("total", time.time()-t0)
'''
with open("step_36.py","w") as f: f.write(code)
exec(code)
with open("progress.log","a") as f: f.write("step_36: eps grid s=3,4 computed\n")
