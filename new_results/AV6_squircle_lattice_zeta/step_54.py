code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time, json
setdps(65)
G = {}
G.update(json.load(open("epsgrid34.json")))
G.update(json.load(open("epsgrid56.json")))
t0=time.time()
for p_ in [160, 192, 256]:
    K = 6*p_
    for s_ in [3,4,5,6]:
        ta=time.time()
        ZC_ = ZC(p_, s_, K=K, R=18)
        G[f"{p_},{s_}"] = str(ZC_)
        print(f"p={p_} s={s_}: {nstr(ZC_, 25)}  ({time.time()-ta:.0f}s)")
json.dump(G, open("epsgrid_all.json","w"))
print("total", time.time()-t0)
'''
exec(code)
with open("progress.log","a") as f: f.write("step_39: eps grid extended to p=256\n")
