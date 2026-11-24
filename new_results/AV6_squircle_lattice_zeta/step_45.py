code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time, json
setdps(60)
def rel(a,b):
    d = fabs(a-b); return d/fabs(b) if d else mpf('1e-999')
results = []
t0=time.time()
for p_ in [0.99, 0.999, 1.001, 1.01]:
    for s_ in [3,4,6]:
        ta=time.time()
        ZA_ = ZA(p_, s_)
        ZC_ = ZC(p_, s_, K=500, R=16)
        r = rel(ZA_, ZC_)
        results.append((p_, s_, str(ZA_), str(r)))
        print(f"p={p_:<6} s={s_}: A-C={nstr(r,4):>9}  ({time.time()-ta:.0f}s)")
json.dump(results, open("sweep3.json","w"))
print("total", time.time()-t0)
'''
exec(code)
with open("progress.log","a") as f: f.write("step_34: sweep part 3 (p near 1)\n")
