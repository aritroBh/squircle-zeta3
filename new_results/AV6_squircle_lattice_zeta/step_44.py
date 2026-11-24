code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time, json
setdps(60)
def rel(a,b):
    d = fabs(a-b); return d/fabs(b) if d else mpf('1e-999')
results = []
ps = [4, 5, 7, 10, 15, 25, 50, 100]
t0=time.time()
for p_ in ps:
    for s_ in [3,4,6]:
        ta=time.time()
        ZA_ = ZA(p_, s_)
        ZC_ = ZC(p_, s_, K=max(400, int(6*p_)), R=16)
        r = rel(ZA_, ZC_)
        results.append((p_, s_, str(ZA_), str(r)))
        print(f"p={p_:<5} s={s_}: A-C={nstr(r,4):>9}  ({time.time()-ta:.0f}s)")
json.dump(results, open("sweep2.json","w"))
print("total", time.time()-t0)
'''
exec(code)
with open("progress.log","a") as f: f.write("step_33: sweep part 2 (p=4..100) ZA vs ZC\n")
