code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time, json
setdps(60)
def rel(a,b):
    d = fabs(a-b); return d/fabs(b) if d else mpf('1e-999')
results = []
ps = [0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5, 2, 2.5, 3]
t0=time.time()
for p_ in ps:
    for s_ in [3,4,5,6]:
        ta=time.time()
        ZA_ = ZA(p_, s_)
        ZC_ = ZC(p_, s_, K=500 if p_<3 else 400, R=14)
        r = rel(ZA_, ZC_)
        results.append((p_, s_, str(ZA_), str(r)))
        print(f"p={p_:<5} s={s_}: A-C={nstr(r,4):>9}  ({time.time()-ta:.0f}s)")
json.dump(results, open("sweep1.json","w"))
print("total", time.time()-t0)
'''
exec(code)
json_dump_note = "saved sweep1.json"
with open("progress.log","a") as f: f.write("step_32: sweep part 1 (p=0.5..3) ZA vs ZC\n")
