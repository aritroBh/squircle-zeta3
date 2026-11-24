from mpmath import mp, mpf
mp.dps = 60
import sys, time, json, os; sys.path.insert(0,'.')
from squircle import bn_der
ckpt='scratchpad/bnder60.json'
res = json.load(open(ckpt)) if os.path.exists(ckpt) else {}
t0=time.time()
for n in [1,3,5,7]:
    for k in [1,2,3]:
        key=f"{n},{k}"
        if key in res: continue
        v = bn_der(n,k)
        res[key]=mp.nstr(v,55)
        json.dump(res, open(ckpt,'w'))
        print(key, mp.nstr(v,32), f"{time.time()-t0:.0f}s", flush=True)
print("CHUNK DONE")
with open('progress.log','a') as f:
    f.write("step_20: recomputed b_n^(k)(2) dps=60 with corrected H3 (n<=7 chunk)\n")
