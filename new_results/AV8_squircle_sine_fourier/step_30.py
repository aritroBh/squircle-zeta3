from mpmath import mp, mpf
mp.dps = 60
import sys, time, json, os; sys.path.insert(0,'.')
from squircle import bn_der
ckpt='scratchpad/bnder60.json'
res = json.load(open(ckpt)) if os.path.exists(ckpt) else {}
t0=time.time()
for n in [1,3,5,7,9,11,13]:
    for k in [1,2]:
        key=f"{n},{k}"
        if key in res: continue
        v = bn_der(n,k)
        res[key]=mp.nstr(v,55)
        json.dump(res, open(ckpt,'w'))
        print(key, mp.nstr(v,32), f"{time.time()-t0:.0f}s", flush=True)
print("CHUNK DONE")
