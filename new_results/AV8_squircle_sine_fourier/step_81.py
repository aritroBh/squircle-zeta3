from mpmath import mp, mpf
mp.dps = 75
import sys, time, json, os; sys.path.insert(0,'.')
from squircle import bn_der
ckpt='scratchpad/bnder_hi.json'
res = json.load(open(ckpt)) if os.path.exists(ckpt) else {}
t0=time.time()
for key,n,k in [("1,3",1,3),("3,3",3,3)]:
    if key in res and len(res[key])>60: continue
    v = bn_der(n,k)
    res[key]=mp.nstr(v,70); json.dump(res,open(ckpt,'w'))
    print(key, mp.nstr(v,50), f"{time.time()-t0:.0f}s", flush=True)
