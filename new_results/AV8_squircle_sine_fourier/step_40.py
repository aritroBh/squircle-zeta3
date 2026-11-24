from mpmath import mp, mpf
mp.dps = 60
import sys, time, json; sys.path.insert(0,'.')
from squircle import bn_der
ckpt='scratchpad/bnder60.json'
res = json.load(open(ckpt))
t0=time.time()
for key,n,k in [("9,3",9,3),("11,3",11,3)]:
    v = bn_der(n,k)
    res[key]=mp.nstr(v,55); json.dump(res,open(ckpt,'w'))
    print(key, mp.nstr(v,35), f"{time.time()-t0:.0f}s", flush=True)
