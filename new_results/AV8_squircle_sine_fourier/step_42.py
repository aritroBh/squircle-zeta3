from mpmath import mp, mpf
mp.dps = 55
import sys, time, json; sys.path.insert(0,'.')
from squircle import bn_der
ckpt='scratchpad/bnder60.json'
res = json.load(open(ckpt))
t0=time.time()
v = bn_der(11,3)
res["11,3"]=mp.nstr(v,50); json.dump(res,open(ckpt,'w'))
print("11,3:", mp.nstr(v,35), f"{time.time()-t0:.0f}s", flush=True)
