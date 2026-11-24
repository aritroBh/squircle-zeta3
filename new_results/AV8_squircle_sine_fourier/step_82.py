from mpmath import mp, mpf
mp.dps = 75
import sys, time, json, os; sys.path.insert(0,'.')
from squircle import bn_der
ckpt='scratchpad/bnder_hi.json'
res = json.load(open(ckpt)) if os.path.exists(ckpt) else {}
t0=time.time()
v = bn_der(5,3)
res["5,3"]=mp.nstr(v,70); json.dump(res,open(ckpt,'w'))
print("5,3:", mp.nstr(v,50), f"{time.time()-t0:.0f}s")
