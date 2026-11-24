from mpmath import mp, mpf
mp.dps = 60
import sys, time, json, os; sys.path.insert(0,'.')
from squircle import bn_der
ckpt='scratchpad/bnder60.json'
res = json.load(open(ckpt))
t0=time.time()
v = bn_der(1,3)
print("1,3:", mp.nstr(v,40), f"{time.time()-t0:.0f}s")
res["1,3"]=mp.nstr(v,55); json.dump(res,open(ckpt,'w'))
t0=time.time()
v = bn_der(9,2)
print("9,2:", mp.nstr(v,40), f"{time.time()-t0:.0f}s")
res["9,2"]=mp.nstr(v,55); json.dump(res,open(ckpt,'w'))
