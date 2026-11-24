from mpmath import mp, mpf
mp.dps = 60
import sys, time, json; sys.path.insert(0,'.')
from squircle import bn_der
res = json.load(open('scratchpad/bnder60.json'))
t0=time.time(); v=bn_der(3,3); print("3,3:", mp.nstr(v,45), f"{time.time()-t0:.0f}s")
print("old :", res["3,3"][:48])
print("diff:", mp.nstr(v-mpf(res["3,3"]),4))
