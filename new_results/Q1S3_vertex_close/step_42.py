src = open('q3_eval.py').read()
src = src.replace("self.gl=gl.calc_nodes(gln, dps)",
                  "self.gl=gl.calc_nodes(gln, int(dps*3.34)+30)")
open('q3_eval.py','w').write(src)

import importlib, q3_engine, q1s_exact, q3_eval
importlib.reload(q1s_exact); importlib.reload(q3_engine); importlib.reload(q3_eval)
from q3_engine import *
from q3_eval import *
import mpmath as mp
def fembc(x): return femb(x,1)
mp.mp.dps=60
def Grec(w, z=1):
    if len(w)==1:
        s=fembc(w[0])
        if feq(w[0],ZERO): return mp.log(z)
        return mp.log(1-z/s)
    s=fembc(w[0])
    return mp.quad(lambda t: Grec(w[1:],t)/(t-s), [0,z])
ge = GEval(60, h='0.012', T='4.6', gln=4)
tests = [(sa,smb,sia),(s0,sa,smb),(sa,s0,sib),(smb,sib,smib),(s0,sib,smia),
         (sa,sa),(s0,sib),(smb,s0),(sia,smb),(sa,s0),(sib,smib),(s0,sa,s0),
         (sib,s0,smb),(smia,sia,s0)]
mx=0
for w in tests:
    wn=tuple(fembc(s) for s in w)
    mx=max(mx,abs(ge.word(wn)-Grec(w)))
print("h=0.012 gln=4 max resid:", mp.nstr(mx,4))
ge2 = GEval(60, h='0.02', T='4.4', gln=4)
mx=0
for w in tests:
    wn=tuple(fembc(s) for s in w)
    mx=max(mx,abs(ge2.word(wn)-Grec(w)))
print("h=0.02 gln=4 max resid:", mp.nstr(mx,4))
