import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp, time
mp.mp.dps=80
def fembc(x): return femb(x,1)
def Grec(w, z=1, tol=None):
    if len(w)==1:
        s=fembc(w[0])
        if feq(w[0],ZERO): return mp.log(z)
        return mp.log(1-z/s)
    s=fembc(w[0])
    return mp.quad(lambda t: Grec(w[1:],t)/(t-s), [0,z], tol=tol) if tol else mp.quad(lambda t: Grec(w[1:],t)/(t-s), [0,z])
ge=GEval(80, M=48, gln=4, subn=4)
for w in [(s0,smb,sa),(s0,sib,smib),(s0,sa,sa)]:
    wn=tuple(fembc(s) for s in w)
    t0=time.time(); v1=ge.word(wn); t1=time.time()
    v0=Grec(w); t2=time.time()
    print(w[0]==ZERO and '0' or 'x', "GEval:", mp.nstr(v1,18))
    print("   Grec:", mp.nstr(v0,18), " resid:", mp.nstr(abs(v1-v0),4), " times:", round(t1-t0,1), round(t2-t1,1))
