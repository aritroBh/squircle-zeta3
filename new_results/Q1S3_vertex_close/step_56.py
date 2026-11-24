import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp
mp.mp.dps=80
def fembc(x): return femb(x,1)
def Grec(w, z=1):
    if len(w)==1:
        s=fembc(w[0])
        if feq(w[0],ZERO): return mp.log(z)
        return mp.log(1-z/s)
    s=fembc(w[0])
    return mp.quad(lambda t: Grec(w[1:],t)/(t-s), [0,z])
ge=GEval(80, M=48, gln=5)
for w in [(s0,sa,smia),(s0,smia,sa),(s0,sa,sa),(s0,smb,sa),(s0,sa,smb)]:
    wn=tuple(fembc(s) for s in w)
    v1=ge.word(wn); v0=Grec(w)
    nm=['0' if feq(s,ZERO) else ('a' if feq(s,sa) else ('-b' if feq(s,smb) else ('-ia' if feq(s,smia) else '?'))) for s in w]
    print(nm, "GEval", mp.nstr(v1,16))
    print("     Grec ", mp.nstr(v0,16), " resid:", mp.nstr(abs(v1-v0),4))
