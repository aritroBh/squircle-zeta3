src = open('q3_eval.py').read()
old = """        for m in range(M):
            h2=mp.mpf(1)/(2*M); mid=(2*m+1)*h2
            for gx,gw in self.gl:
                self.pts.append((mid+h2*gx, gw*h2))"""
new = """        for m in range(M):
            h2=mp.mpf(1)/(2*M); mid=(2*m+1)*h2
            for gx,gw in sorted(self.gl, key=lambda p: p[0]):
                self.pts.append((mid+h2*gx, gw*h2))"""
assert old in src
open('q3_eval.py','w').write(src.replace(old,new))

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
for w in [(s0,sa,sa),(s0,sa,smia),(s0,smia,sa),(s0,smb,sa),(s0,sa,smb),(s0,sib,smib)]:
    wn=tuple(fembc(s) for s in w)
    v1=ge.word(wn); v0=Grec(w)
    print(mp.nstr(abs(v1-v0),4), " ", end="")
print()
# full calibration set again
r2=mp.sqrt(2)
cases = []
pairs = [(1j*(r2-1),'ib'),(-1j*(r2-1),'-ib'),(1+r2,'a'),(-(r2-1),'-b'),(1j*(1+r2),'ia'),(-1j*(1+r2),'-ia')]
for v,name in pairs:
    cases.append(((mp.mpc(0),mp.mpc(v)), -mp.polylog(2,1/v)))
    cases.append(((mp.mpc(v),mp.mpc(v)), mp.log(1-1/v)**2/2))
mx=0
for w,val in cases:
    mx=max(mx,abs(ge.word(w)-val))
print("calibration max resid:", mp.nstr(mx,4))
