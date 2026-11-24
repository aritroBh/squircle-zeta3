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
ge=GEval(80, M=48, gln=4, subn=3)
print("gl nodes:", len(ge.gl), " glsub nodes:", len(ge.glsub))
bn=fembc(smb); an=fembc(sa)
F,tot=ge.rest(bn,an)
mx=0
for x in [mp.mpf('0.87'), mp.mpf('0.31'), mp.mpf('0.045'), mp.mpf('0.995')]:
    r=abs(F(x)-Grec((smb,sa),x))
    mx=max(mx,r)
print("rest(-b,a) max err:", mp.nstr(mx,4))
# outer integral with exact rest injected:
Fex=lambda x: Grec((smb,sa),x)
num = mp.mpc(0)
for m in range(ge.M):
    num+=sum(ww*Fex(x)/x for x,ww in ge.pan[m])
print("outer with exact rest:", mp.nstr(abs(num-Grec((s0,smb,sa))),4))
# so where is the loss? check F on the FIRST panel where 1/x is large
for x in [mp.mpf('0.001'), mp.mpf('0.01')]:
    print("x=",mp.nstr(x,5),"F err:", mp.nstr(abs(F(x)-Grec((smb,sa),x)),4))
