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
an=fembc(sa)
arr,tot=ge.rest_all(an,an)
mx=0; bad=0
for j in [0,3,17,100,500,1200,2000]:
    x=ge.pts[j][0]
    r=abs(arr[j]-Grec((sa,sa),x))
    print("x=",mp.nstr(x,8),"resid:",mp.nstr(r,4))
# outer integral
num = sum(w*arr[j]/x for j,(x,w) in enumerate(ge.pts))
print("outer (0,a,a):", mp.nstr(num,16), " Grec:", mp.nstr(Grec((s0,sa,sa)),16))
# investigate integrand near t=0: G(a,a;t)/t = ~ -t/a.../t fine
for j in [0,1,2,3]:
    x=ge.pts[j][0]
    print("j",j,"x",mp.nstr(x,10),"arr/x:",mp.nstr(arr[j]/x,12))
