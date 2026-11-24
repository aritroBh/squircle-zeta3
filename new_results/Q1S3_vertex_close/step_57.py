import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp
mp.mp.dps=80
def fembc(x): return femb(x,1)
an = fembc(sa)
ge=GEval(80, M=48, gln=5)
arr, tot = ge.rest_all(an, an)
print("G(a,a;1) total:", mp.nstr(tot,18), " exact:", mp.nstr(mp.log(1-1/an)**2/2,18))
# check at a mid point
def Grec(w, z=1):
    if len(w)==1:
        s=fembc(w[0])
        if feq(w[0],ZERO): return mp.log(z)
        return mp.log(1-z/s)
    s=fembc(w[0])
    return mp.quad(lambda t: Grec(w[1:],t)/(t-s), [0,z])
for j in [0, 5, 100, 700, 1300]:
    x = ge.pts[j][0]
    print("x=",mp.nstr(x,6), " arr:", mp.nstr(arr[j],16), " Grec:", mp.nstr(Grec((sa,sa),x),16), " resid:", mp.nstr(abs(arr[j]-Grec((sa,sa),x)),4))
# now outer integral of G(a,a;t)/t
num = sum(w*arr[j]/x for j,(x,w) in enumerate(ge.pts))
print("outer:", mp.nstr(num,16), " Grec(0,a,a):", mp.nstr(Grec((s0,sa,sa)),16))
# try direct quad of G(a,a;t)/t with Grec inner
direct = mp.quad(lambda t: Grec((sa,sa),t)/t, [0, mp.mpf('0.4'), 1])
print("direct quad:", mp.nstr(direct,16))
