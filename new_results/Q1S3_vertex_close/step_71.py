import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp, time
mp.mp.dps=90
r2=mp.sqrt(2); a=1+r2; b=r2-1
ge=GEval(90, M=48, gln=4, subn=4)
_gw={}
def acc_value(acc):
    tot=mp.mpc(0)
    ws=set(tuple(femb(s,1) for s in w) for (w,lvm) in acc.W if len(w)>=2)
    new=[w for w in ws if w not in _gw]
    if new:
        _gw.update(ge.word_batch([w for w in new if len(w)==3]))
        for w in new:
            if len(w)==2: _gw[w]=ge.rest_total(*w)
    for (w,lvm),c in acc.W.items():
        wn=tuple(femb(s,1) for s in w)
        v=ge.word(wn) if len(w)==1 else _gw[wn]
        for f in lvm: v*=mp.log(femb(f[1],1))
        tot+=femb(c,1)*v
    for lvm,c in acc.P.items():
        v=mp.mpc(1)
        for f in lvm: v*=mp.log(femb(f[1],1))
        tot+=femb(c,1)*v
    tot+=femb(acc.R,1)
    return tot
def factor_val(f,y):
    if f[0]=='C': return mp.log(2)
    if feq(f[1],ZERO): return mp.log(y)
    return mp.log(1-y/femb(f[1],1))
def hg(y,E): return sum(femb(c,1)*factor_val(f,y) for f,c in E.items())
P12 = le_prod(g1,g2)
hr2 = fscale(r2, Fr(1,2))
t0=time.time()
accA=Acc(); int_simple_kernel(accA, hr2, sa, P12)
dA = mp.quad(lambda y: (r2/2)/(y-a)*hg(y,g1)*hg(y,g2), [0,1])
print("kA: engine", mp.nstr(acc_value(accA),18), " direct", mp.nstr(dA,18))
accB=Acc(); int_simple_kernel(accB, fneg(hr2), smb, P12)
dB = mp.quad(lambda y: -(r2/2)/(y+b)*hg(y,g1)*hg(y,g2), [0,1])
print("kB: engine", mp.nstr(acc_value(accB),18), " direct", mp.nstr(dB,18))
print("time:",round(time.time()-t0,1))
