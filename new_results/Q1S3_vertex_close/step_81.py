import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp, time
mp.mp.dps=100
r2n=mp.sqrt(2); an=1+r2n; bn=r2n-1
def fembc(x): return femb(x,1)
ge=GEval(100, M=48, gln=4, subn=4)
_gw={}
def acc_value(acc):
    tot=mp.mpc(0)
    ws=set(tuple(fembc(s) for s in w) for (w,lvm) in acc.W if len(w)>=2)
    new=[w for w in ws if w not in _gw]
    if new:
        _gw.update(ge.word_batch([w for w in new if len(w)==3]))
        for w in new:
            if len(w)==2: _gw[w]=ge.rest_total(*w)
    for (w,lvm),c in acc.W.items():
        wn=tuple(fembc(s) for s in w)
        v=ge.word(wn) if len(w)==1 else _gw[wn]
        for f in lvm: v*=mp.log(fembc(f[1]))
        tot+=fembc(c)*v
    for lvm,c in acc.P.items():
        v=mp.mpc(1)
        for f in lvm: v*=mp.log(fembc(f[1]))
        tot+=fembc(c)*v
    tot+=fembc(acc.R)
    return tot
def factor_val(f,y):
    if f[0]=='C': return mp.log(2)
    if feq(f[1],ZERO): return mp.log(y)
    return mp.log(1-y/fembc(f[1]))
def hg(y,E): return sum(fembc(c)*factor_val(f,y) for f,c in E.items())
P12 = le_prod(g1,g2); dP12=prod_deriv(P12)
a_2=fmul(a_,a_); b_2=fmul(b_,b_)
t0=time.time()
G2=Acc(); int_formP(G2, a_2, sa, dP12); boundary_P_at1(G2, fdiv(ONE,b_), P12)
d2 = mp.quad(lambda y: hg(y,g1)*hg(y,g2)/(1-bn*y)**2, [0, mp.mpf('0.5'),1])
print("G2: engine", mp.nstr(acc_value(G2),16), " direct", mp.nstr(d2,16))
G3=Acc(); int_formP(G3, b_2, smb, dP12); boundary_P_at1(G3, fneg(fdiv(ONE,a_)), P12)
d3 = mp.quad(lambda y: hg(y,g1)*hg(y,g2)/(1+an*y)**2, [0, mp.mpf('0.5'),1])
print("G3: engine", mp.nstr(acc_value(G3),16), " direct", mp.nstr(d3,16))
print("time:", round(time.time()-t0,1))
