import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp
mp.mp.dps=80
r2=mp.sqrt(2); a=1+r2; b=r2-1
ge=GEval(80, M=48, gln=5)
def acc_value(acc):
    tot=mp.mpc(0)
    for (w,lvm),c in acc.W.items():
        wn=tuple(femb(s,1) for s in w)
        v=ge.word(wn)
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
P12 = le_prod(g1,g2)
totE=mp.mpc(0); totD=mp.mpc(0)
for (c,fA,fB) in P12:
    acc=Acc(); int_simple_kernel(acc, ONE, s0, [(c,fA,fB)])
    vE=acc_value(acc)
    cc=femb(c,1)
    vD=mp.quad(lambda y: cc*factor_val(fA,y)*factor_val(fB,y)/y, [0, mp.mpf('0.5'), 1])
    totE+=vE; totD+=vD
    flag = "" if abs(vE-vD)<1e-60 else "  <-- MISMATCH"
    if flag: print("term", fA, fB, "engine", mp.nstr(vE,14), "direct", mp.nstr(vD,14), flag)
print("totals: engine", mp.nstr(totE,15), " direct", mp.nstr(totD,15))
