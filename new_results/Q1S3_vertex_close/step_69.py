import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp, time
mp.mp.dps=90
r2=mp.sqrt(2); a=1+r2; b=r2-1
ge=GEval(90, M=48, gln=4, subn=4)
def acc_value(acc):
    tot=mp.mpc(0)
    ws=[tuple(femb(s,1) for s in w) for (w,lvm) in acc.W if len(w)>=2]
    vals=ge.word_batch([w for w in ws if len(w)==3])
    v2={w:ge.rest_total(*w) for w in ws if len(w)==2}
    for (w,lvm),c in acc.W.items():
        wn=tuple(femb(s,1) for s in w)
        v=ge.word(wn) if len(w)==1 else (v2[wn] if len(w)==2 else vals[wn])
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
t0=time.time()
bad=[]
totE=mp.mpc(0); totD=mp.mpc(0)
for (c,fA,fB) in P12:
    acc=Acc(); int_simple_kernel(acc, ONE, s0, [(c,fA,fB)])
    vE=acc_value(acc)
    cc=femb(c,1)
    vD=mp.quad(lambda y: cc*factor_val(fA,y)*factor_val(fB,y)/y, [0, mp.mpf('0.5'), 1])
    totE+=vE; totD+=vD
    if abs(vE-vD)>1e-70:
        bad.append((fA,fB,vE,vD))
for fA,fB,vE,vD in bad:
    print("MISMATCH", fA[1] if fA[0]=='G' else 'C', fB[1] if fB[0]=='G' else 'C')
    print("   engine", mp.nstr(vE,15))
    print("   direct", mp.nstr(vD,15))
print("bad:",len(bad)," totals:", mp.nstr(totE,14), mp.nstr(totD,14), " time:", round(time.time()-t0,1))
