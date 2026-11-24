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
# check letter assignments of g2: g2 = log(1-iay)+log(1+iay)+log(1-iby)+log(1+iby)-2log(1-by)-2log(1+ay)
# letters: log(1-iay): c=ia -> 1/c=-ib=smib ; log(1+iay): c=-ia -> 1/c=ib=sib
#          log(1-iby): c=ib -> 1/c=-ia=smia ; log(1+iby): c=-ib -> 1/c=ia=sia
print("g2 letters:", sorted(set(k[1] for k in g2 if k[0]=='G'), key=repr))
print("expected: smib,sib,smia,sia,a,b-negative")
print("g2 as stored:", {k: fshow(v) for k,v in g2.items()})
# verify g2 numerically as a function
def g2num(y):
    return (mp.log(1-1j*a*y)+mp.log(1+1j*a*y)+mp.log(1-1j*b*y)+mp.log(1+1j*b*y)
            -2*mp.log(1-b*y)-2*mp.log(1+a*y))
def lexpr_val(E,y):
    tot=mp.mpc(0)
    for k,c in E.items():
        cc=femb(c,1)
        if k[0]=='C': tot+=cc*mp.log(2)
        elif feq(k[1],ZERO): tot+=cc*mp.log(y)
        else: tot+=cc*mp.log(1-y/femb(k[1],1))
    return tot
for y in [mp.mpf('0.2'), mp.mpf('0.7'), mp.mpf('0.95')]:
    print("g2 check y=%s: expr %s direct %s"%(y, mp.nstr(lexpr_val(g2,y),10), mp.nstr(g2num(y),10)))
