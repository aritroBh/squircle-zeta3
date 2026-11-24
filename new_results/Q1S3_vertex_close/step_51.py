import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp
DPS=90
mp.mp.dps=DPS
ge=GEval(DPS, M=48, gln=5)
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

r2_ = r2; a_2 = fmul(a_,a_); b_2 = fmul(b_,b_)
hr2 = fscale(r2_, Fr(1,2))
P12 = le_prod(g1,g2); dP12 = prod_deriv(P12)
G1=Acc()
int_simple_kernel(G1, ONE, s0, P12)
int_simple_kernel(G1, hr2, sa, P12)
int_simple_kernel(G1, fneg(hr2), smb, P12)
G2=Acc()
int_formP(G2, a_2, sa, dP12)
boundary_P_at1(G2, fdiv(ONE,b_), P12)
G3=Acc()
int_formP(G3, b_2, smb, dP12)
boundary_P_at1(G3, fneg(fdiv(ONE,a_)), P12)
print("G1 engine:", mp.nstr(acc_value(G1),20), " target -1.34483248386795")
print("G2 engine:", mp.nstr(acc_value(G2),20), " target -0.184917904781547")
print("G3 engine:", mp.nstr(acc_value(G3),20), " target -0.0846038124400114")
