import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp, time
mp.mp.dps=90
r2n=mp.sqrt(2); an=1+r2n; bn=r2n-1
def fembc(x): return femb(x,1)
ge=GEval(90, M=48, gln=4, subn=4)
def acc_value(acc):
    tot=mp.mpc(0)
    for (w,lvm),c in acc.W.items():
        wn=tuple(fembc(s) for s in w)
        v=ge.word(wn)
        for f in lvm: v*=mp.log(fembc(f[1]))
        tot+=fembc(c)*v
    for lvm,c in acc.P.items():
        v=mp.mpc(1)
        for f in lvm: v*=mp.log(fembc(f[1]))
        tot+=fembc(c)*v
    tot+=fembc(acc.R)
    return tot
# P = log(1-by)*log(1-iay): letters a (c=b), -ib (c=ia)
P = le_prod({Gf(sa):ONE},{Gf(smib):ONE})
dP = prod_deriv(P)
G2=Acc(); int_formP(G2, fmul(a_,a_), sa, dP); boundary_P_at1(G2, fdiv(ONE,b_), P)
d2 = mp.quad(lambda y: mp.log(1-bn*y)*mp.log(1-1j*an*y)/(1-bn*y)**2, [0,1])
print("G2(P) engine:", mp.nstr(acc_value(G2),18))
print("G2(P) direct:", mp.nstr(d2,18))
# check the pieces of the engine output
print("words:", {k:fshow(v) for k,v in G2.W.items()})
print("prods:", {k:fshow(v) for k,v in G2.P.items()})
print("R:", fshow(G2.R))
