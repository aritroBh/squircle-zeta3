import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp
mp.mp.dps=90
r2n=mp.sqrt(2); an=1+r2n; bn=r2n-1
def fembc(x): return femb(x,1)
ge=GEval(90, M=48, gln=4, subn=4)
P = le_prod({Gf(sa):ONE},{Gf(smib):ONE})
dP = prod_deriv(P)
G2=Acc(); int_formP(G2, fmul(a_,a_), sa, dP); boundary_P_at1(G2, fdiv(ONE,b_), P)
tot=mp.mpc(0)
for (w,lvm),c in G2.W.items():
    wn=tuple(fembc(s) for s in w)
    v=ge.word(wn); tot+=fembc(c)*v
    print("word", ["0" if feq(s,ZERO) else mp.nstr(fembc(s),3) for s in w], "coeff", mp.nstr(fembc(c),6), "G:", mp.nstr(v,10))
for lvm,c in G2.P.items():
    v=mp.mpc(1)
    for f in lvm: v*=mp.log(fembc(f[1]))
    tot+=fembc(c)*v
    print("prod coeff", mp.nstr(fembc(c),6), " factors:", [mp.nstr(fembc(f[1]),4) for f in lvm], " value:", mp.nstr(v,10))
print("engine total:", mp.nstr(tot,14))
# direct independent computation of the same identity: -(1/b) int P' kb + (1/b)[P]01
t1 = mp.quad(lambda y: (mp.log(1-bn*y)*(-1j*an/(1-1j*an*y)) + (-bn/(1-bn*y))*mp.log(1-1j*an*y))*(-an/(y-an)), [0,1])
bd = mp.log(1-bn)*mp.log(1-1j*an)
print("-(1/b)*intP'kb + (1/b)bd =", mp.nstr(-(1/bn)*t1 + (1/bn)*bd, 14))
d2 = mp.quad(lambda y: mp.log(1-bn*y)*mp.log(1-1j*an*y)/(1-bn*y)**2, [0,1])
print("direct:", mp.nstr(d2,14))
