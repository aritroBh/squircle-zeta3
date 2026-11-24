src = open('q3_engine.py').read()
old = """            A=fdiv(u,fsub(u,s)); Bm=fneg(fdiv(s,fsub(u,s)))
            acc.add_prod((lv_G1(u),), fmul(fdiv(cc,s),A))
            acc.add_prod((lv_G1(s),), fmul(fdiv(cc,s),fneg(Bm)))"""
new = """            A=fdiv(u,fsub(u,s)); Bm=fneg(fdiv(s,fsub(u,s)))
            acc.add_prod((lv_G1(u),), fmul(fdiv(cc,s),A))
            acc.add_prod((lv_G1(s),), fmul(fdiv(cc,s),Bm))"""
assert old in src
open('q3_engine.py','w').write(src.replace(old,new))
print("fixed")

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
        v=ge.word(tuple(fembc(s) for s in w))
        for f in lvm: v*=mp.log(fembc(f[1]))
        tot+=fembc(c)*v
    for lvm,c in acc.P.items():
        v=mp.mpc(1)
        for f in lvm: v*=mp.log(fembc(f[1]))
        tot+=fembc(c)*v
    tot+=fembc(acc.R)
    return tot
# single term check vs the G2_X-form identity -(1/b)int P'kb + (1/b)P(1):
P = le_prod({Gf(sa):ONE},{Gf(smib):ONE})
dP = prod_deriv(P)
G2=Acc(); int_formP(G2, fmul(a_,a_), sa, dP); boundary_P_at1(G2, fdiv(ONE,b_), P)
t1 = mp.quad(lambda y: (mp.log(1-bn*y)*(-1j*an/(1-1j*an*y)) + (-bn/(1-bn*y))*mp.log(1-1j*an*y))*(-an/(y-an)), [0,1])
bd = mp.log(1-bn)*mp.log(1-1j*an)
target = -(1/bn)*t1 + (1/bn)*bd
print("single term: engine", mp.nstr(acc_value(G2),15))
print("             target", mp.nstr(target,15))
