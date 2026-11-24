import importlib, q3_engine, q1s_exact
importlib.reload(q1s_exact); importlib.reload(q3_engine)
from q3_engine import *
from q1s_exact import *
import mpmath as mp
mp.mp.dps = 60
def fembc(x): return femb(x,1)
def lv(x): return mp.log(fembc(x))
def Grec(w, z=1):
    if len(w)==0: return mp.mpc(1)
    s=fembc(w[0])
    if len(w)==1:
        if feq(w[0],ZERO): return mp.log(z)
        return mp.log(1-z/s)
    return mp.quad(lambda t: Grec(w[1:],t)/(t-s), [0,z]) if z!=0 else mp.mpc(0)
def acc_value(acc):
    tot=mp.mpc(0)
    for (w,lvm),c in acc.W.items():
        v=Grec(w)
        for f in lvm: v*=lv(f[1])
        tot+=fembc(c)*v
    for lvm,c in acc.P.items():
        v=mp.mpc(1)
        for f in lvm: v*=lv(f[1])
        tot+=fembc(c)*v
    tot+=fembc(acc.R)
    return tot

b = mp.sqrt(2)-1; a=1+mp.sqrt(2)
P = le_prod({Gf(s0):ONE},{Gf(sa):ONE})
acc=Acc(); dt=prod_deriv(P)
int_formP(acc, fmul(a_,a_), sa, dt)
boundary_P_at1(acc, fdiv(ONE,b_), P)
val_direct = mp.quad(lambda y: mp.log(y)*mp.log(1-b*y)/(1-b*y)**2, [0,1])
print("control1 resid:", mp.nstr(abs(acc_value(acc)-val_direct),5))

P2 = le_prod({Gf(sa):ONE},{Gf(smb):ONE})
acc2=Acc(); int_simple_kernel(acc2, ONE, s0, P2)
val_d2 = mp.quad(lambda y: mp.log(1-b*y)*mp.log(1+a*y)/y, [0,1])
print("control2 resid:", mp.nstr(abs(acc_value(acc2)-val_d2),5))

P3 = le_prod({Gf(s0):ONE},{Gf(smb):ONE})
acc3=Acc(); int_simple_kernel(acc3, fneg(a_), sa, P3)
val_d3 = mp.quad(lambda y: mp.log(y)*mp.log(1+a*y)/(1-b*y), [0,1])
print("control3 resid:", mp.nstr(abs(acc_value(acc3)-val_d3),5))

acc4=Acc()
for f,c in g1.items():
    for lv1 in atom_at1(f):
        acc4.add_prod((lv1,), fdiv(c,fsub(ONE,b_)))
for c,sF in le_deriv_terms(g1):
    int_y_rr(acc4, fmul(sa,c), sF, sa)
val_d4 = mp.quad(lambda y: mp.log((1-b*y)*(1+a*y)/(2*y))/(1-b*y)**2, [0,1])
print("control4 (H4) resid:", mp.nstr(abs(acc_value(acc4)-val_d4),5), " val:", mp.nstr(val_d4,15))
