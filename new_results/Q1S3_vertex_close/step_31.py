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
    return mp.quad(lambda t: Grec(w[1:],t)/(t-s), [0,z])
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
# int_P_dy consistency: int P' dy = [P]_0^1
acc5=Acc(); int_P_dy(acc5, ONE, dt)
acc6=Acc(); boundary_P_at1(acc6, ONE, P)
print("int_P_dy vs boundary:", mp.nstr(acc_value(acc5),15), mp.nstr(acc_value(acc6),15))
