import importlib, q3_engine, q1s_exact
importlib.reload(q1s_exact); importlib.reload(q3_engine)
from q3_engine import *
from q1s_exact import *
import mpmath as mp
mp.mp.dps = 60
def fembc(x): return femb(x,1)
def lv(x): return mp.log(fembc(x))
def Grec(w):
    """slow recursive G(w;1) for controls."""
    if len(w)==0: return mp.mpc(1)
    s=fembc(w[0])
    if len(w)==1:
        return mp.log(1-1/s)
    return mp.quad(lambda t: Grec(w[1:])/(t-s), [0,1])
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

b = mp.sqrt(2)-1
# Control 1: int_0^1 log y * log(1-b y)/(1-b y)^2 dy
P = le_prod({Gf(s0):ONE},{Gf(sa):ONE})
acc=Acc()
dt=prod_deriv(P)
int_formP(acc, fmul(a_,a_), sa, dt)          # -(1/b)*(-a) = a^2 factor, kernel letter a
boundary_P_at1(acc, fdiv(ONE,b_), P)
val_engine = acc_value(acc)
val_direct = mp.quad(lambda y: mp.log(y)*mp.log(1-b*y)/(1-b*y)**2, [0,1])
print("control1:", mp.nstr(val_engine,20), mp.nstr(val_direct,20), mp.nstr(abs(val_engine-val_direct),5))

# Control 2: int_0^1 [1/y] log(1-by) log(1+ay) dy  (simple kernel k0)
P2 = le_prod({Gf(sa):ONE},{Gf(smb):ONE})
acc2=Acc()
int_simple_kernel(acc2, ONE, s0, P2)
val_e2 = acc_value(acc2)
val_d2 = mp.quad(lambda y: mp.log(1-b*y)*mp.log(1+b... if False else 0, [0,1]) if False else mp.quad(lambda y: mp.log(1-b*y)*mp.log(1+a... ,[0,1]) if False else 0
a=1+mp.sqrt(2)
val_d2 = mp.quad(lambda y: mp.log(1-b*y)*mp.log(1+a*y)/y, [0,1])
print("control2:", mp.nstr(val_e2,20), mp.nstr(val_d2,20), mp.nstr(abs(val_e2-val_d2),5))

# Control 3: simple kernel 1/(1-by): int_0^1 log(y)log(1+ay)/(1-by) dy
P3 = le_prod({Gf(s0):ONE},{Gf(smb):ONE})
acc3=Acc()
int_simple_kernel(acc3, fneg(a_), sa, P3)   # 1/(1-by) = -a/(y-a)
val_e3 = acc_value(acc3)
val_d3 = mp.quad(lambda y: mp.log(y)*mp.log(1+a*y)/(1-b*y), [0,1])
print("control3:", mp.nstr(val_e3,20), mp.nstr(val_d3,20), mp.nstr(abs(val_e3-val_d3),5))
