import importlib, q3_engine, q1s_exact
importlib.reload(q1s_exact); importlib.reload(q3_engine)
from q3_engine import *
from q1s_exact import *
import mpmath as mp
mp.mp.dps = 60

# numeric helpers for control checks
def fembc(x): return femb(x,1)
def acc_value(acc, dps=60):
    """numeric value of Acc via direct word evaluation (mp.quad recursion, low order OK for small tests)."""
    mp.mp.dps=dps
    def G(w):
        if len(w)==1:
            s=w[0]
            return mp.log(1-1/fembc(s))
        s=fembc(w[0])
        return mp.quad(lambda t: G(w[1:])/(t-s) if not feq(w[0],ZERO) else G(w[1:])/t, [0,1])
    def lv(x): return mp.log(fembc(x))
    tot=mp.mpc(0)
    for (w,lvm),c in acc.W.items():
        v=G(w)
        for f in lvm: v*=lv(f[1])
        tot+=fembc(c)*v
    for lvm,c in acc.P.items():
        v=mp.mpc(1)
        for f in lvm: v*=lv(f[1])
        tot+=fembc(c)*v
    tot+=fembc(acc.R)
    return tot

b = mp.sqrt(2)-1
# Control 1: int_0^1 log y * log(1-b y)/(1-b y)^2 dy  via engine rule (square kernel)
P = le_prod({Gf(s0):ONE},{Gf(sa):ONE})   # log y * log(1-by): letter of log(1-by) is a
acc=Acc()
# int P/(1-by)^2 = -(1/b) int P' kb + (1/b)[P]_0^1, kb=(cK=-a,sK=a)
int_formP(acc, ONE, ZERO, fdiv(fneg(ONE),b_), sa, [])  # placeholder no
acc=Acc()
dt=prod_deriv(P)
int_formP(acc, ONE,ZERO, fdiv(fneg(ONE),b_), sa, dt)  # wrong signature; fix below
print("skip")
