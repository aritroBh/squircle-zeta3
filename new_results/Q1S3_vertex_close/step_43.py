import importlib, q3_engine, q1s_exact, q3_eval
importlib.reload(q1s_exact); importlib.reload(q3_engine); importlib.reload(q3_eval)
from q3_engine import *
from q3_eval import *
import mpmath as mp
def fembc(x): return femb(x,1)
mp.mp.dps=130
r2=mp.sqrt(2)
exact = []
for w,val in [
    ((0,1j*(r2-1)), -mp.polylog(2,-1j*(1+r2))),
    ((0,1j*(1+r2)), -mp.polylog(2,-1j*(r2-1))),
    ((0,r2-1), -mp.polylog(2,1+r2)),          # G(0,-b;1) = -Li2(1/(-b)) = -Li2(-a)... careful
]:
    pass
# build exact cases properly: G(0,v;1) = -Li2(1/v); G(s,s;1)=1/2 log(1-1/s)^2; G(0,v,v;1) = ?
cases = []
for v,name in [(1j*(r2-1),'ib'),(-1j*(r2-1),'-ib'),(1+r2,'a'),-(r2-1),'-b' if False else '-b',(1j*(1+r2),'ia'),(-1j*(1+r2),'-ia')]:
    cases.append(((mp.mpc(0),mp.mpc(v)), -mp.polylog(2,1/v), 'G(0,%s)'%name))
    cases.append(((mp.mpc(v),mp.mpc(v)), mp.log(1-1/v)**2/2, 'G(%s,%s)'%(name,name)))
for h,T in [('0.012','4.6'),('0.008','4.8'),('0.006','5.0')]:
    ge=GEval(130,h=h,T=T,gln=4)
    mx=0
    for w,val,_ in cases:
        mx=max(mx,abs(ge.word(w)-val))
    print("h=%s: max resid vs exact: %s   (N=%d)"%(h, mp.nstr(mx,4), ge.N))
