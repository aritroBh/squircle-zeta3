import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from flint import arb, acb, ctx
ctx.dps=60
from arb_ei import (iv, find_t0_range, Phi, Phip, Phipp, c2f, RePhi, RePhip, Phi3b, Phi4b, gmoment, mom, preF, ps_real_cv, tail_RePhi_cv, alpha, pi, rt2, b)
def f(x): return float(arb(x).mid())
Ncell=iv("20000","20142")
J=48
t0=find_t0_range(Ncell)
print("t0:",t0)
om=Phip(t0,Ncell).imag; p2=Phipp(t0,Ncell); c2=c2f(t0,Ncell)
print("om:",f(om),"c2:",f(c2),"sig:",f(1/c2.sqrt()))
sig=1/c2.sqrt(); w=(Ncell.log()+arb("1.0")).sqrt(); delta=w*sig
t03=t0/3
print("w*sig:",f(delta)," t0/3:",f(t03))
delta=iv(min(arb(delta).lower(),arb(t03).lower()), min(arb(delta).upper(),arb(t03).upper()))
print("delta used:",f(delta)," (t0/3 capped?)", f(delta)==f(t03))
pre=preF(Ncell); P0=abs(pre)*(ps_real_cv(Ncell)).exp()
print("P0:",P0)
W1=P0*abs(om)*(Phi3b(t0,Ncell)/6*mom(4,c2))
print("W1:",W1)
W2=arb(0); W3=arb(0)
for side in [1,-1]:
    for jj in range(J):
        u0=delta*jj/J; u1=delta*(jj+1)/J
        xsegL=t0-u1 if side==-1 else t0
        M3=Phi3b(xsegL,Ncell); M4=Phi4b(xsegL,Ncell)
        R4b=M4*u1**4/24; V=M3*u1**3/6; gw=(-c2*u0**2/2).exp(); du=u1-u0
        W2=W2+P0*gw*R4b*du; W3=W3+P0*gw*(V*V/2)*du
print("W2:",W2," W3:",W3)
reR,rpR=tail_RePhi_cv(Ncell,1); reL,rpL=tail_RePhi_cv(Ncell,-1)
TR=abs(pre)*(reR).exp()/abs(rpR)
print("TR:",TR)
TLm=abs(pre)*(reL).exp()/rpL
print("TLm:",TLm)
x=delta/sig; TM=P0*2*sig*(-x**2/2).exp()/x
print("TM:",TM)
