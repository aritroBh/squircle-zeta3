import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); alpha=2*pi/3
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
def RePhip2(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    num=arb("0.5")*t**2 + alpha*((N-1)/rt2 + arb("0.5")/rt2)*t + (N-1)*alpha**2
    return -1/rt2 + num/(t*q)
# test on thin cells
for (n1,n2) in [(400,404),(400,408),(400,420)]:
    Ncell=iv(str(n1),str(n2))
    # xr = t0(N)+delta; use t0 range
    s1=(arb(n1)*alpha).sqrt(); s2=(arb(n2)*alpha).sqrt()
    # t0 ~ 1.008 sqrt(Na); use xr range
    xr=iv(s1*arb("1.005"),s2*arb("1.01"))
    rp=RePhip2(xr,Ncell)
    print(f"[{n1},{n2}]: RePhip2(xr)=[{float(arb(rp).lower()):.4f},{float(arb(rp).upper()):.4f}]  sign-def neg? {arb(rp).upper()<0}")
