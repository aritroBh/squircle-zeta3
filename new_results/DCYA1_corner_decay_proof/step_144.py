import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); alpha=2*pi/3
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
def RePhi(t,N):
    u=rt2*alpha/t + alpha**2/t**2
    return -t/rt2 + arb("0.5")*t.log() - (N-arb("1.5"))/2*(1+u).log()
def RePhip(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    num=arb("0.5")*t**2 + alpha*((N-1)/rt2 + arb("0.5")/rt2)*t + (N-1)*alpha**2
    return -1/rt2 + num/(t*q)
# center+variation for ReΦ(xr) and ReΦ'(xr)
def tail_cv(Ncell, side):
    n1=arb(Ncell).lower(); n2=arb(Ncell).upper(); nc=(n1+n2)/2
    # xr = t0(N)+delta(N); approximate delta by w*sig (but we need c2, t0 at center)
    # For the CV, compute xr range via t-derivative of RePhi
    # First need t0(nc), c2(nc), sig(nc), w(nc), delta(nc)
    # (import from a helper) -- here inline:
    return None
# Debug: check interval widths of ReΦ(xr) and ReΦ'(xr) for [1600,1616]
Ncell=iv("1600","1616")
# xr range: t0 range [58.18,58.76], delta ~ w*sig ~ sqrt(logN+1)/sqrt(c2) ~ let me estimate ~ 13
# so xr ~ 58.5+13 = 71.5, xl ~ 58.5-13 = 45.5
xr=iv("70.0","73.0")
rp=RePhip(xr,Ncell)
print("RePhip(xr) interval:",rp,"width:",float(arb(rp).upper()-arb(rp).lower()))
rpv=RePhi(xr,Ncell)
print("RePhi(xr) interval:",rpv,"width:",float(arb(rpv).upper()-arb(rpv).lower()))
