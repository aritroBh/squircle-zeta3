import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi()
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
alpha=2*pi/3
# saddle cubic: t^3 - rt2(1/2-alpha) t^2 - alpha(N-1/2-alpha) t - rt2(N-1) alpha^2 = 0
def cubic(t,N):
    t=arb(t); N=arb(N)
    return t**3 - rt2*(arb("0.5")-alpha)*t**2 - alpha*(N-arb("0.5")-alpha)*t - rt2*(N-1)*alpha**2
def find_t0(N):
    # enclose the positive root of cubic near sqrt(N alpha)
    s=(arb(N)*alpha).sqrt()
    # scan for sign change of cubic from s to 2s
    import numpy as np
    pts=[s*(1+i*0.01) for i in range(120)]
    prev=cubic(pts[0],N)
    for i in range(1,len(pts)):
        cur=cubic(pts[i],N)
        if arb(prev).lower()<=0<=arb(cur).upper() or arb(prev).upper()>=0>=arb(cur).lower():
            # root in [pts[i-1],pts[i]]
            return iv(pts[i-1],pts[i])
        prev=cur
    return iv(s,2*s)
for N in [100,400,1600]:
    t0=find_t0(N)
    s=(arb(N)*alpha).sqrt()
    print(f"N={N}: t0 interval = {t0}  (width={arb(t0).upper()-arb(t0).lower()}, sqrt(Na)={s}, t0/sqrt={arb(t0).mid()/s})")
    print(f"    cubic(t0 interval) = {cubic(t0,N)}")
