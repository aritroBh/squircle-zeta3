import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi()
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
alpha=2*pi/3
def cubic(t,N):
    t=arb(t); N=arb(N)
    return t**3 - rt2*(arb("0.5")-alpha)*t**2 - alpha*(N-arb("0.5")-alpha)*t - rt2*(N-1)*alpha**2
def find_t0(N, nsteps=200, nbisect=40):
    s=(arb(N)*alpha).sqrt()
    # coarse scan
    lo,hi=None,None
    prev=cubic(s,N)
    for i in range(1,nsteps):
        x=s*(1+i*0.005)
        cur=cubic(x,N)
        if arb(prev).lower()<=0<=arb(cur).upper() or arb(prev).upper()>=0>=arb(cur).lower():
            lo,hi=s*(1+(i-1)*0.005),x; break
        prev=cur
    if lo is None: return iv(s,2*s)
    # bisect
    for _ in range(nbisect):
        mid=(lo+hi)/2
        fm=cubic(mid,N)
        flo=cubic(lo,N)
        if (arb(flo).lower()<=0<=arb(fm).upper()) or (arb(flo).upper()>=0>=arb(fm).lower()):
            hi=mid
        else:
            lo=mid
    return iv(lo,hi)
for N in [100,400,1600]:
    t0=find_t0(N)
    print(f"N={N}: t0=[{arb(t0).lower()}, {arb(t0).upper()}] width={arb(t0).upper()-arb(t0).lower():.2e} cubic={cubic(t0,N)}")
