import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); b=(acb(0,1)*pi/4).exp(); ialpha=acb(0,2*pi/3); alpha=2*pi/3
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
def cubic(t,N):
    t=arb(t); N=arb(N)
    return t**3 - rt2*(arb("0.5")-alpha)*t**2 - alpha*(N-arb("0.5")-alpha)*t - rt2*(N-1)*alpha**2
def find_t0(N, nbisect=50):
    s=(arb(N)*alpha).sqrt(); slo=arb(s).lower(); shi=arb(s).upper()
    lo,hi=None,None; prev=cubic(slo,N); x=slo
    for i in range(1,200):
        x2=arb(slo)+(arb(shi)*arb("1.03")-arb(slo))*arb(i)/199
        cur=cubic(x2,N)
        if arb(prev).lower()<=0<=arb(cur).upper() or arb(prev).upper()>=0>=arb(cur).lower():
            lo,hi=x,x2; break
        prev=cur; x=x2
    if lo is None: return iv(slo,shi*arb("1.03"))
    for _ in range(nbisect):
        mid=(lo+hi)/2; fm=cubic(mid,N); flo=cubic(lo,N)
        if (arb(flo).lower()<=0<=arb(fm).upper()) or (arb(flo).upper()>=0>=arb(fm).lower()): hi=mid
        else: lo=mid
    return iv(lo,hi)
def find_t0_range(Ncell):
    n1=arb(Ncell).lower(); n2=arb(Ncell).upper()
    return iv(find_t0(n1).lower(), find_t0(n2).upper())
def RePhipp(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return -(N-1)/t**2+(N-arb("1.5"))*(t**2+rt2*alpha*t)/q**2
def RePhip(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return -1/rt2+(N-1)/t-(N-arb("1.5"))*(t+alpha/rt2)/q
def Phip(t,N): return -b + (N-1)*ialpha/(t*(b*t+ialpha)) + arb("0.5")*b/(b*t+ialpha)
for (n1,n2) in [(400,404),(400,408),(400,420),(1600,1616)]:
    Ncell=iv(str(n1),str(n2)); t0=find_t0_range(Ncell)
    c2=-RePhipp(t0,Ncell); om=Phip(t0,Ncell).imag
    xr=t0+3; xl=t0-3
    print(f"[{n1},{n2}]: t0w={float(arb(t0).upper()-arb(t0).lower()):.3f} c2=[{float(arb(c2).lower()):.4f},{float(arb(c2).upper()):.4f}] om=[{float(arb(om).lower()):.4f},{float(arb(om).upper()):.4f}] RePhip(xr)=[{float(arb(RePhip(xr,Ncell)).lower()):.3f},{float(arb(RePhip(xr,Ncell)).upper()):.3f}]")
