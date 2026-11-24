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
def Phi(t,N): return -b*t+(N-1)*(b*t).log()-(N-arb("1.5"))*(b*t+ialpha).log()
def Phip(t,N): return -b+(N-1)/t-(N-arb("1.5"))*b/(b*t+ialpha)
def Phipp(t,N): return -(N-1)/t**2+(N-arb("1.5"))*b**2/(b*t+ialpha)**2
def RePhipp(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return -(N-1)/t**2+(N-arb("1.5"))*(t**2+rt2*alpha*t)/q**2
N=iv("400","420"); t0=find_t0_range(N)
print("t0:", arb(t0).lower(), arb(t0).upper())
ps=Phi(t0,N)
print("ps:", ps)
print("ps.real:", ps.real)
print("(ps.real).exp():", (ps.real).exp())
om=Phip(t0,N).imag
print("om:", om)
p2=Phipp(t0,N)
print("p2:", p2)
c2=-RePhipp(t0,N)
print("c2:", c2)
print("sig:", 1/c2.sqrt())
