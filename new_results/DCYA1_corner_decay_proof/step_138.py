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
def Phi(t,N): return -b*t - (N-1)*(1+ialpha/(b*t)).log() + arb("0.5")*(b*t+ialpha).log()
def Phip(t,N): return -b + (N-1)*ialpha/(t*(b*t+ialpha)) + arb("0.5")*b/(b*t+ialpha)
def Phipp(t,N): return (N-1)*(-2*ialpha*b*t-ialpha**2*-1)/(t**2*(b*t+ialpha)**2) + arb("0.5")*b**2/(b*t+ialpha)**2
def RePhipp(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return -(N-1)/t**2+(N-arb("1.5"))*(t**2+rt2*alpha*t)/q**2
def RePhi(t,N):
    u=rt2*alpha/t + alpha**2/t**2
    return -t/rt2 + arb("0.5")*t.log() - (N-arb("1.5"))/2*(1+u).log()
def RePhip(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    num=arb("0.5")*t**2 + alpha*((N-1)/rt2 + arb("0.5")/rt2)*t + (N-1)*alpha**2
    return -1/rt2 + num/(t*q)
Ncell=iv("1600","1616"); t0=find_t0_range(Ncell)
print("t0:",t0)
om=Phip(t0,Ncell).imag; p2=Phipp(t0,Ncell); ps=Phi(t0,Ncell); c2=-RePhipp(t0,Ncell)
print("om:",om,"c2:",c2)
sig=1/c2.sqrt()
print("sig:",sig)
w=(Ncell.log()+arb("1.0")).sqrt(); delta=w*sig
print("delta:",delta," t03:",t0/3)
deltamin=iv(min(arb(delta).lower(),arb(t03).lower()), min(arb(delta).upper(),arb(t03).upper()))
print("deltamin:",deltamin)
xr=t0+deltamin; xl=t0-deltamin
print("RePhip(xr):",RePhip(xr,Ncell))
print("RePhip(xl):",RePhip(xl,Ncell))
