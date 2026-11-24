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
def RePhi(t,N):
    u=rt2*alpha/t + alpha**2/t**2
    return -t/rt2 + arb("0.5")*t.log() - (N-arb("1.5"))/2*(1+u).log()
def ps_real_cv(Ncell):
    # Re Phi(t0(N),N) at center + variation via d/dN = log t0 - (1/2) log(t0^2+rt2 a t0+a^2)
    n1=arb(Ncell).lower(); n2=arb(Ncell).upper(); nc=(n1+n2)/2
    t0c=find_t0(nc)
    u=rt2*alpha/t0c + alpha**2/t0c**2
    psc=-t0c/rt2 + arb("0.5")*t0c.log() - (nc-arb("1.5"))/2*(1+u).log()
    # derivative d(ps)/dN = log t0 - (1/2) log(t0^2+rt2 a t0+a^2)
    deriv=t0c.log() - arb("0.5")*(t0c**2+rt2*alpha*t0c+alpha**2).log()
    dN=n2-n1
    return iv(psc-abs(deriv)*dN, psc+abs(deriv)*dN)
import mpmath as mp
mp.mp.dps=30
bmp=mp.e**(1j*mp.pi/4)
def RePhi_m(t,N): return -t/mp.sqrt(2)+(N-1)*mp.log(t)-(N-mp.mpf('1.5'))/2*mp.log(t**2+mp.sqrt(2)*2.0944*t+2.0944**2)
def t0_m(N):
    return float(mp.findroot(lambda t: -1/mp.sqrt(2)+(N-1)/t-(N-mp.mpf('1.5'))*(t+2.0944/mp.sqrt(2))/(t**2+mp.sqrt(2)*2.0944*t+2.0944**2), mp.sqrt(N*2.0944)))
for (n1,n2) in [(400,404),(1600,1616)]:
    Ncell=iv(str(n1),str(n2)); pr=ps_real_cv(Ncell)
    # point values
    pts=[RePhi_m(t0_m(n),n) for n in [n1,(n1+n2)//2,n2]]
    print(f"[{n1},{n2}]: ps_real_cv=[{float(arb(pr).lower()):.4f},{float(arb(pr).upper()):.4f}] points={[round(p,3) for p in pts]}")
