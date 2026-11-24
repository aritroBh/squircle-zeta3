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
def c2f(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    num=arb("0.5")*t**4 + rt2*alpha*(N-arb("0.5"))*t**3 + 4*alpha**2*(N-1)*t**2 + 2*rt2*alpha**3*(N-1)*t + alpha**4*(N-1)
    return num/(t**2*q**2)
def RePhi(t,N):
    u=rt2*alpha/t + alpha**2/t**2
    return -t/rt2 + arb("0.5")*t.log() - (N-arb("1.5"))/2*(1+u).log()
def RePhip(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    num=arb("0.5")*t**2 + alpha*((N-1)/rt2 + arb("0.5")/rt2)*t + (N-1)*alpha**2
    return -1/rt2 + num/(t*q)
# center+variation for ReΦ(xr,N) where xr=t0+delta (side=+1) or t0-delta (side=-1)
def tail_RePhi_cv(Ncell, side):
    n1=arb(Ncell).lower(); n2=arb(Ncell).upper(); nc=(n1+n2)/2
    t0c=find_t0(nc); c2c=c2f(t0c,nc); sigc=1/c2c.sqrt(); wc=(nc.log()+arb("1.0")).sqrt()
    deltac=wc*sigc; t03c=t0c/3
    deltac=iv(min(arb(deltac).lower(),arb(t03c).lower()), min(arb(deltac).upper(),arb(t03c).upper()))
    xrc=t0c+side*deltac
    # value at center
    u=rt2*alpha/xrc + alpha**2/xrc**2
    val=-xrc/rt2 + arb("0.5")*xrc.log() - (nc-arb("1.5"))/2*(1+u).log()
    # t-variation: |RePhi'(xr)| * Delta(xr over cell). Delta xr ~ (dt0+ddelta) ~ from t0_range and delta range
    rp=RePhip(xrc,nc)
    # N-variation of the -(N-1.5)/2 log(1+u) term: |-(1/2) log(1+u)| * dN
    dNdN=arb("0.5")*(1+u).log()
    dN=n2-n1
    # delta xr over cell: from t0 range and c2 range
    t0r=find_t0_range(Ncell)
    c2r_lo=c2f(find_t0(n1),n1); c2r_hi=c2f(find_t0(n2),n2)
    sig1=1/c2r_lo.sqrt(); sig2=1/c2r_hi.sqrt()
    w1=(n1.log()+arb("1.0")).sqrt(); w2=(n2.log()+arb("1.0")).sqrt()
    dxr=(arb(t0r).upper()-arb(t0r).lower()) + (w2*sig2-w1*sig1)
    var=abs(rp)*abs(dxr) + dNdN*dN
    return iv(val-var, val+var)
def f(x): return float(arb(x).mid())
def fw(x): return float(arb(x).upper()-arb(x).lower())
for (n1,n2,side) in [(400,404,1),(1600,1616,1),(400,404,-1),(1600,1616,-1)]:
    Ncell=iv(str(n1),str(n2))
    tr=tail_RePhi_cv(Ncell,side)
    print(f"[{n1},{n2}] side={side}: RePhi(tail)=[{f(arb(tr).lower()):.3f},{f(arb(tr).upper()):.3f}] width={fw(tr):.3f}")
