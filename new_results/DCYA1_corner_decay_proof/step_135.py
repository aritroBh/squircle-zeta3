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
def RePhi(t,N): return -t/rt2+(N-1)*t.log()-(N-arb("1.5"))/2*(t**2+rt2*alpha*t+alpha**2).log()
def RePhip(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    num=arb("0.5")*t**2 + alpha*((N-1)/rt2 + arb("0.5")/rt2)*t + (N-1)*alpha**2
    return -1/rt2 + num/(t*q)
def Phi3b(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return (t**3+6*alpha*(N-1)*t**2+6*alpha**2*(N-1)*t+2*alpha**3*(N-1))/(t**3*q**arb("1.5"))
def Phi4b(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return (3*t**4+24*alpha*(N-1)*t**3+36*alpha**2*(N-1)*t**2+24*alpha**3*(N-1)*t+6*alpha**4*(N-1))/(t**4*q**2)
def gmoment(n,p2,om):
    n=int(n); P=[acb(1)]
    for kk in range(1,n+1): P.append((om*P[kk-1]+(kk-1)*P[kk-2] if kk>=2 else om*P[kk-1])/p2)
    return (-acb(0,1))**n*(2*pi/(-p2)).sqrt()*(om**2/(2*p2)).exp()*P[n]
def mom(j,c): return (arb(j+1)/2).gamma()/c**((arb(j)+1)/2)
gplus=(2*alpha.sin()).sqrt()*(acb(0,1)*(pi/4-alpha/2)).exp()
def gamma_ratio_interval(Ncell):
    n1=arb(Ncell).lower(); n2=arb(Ncell).upper(); nc=(n1+n2)/2
    lg1=(nc-arb("1.5")).lgamma(); lg2=nc.lgamma()
    logr_c=lg1-lg2; dn1=arb(n1).lower()
    dg=arb("1.6")/(dn1-arb("1.5")); dlogr=dg*(n2-n1)
    return iv(logr_c-dlogr, logr_c+dlogr).exp()
def preF(Ncell): return gplus*(acb(0,1)/pi)*(arb("1.5").gamma())*gamma_ratio_interval(Ncell)
def f(x): return float(arb(x).mid())
def fw(x): return float(arb(x).upper()-arb(x).lower())
Ncell=iv("400","404"); J=48
t0=find_t0_range(Ncell)
om=Phip(t0,Ncell).imag; p2=Phipp(t0,Ncell); ps=Phi(t0,Ncell); c2=-RePhipp(t0,Ncell)
sig=1/c2.sqrt(); w=(Ncell.log()+arb("1.0")).sqrt(); delta=w*sig
t03=t0/3; delta=iv(min(arb(delta).lower(),arb(t03).lower()), min(arb(delta).upper(),arb(t03).upper()))
pre=preF(Ncell); P0=abs(pre)*(ps.real).exp()
print("P0:",f(P0),"width:",fw(P0))
W1=P0*abs(om)*(Phi3b(t0,Ncell)/6*mom(4,c2))
print("W1:",f(W1),"width:",fw(W1))
W2=arb(0); W3=arb(0)
for side in [1,-1]:
    for jj in range(J):
        u0=delta*jj/J; u1=delta*(jj+1)/J
        xsegL=t0-u1 if side==-1 else t0
        M3=Phi3b(xsegL,Ncell); M4=Phi4b(xsegL,Ncell)
        R4b=M4*u1**4/24; V=M3*u1**3/6; gw=(-c2*u0**2/2).exp(); du=u1-u0
        W2=W2+P0*gw*R4b*du; W3=W3+P0*gw*(V*V/2)*du
print("W2:",f(W2),"width:",fw(W2)," W3:",f(W3),"width:",fw(W3))
xr=t0+delta; xl=t0-delta
TR=abs(pre)*(RePhi(xr,Ncell)).exp()/abs(RePhip(xr,Ncell))
print("TR:",f(TR),"width:",fw(TR))
TLm=abs(pre)*(RePhi(xl,Ncell)).exp()/RePhip(xl,Ncell)
print("TLm:",f(TLm),"width:",fw(TLm))
x=delta/sig; TM=P0*2*sig*(-x**2/2).exp()/x
print("TM:",f(TM),"width:",fw(TM))
