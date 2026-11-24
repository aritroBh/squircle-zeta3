import math
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); b=acb(0,1)*pi/4; b=(b).exp()  # e^{i pi/4}
ialpha=acb(0,2*pi/3)  # i*alpha
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
alpha=2*pi/3
def cubic(t,N):
    t=arb(t); N=arb(N)
    return t**3 - rt2*(arb("0.5")-alpha)*t**2 - alpha*(N-arb("0.5")-alpha)*t - rt2*(N-1)*alpha**2
def find_t0(N, nbisect=60):
    s=(arb(N)*alpha).sqrt(); slo=arb(s).lower(); shi=arb(s).upper()
    # scan over [slo, 1.03*shi]
    lo,hi=None,None; prev=cubic(slo,N); x=slo
    for i in range(1,300):
        x2=arb(slo)+(arb(shi)*arb("1.03")-arb(slo))*arb(i)/299
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
def Phi(t,N): return -b*t+(N-1)*(b*t).log()-(N-arb("1.5"))*(b*t+ialpha).log()
def Phip(t,N): return -b+(N-1)/t-(N-arb("1.5"))*b/(b*t+ialpha)
def Phipp(t,N): return -(N-1)/t**2+(N-arb("1.5"))*b**2/(b*t+ialpha)**2
def RePhipp(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return -(N-1)/t**2+(N-arb("1.5"))*(t**2+rt2*alpha*t)/q**2
def RePhi(t,N): return -t/rt2+(N-1)*t.log()-(N-arb("1.5"))/2*(t**2+rt2*alpha*t+alpha**2).log()
def RePhip(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return -1/rt2+(N-1)/t-(N-arb("1.5"))*(t+alpha/rt2)/q
def Phi3b(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return (t**3+6*alpha*(N-1)*t**2+6*alpha**2*(N-1)*t+2*alpha**3*(N-1))/(t**3*q**arb("1.5"))
def Phi4b(t,N):
    q=t**2+rt2*alpha*t+alpha**2
    return (3*t**4+24*alpha*(N-1)*t**3+36*alpha**2*(N-1)*t**2+24*alpha**3*(N-1)*t+6*alpha**4*(N-1))/(t**4*q**2)
def gmoment(n,p2,om):
    # complex gaussian moments, om real
    n=int(n)
    P=[acb(1)]
    for kk in range(1,n+1):
        P.append((om*P[kk-1]+(kk-1)*P[kk-2] if kk>=2 else om*P[kk-1])/p2)
    return (-acb(0,1))**n*(2*pi/(-p2)).sqrt()*(om**2/(2*p2)).exp()*P[n]
def mom(j,c): return (arb(j+1)/2).gamma()/c**((arb(j)+1)/2)
def f(x): return float(arb(x).mid())

# verify on a thin cell N=400
N=iv("400","400")
t0=find_t0(N)
om=acb(Phip(t0,N)).imag(); p2=Phipp(t0,N); ps=Phi(t0,N); c2=-RePhipp(t0,N)
print("t0:",f(arb(t0).mid()),"om:",f(om),"c2:",f(c2),"Reps:",f(acb(ps).real()))
pre=(acb(2*(2*pi/3)*((2*pi/3).sin())).sqrt())*(acb(0,1)*(pi/4-(2*pi/3)/2)).exp()*(acb(0,1)/pi)*(arb("1.5").gamma())*(N-arb("1.5")).gamma()/N.gamma()
# careful: gamma_+ = sqrt(2 sin a) e^{i(pi/4 - a/2)}, pre = gamma_+ * (i/pi) * Gamma(3/2) * Gamma(N-3/2)/Gamma(N)
gplus=(2*(alpha).sin()).sqrt()*(acb(0,1)*(pi/4-alpha/2)).exp()
pre=gplus*(acb(0,1)/pi)*(arb("1.5").gamma())*(N-arb("1.5")).gamma()/N.gamma()
MV0=pre*(ps).exp()*gmoment(0,p2,om)*b
print("2Re MV0:",f(acb(MV0).real()))
# lead
th=(2*alpha*N).sqrt(); lead=-(acb(2*alpha*alpha.sin()).sqrt()/N)*(-th).exp()*(th).cos()
print("lead:",f(acb(lead).real()))
