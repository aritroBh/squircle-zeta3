import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
b_mp=mp.e**(1j*mp.pi/4)
def Phi_m(t,N): return -b_mp*t+(N-1)*mp.log(b_mp*t)-(N-mp.mpf('1.5'))*mp.log(b_mp*t+1j*a)
def Phipp_m(t,N): return -(N-1)/t**2+(N-mp.mpf('1.5'))*b_mp**2/(b_mp*t+1j*a)**2
def Phip_m(t,N): return -b_mp+(N-1)/t-(N-mp.mpf('1.5'))*b_mp/(b_mp*t+1j*a)
def gmoment0_m(p2,om): return mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))
a=2*mp.pi/3; N=400
t0=mp.findroot(lambda t: mp.re(Phip_m(t,N)), mp.sqrt(N*a))
om=mp.im(Phip_m(t0,N)); p2=Phipp_m(t0,N); ps=Phi_m(t0,N)
pre_mp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
MV0_mp=pre_mp*mp.e**ps*gmoment0_m(p2,om)*b_mp
print("mpmath MV0 =", mp.nstr(MV0_mp,10))
# arb
from flint import arb, acb, ctx
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); b=(acb(0,1)*pi/4).exp(); ialpha=acb(0,2*pi/3); alpha=2*pi/3
def Phi(t,N): return -b*t+(N-1)*(b*t).log()-(N-arb("1.5"))*(b*t+ialpha).log()
def Phipp(t,N): return -(N-1)/t**2+(N-arb("1.5"))*b**2/(b*t+ialpha)**2
def Phip(t,N): return -b+(N-1)/t-(N-arb("1.5"))*b/(b*t+ialpha)
def gmoment(n,p2,om):
    n=int(n); P=[acb(1)]
    for kk in range(1,n+1): P.append((om*P[kk-1]+(kk-1)*P[kk-2] if kk>=2 else om*P[kk-1])/p2)
    return (-acb(0,1))**n*(2*pi/(-p2)).sqrt()*(om**2/(2*p2)).exp()*P[n]
gplus=(2*alpha.sin()).sqrt()*(acb(0,1)*(pi/4-alpha/2)).exp()
Nar=arb("400"); t0a=arb("29.178121184894174")
oma=Phip(t0a,Nar).imag; p2a=Phipp(t0a,Nar); psa=Phi(t0a,Nar)
prea=gplus*(acb(0,1)/pi)*(arb("1.5").gamma())*(Nar-arb("1.5")).gamma()/Nar.gamma()
MV0_arb=prea*(psa).exp()*gmoment(0,p2a,oma)*b
print("arb MV0 =", MV0_arb)
