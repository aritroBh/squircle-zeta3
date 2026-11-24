import math, sys
from flint import arb, acb, ctx
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
ctx.dps=60
rt2=arb(2).sqrt(); pi=arb.pi(); b=(acb(0,1)*pi/4).exp(); ialpha=acb(0,2*pi/3)
alpha=2*pi/3
def Phi(t,N): return -b*t+(N-1)*(b*t).log()-(N-arb("1.5"))*(b*t+ialpha).log()
def Phipp(t,N): return -(N-1)/t**2+(N-arb("1.5"))*b**2/(b*t+ialpha)**2
def Phip(t,N): return -b+(N-1)/t-(N-arb("1.5"))*b/(b*t+ialpha)
def gmoment(n,p2,om):
    n=int(n); P=[acb(1)]
    for kk in range(1,n+1): P.append((om*P[kk-1]+(kk-1)*P[kk-2] if kk>=2 else om*P[kk-1])/p2)
    return (-acb(0,1))**n*(2*pi/(-p2)).sqrt()*(om**2/(2*p2)).exp()*P[n]
def f(x): return float(arb(x).mid())
gplus=(2*alpha.sin()).sqrt()*(acb(0,1)*(pi/4-alpha/2)).exp()
N=arb("400"); t0=arb("29.178121184894174")
om=Phip(t0,N).imag; p2=Phipp(t0,N); ps=Phi(t0,N)
pre=gplus*(acb(0,1)/pi)*(arb("1.5").gamma())*(N-arb("1.5")).gamma()/N.gamma()
# mpmath
bmp=mp.e**(1j*mp.pi/4)
def Phipp_m(t): return -(399)/t**2+(399-mp.mpf('0.5'))*bmp**2/(bmp*t+1j*alpha)**2
def Phip_m(t): return -bmp+(399)/t-(399-mp.mpf('0.5'))*bmp/(bmp*t+1j*alpha)
p2m=Phipp_m(29.178121184894174); omm=mp.im(Phip_m(29.178121184894174))
def gmoment0_m(p2,om): return mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))
print("om: arb=",f(om)," mp=",float(omm))
print("gmoment0: arb=",f(abs(gmoment(0,p2,om)))," mp=",float(abs(gmoment0_m(p2m,omm))))
print("|pre|: arb=",f(abs(pre)))
gpm=mp.sqrt(2*mp.sin(alpha))*mp.e**(1j*(mp.pi/4-alpha/2))
prem=gpm*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(399-mp.mpf('0.5'))/mp.gamma(399)
print("|pre|: mp=",float(abs(prem)))
print("e^{Re ps}: arb=",f((ps.real).exp()))
def Phi_m(t): return -bmp*t+(399)*mp.log(bmp*t)-(399-mp.mpf('0.5'))*mp.log(bmp*t+1j*alpha)
print("e^{Re ps}: mp=",float(mp.e**mp.re(Phi_m(29.178121184894174))))
