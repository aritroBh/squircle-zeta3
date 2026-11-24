import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0):
    return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0):
    return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0):
    return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def mainterm_true(k,a):
    N=k+1; gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a); J=mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
    return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*J
def saddle_ray(N,a,j=0):
    # max of Re Phi on real ray: solve Re Phip=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j))
    return t0,om
def model_value(k,a,j=0):
    N=k+1; t0,om=saddle_ray(N,a,j)
    ps=Phi(t0,N,a,j); p2=Phipp(t0,N,a,j)
    # model with linear phase: int e^{ps + .5 p2 u^2 + i om u} du = e^{ps} sqrt(2pi/(-p2)) e^{om^2/(2 p2)} * ... 
    # exact complex gaussian: int e^{.5 p2 u^2 + i om u} du = sqrt(2 pi/(-p2)) e^{(i om)^2/(2(-p2))}... 
    # use formula: int_RR e^{.5 p2 u^2 + i om u} du = sqrt(2 pi/(-p2)) * e^{om^2/(2 p2)} (check sign)
    MV=mp.e**ps*mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))*b
    gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*MV, t0, om, p2
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [400,1600]:
    N=k+1; th=mp.sqrt(2*a*N); lead=-(A/N)*mp.e**(-th)*mp.cos(th); sc=N**(-1.5)*mp.e**(-th)
    I0=mainterm_true(k,a); mv,t0,om,p2=model_value(k,a)
    print(f"k={k}: t0={mp.nstr(t0,7)} om={mp.nstr(om,5)} |2ReI0-2ReMV|/sc={mp.nstr(abs(2*mp.re(I0)-2*mp.re(mv))/sc,4)} |2ReMV-lead|/sc={mp.nstr(abs(2*mp.re(mv)-lead)/sc,4)}")
