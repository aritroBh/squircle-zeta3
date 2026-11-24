import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def mainterm_I0(k,a):
    N=k+1; gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a); J=mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
    return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*J
def gauss_moment0(p2,om):
    return mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))
def exact_window_bound(k,a):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j); Reps=mp.re(ps)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    # model value
    MV=pre*mp.e**ps*gauss_moment0(p2,om)*b
    # exact window bound: |pre| int e^{Re Phi_m(u)} |e^{R(u)}-1| du
    def integrand(u):
        Phim=ps+0.5*p2*u**2+1j*om*u
        R=Phi(t0+u,N,a,j)-Phim
        return mp.e**mp.re(Phim)*abs(mp.e**R-1)
    W=abs(pre)*mp.quad(integrand,[-6/mp.sqrt(-mp.re(p2)),0,6/mp.sqrt(-mp.re(p2))])
    return MV,W
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [400,1600,6400]:
    N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th); lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    MV,W=exact_window_bound(k,a); I0=mainterm_I0(k,a)
    true_res=abs(2*mp.re(I0)-2*mp.re(MV))
    print(f"k={k}: |2ReI0-2ReMV|/sc={mp.nstr(true_res/sc,4)}  exactWbound/sc={mp.nstr(W/sc,4)}  margin={mp.nstr(W/true_res,3)}")
