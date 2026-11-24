import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
print("N : t0/sqrt(N a) : RePhi(t0)+theta : -RePhi''(t0)*N : omega*sqrt(N) : (model/lead ratio)")
for N in [100,400,1600,6400,25600]:
    j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    th=mp.sqrt(2*a*N)
    # model value MV_0 (quadratic)
    def gmoment0(p2,om): return mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    MV0=pre*mp.e**ps*gmoment0(p2,om)*b
    lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    print(f"N={N}: {mp.nstr(t0/mp.sqrt(N*a),6)} : {mp.nstr(mp.re(ps)+th,6)} : {mp.nstr(-mp.re(p2)*N,6)} : {mp.nstr(om*mp.sqrt(N),6)} : {mp.nstr(2*mp.re(MV0)/lead,6)}")
