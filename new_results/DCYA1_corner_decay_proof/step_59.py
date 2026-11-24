import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
def Phip(t,N,a):
    b=mp.e**(1j*mp.pi/4); return -b+(N-1)/t-(N-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a):
    b=mp.e**(1j*mp.pi/4); return -(N-1)/t**2+(N-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def Phi(t,N,a):
    b=mp.e**(1j*mp.pi/4); return -b*t+(N-1)*mp.log(b*t)-(N-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def I0_model(k,a):
    N=k+1; b=mp.e**(1j*mp.pi/4)
    ts=mp.findroot(lambda t: Phip(t,N,a), mp.sqrt(N*a))
    gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    gauss=mp.sqrt(2*mp.pi/(-Phipp(ts,N,a)))
    J_model=mp.e**Phi(ts,N,a)*gauss*b  # *b from dx=b dt
    return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*J_model
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [400,1600,6400]:
    N=k+1; th=mp.sqrt(2*a*N); lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    mv=2*mp.re(I0_model(k,a))
    print(f"k={k}: model 2ReI0={mp.nstr(mv,10)}  lead={mp.nstr(lead,10)}  ratio={mp.nstr(mv/lead,7)}")
