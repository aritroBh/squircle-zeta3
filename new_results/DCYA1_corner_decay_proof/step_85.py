import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def RePhi(t,N,a,j=0): return mp.re(Phi(t,N,a,j))
def tail_parts(k,a):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    c2=-mp.re(Phipp(t0,N,a,j)); sig=1/mp.sqrt(2*c2); w=mp.sqrt(mp.log(N)+0.5); delta=min(w*sig,t0/3)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    ap=abs(pre); sc=N**(-1.5)*mp.e**(-mp.sqrt(2*a*N))
    # right tail (t>t0+delta), left-mid (t0/2,t0-delta), far-left (1,t0/2), near 0 (0,1)
    TR=ap*mp.quad(lambda t: mp.e**RePhi(t,N,a,j),[t0+delta,3*t0,10*t0,mp.inf])
    TLm=ap*mp.quad(lambda t: mp.e**RePhi(t,N,a,j),[t0/2,t0-delta])
    TLf=ap*mp.quad(lambda t: mp.e**RePhi(t,N,a,j),[mp.mpf(1),t0/2])
    TL0=ap*mp.quad(lambda t: mp.e**RePhi(t,N,a,j),[mp.mpf('0.001'),mp.mpf(1)])
    print(f"k={k}: TR/sc={mp.nstr(TR/sc,3)} TLm/sc={mp.nstr(TLm/sc,3)} TLf/sc={mp.nstr(TLf/sc,3)} TL0/sc={mp.nstr(TL0/sc,3)}")
a=2*mp.pi/3
for k in [100,400,1600,6400]:
    tail_parts(k,a)
