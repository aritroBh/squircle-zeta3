import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
for a in [2*mp.pi/3, mp.pi/3]:
    for N in [400,1600,6400]:
        j=0
        t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
        om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
        c2=-mp.re(p2); sig=1/mp.sqrt(c2)
        print(f"a={mp.nstr(a,4)} N={N}: sig={mp.nstr(sig,3)}")
        for m in [1,2,3,4,6]:
            u=m*sig
            R=Phi(t0+u,N,a,j)-ps-0.5*p2*u**2-1j*om*u
            print(f"   u={m}sig: Re R={mp.nstr(mp.re(R),4)} Im R={mp.nstr(mp.im(R),4)}  |e^R|=e^{mp.nstr(mp.re(R),3)}")
