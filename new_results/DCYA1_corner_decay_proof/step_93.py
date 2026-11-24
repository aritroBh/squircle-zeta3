import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
b=mp.e**(1j*mp.pi/4)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def RePhipp(t,N,a,j=0):
    q=t**2+mp.sqrt(2)*a*t+a**2
    return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*(t**2+mp.sqrt(2)*a*t)/q**2
a=2*mp.pi/3
for N in [100,400,1600]:
    j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    c2_0=-RePhipp(t0,N,a,j)
    sig=1/mp.sqrt(c2_0)
    print(f"N={N}: t0={mp.nstr(t0,6)} c2(t0)={mp.nstr(c2_0,5)} sig={mp.nstr(sig,4)}")
    for mult in [0.5,0.7,1.0,1.3,1.5,2.0,3.0]:
        t=t0*mult
        c=-RePhipp(t,N,a,j)
        print(f"   t={mult}*t0: c2={mp.nstr(c,5)}  ratio to c2(t0)={mp.nstr(c/c2_0,4)}")
