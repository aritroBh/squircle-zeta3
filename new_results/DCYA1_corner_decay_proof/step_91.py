import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def RePhi(t,N,a,j=0): return -t/mp.sqrt(2)+(N-1)*mp.log(t)-(N-j-mp.mpf(3)/2)/2*mp.log(t**2+mp.sqrt(2)*a*t+a**2)
def RePhip(t,N,a,j=0):
    q=t**2+mp.sqrt(2)*a*t+a**2; return -1/mp.sqrt(2)+(N-1)/t-(N-j-mp.mpf(3)/2)*(t+a/mp.sqrt(2))/q
a=2*mp.pi/3
for N in [100,400,1600,6400,25600,102400,10**6]:
    j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    c2=-mp.re(Phipp(t0,N,a,j)); sig=1/mp.sqrt(c2); delta=min(sig*mp.sqrt(2*mp.log(N)),t0/3)
    xr=t0+delta
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    ap=abs(pre); th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th)
    ps=RePhi(t0,N,a,j); pxr=RePhi(xr,N,a,j); rp=RePhip(xr,N,a,j)
    dec=ps-pxr  # decrease over window
    TRtangent=ap*mp.e**pxr/abs(rp)
    TRquad=ap*mp.quad(lambda t: mp.e**RePhi(t,N,a,j),[xr,3*t0,10*t0,mp.inf])
    print(f"N={N}: dec={mp.nstr(dec,4)} (vs c2d2/2={mp.nstr(c2*delta**2/2,4)}) e^dec/N^(1/4)={mp.nstr(mp.e**(-dec)/N**(-0.25),4)}")
    print(f"     RePhi'(xr)={mp.nstr(rp,4)} TRtang/sc={mp.nstr(TRtangent/sc,4)} TRquad/sc={mp.nstr(TRquad/sc,4)}")
