import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
b=mp.e**(1j*mp.pi/4)
def RePhi(t,N,a,j=0): 
    return -t/mp.sqrt(2)+(N-1)*mp.log(t)-(N-j-mp.mpf(3)/2)/2*mp.log(t**2+mp.sqrt(2)*a*t+a**2)
def RePhip(t,N,a,j=0):
    q=t**2+mp.sqrt(2)*a*t+a**2
    return -1/mp.sqrt(2)+(N-1)/t-(N-j-mp.mpf(3)/2)*(t+a/mp.sqrt(2))/q
def RePhipp(t,N,a,j=0):
    q=t**2+mp.sqrt(2)*a*t+a**2
    return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*(t**2+mp.sqrt(2)*a*t)/q**2
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
a=2*mp.pi/3
for N in [100,400,1600]:
    j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    # inflection: where RePhipp=0 for t>t0
    try:
        tc=mp.findroot(lambda t: RePhipp(t,N,a,j), (t0*1.5, t0*3))
    except: tc=None
    print(f"N={N}: t0={mp.nstr(t0,6)}, RePhip(t0)={mp.nstr(RePhip(t0,N,a,j),4)}")
    # check RePhip sign and RePhipp sign at various t>t0
    for mult in [1.2,1.5,2,3,5]:
        t=t0*mult
        print(f"   t={mp.nstr(mult,2)}*t0: RePhip={mp.nstr(RePhip(t,N,a,j),4)} RePhipp={mp.nstr(RePhipp(t,N,a,j),5)}")
    print(f"   inflection tc/t0={mp.nstr(tc/t0 if tc else None,4)}")
