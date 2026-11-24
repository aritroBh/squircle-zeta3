import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
# Phi(t) on ray x=e^{i pi/4} t for J0:
# Phi(t) = -e^{i pi/4} t + (N-1) Log(e^{i pi/4} t) - (N-3/2) Log(e^{i pi/4} t + i a)
# plus outer d/dx-> e^{i pi/4}. We analyze Re/Im Phi''', Phi'''' at saddle.
def Phi(t,N,a):
    return -mp.e**(1j*mp.pi/4)*t+(N-1)*mp.log(mp.e**(1j*mp.pi/4)*t)-(N-mp.mpf(3)/2)*mp.log(mp.e**(1j*mp.pi/4)*t+1j*a)
def Phip(t,N,a):
    b=mp.e**(1j*mp.pi/4)
    return -b+(N-1)/t-(N-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a):
    b=mp.e**(1j*mp.pi/4)
    return -(N-1)/t**2+(N-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def Phippp(t,N,a):
    b=mp.e**(1j*mp.pi/4)
    return 2*(N-1)/t**3-2*(N-mp.mpf(3)/2)*b**3/(b*t+1j*a)**3
def Phipppp(t,N,a):
    b=mp.e**(1j*mp.pi/4)
    return -6*(N-1)/t**4+6*(N-mp.mpf(3)/2)*b**4/(b*t+1j*a)**4
for a in [mp.pi/3, 2*mp.pi/3, mp.pi/2]:
    N=400
    ts=mp.findroot(lambda t: Phip(t,N,a), mp.sqrt(N*a))
    p2=Phipp(ts,N,a); p3=Phippp(ts,N,a); p4=Phipppp(ts,N,a)
    print(f"a={mp.nstr(a,5)}: ts={mp.nstr(mp.re(ts),7)}")
    print(f"   RePhi''={mp.nstr(mp.re(p2),5)} ImPhi''={mp.nstr(mp.im(p2),5)}")
    print(f"   RePhi'''={mp.nstr(mp.re(p3),6)} ImPhi'''={mp.nstr(mp.im(p3),6)}  |Phi'''|*sigma^3={mp.nstr(abs(p3)*(1/mp.sqrt(-mp.re(p2)))**3,5)}")
    print(f"   |Phi''''|={mp.nstr(abs(p4),6)}  sigma={mp.nstr(1/mp.sqrt(-mp.re(p2)),5)}")
