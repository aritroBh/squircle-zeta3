import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def RePhipp(t,N,a,j=0):
    q=t**2+mp.sqrt(2)*a*t+a**2; return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*(t**2+mp.sqrt(2)*a*t)/q**2
def gmoment(n,p2,om):
    P=[mp.mpc(1)]
    for kk in range(1,n+1): P.append((om*P[kk-1]+(kk-1)*P[kk-2] if kk>=2 else om*P[kk-1])/p2)
    return (-1j)**n*mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))*P[n]
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
def MV0(k,a):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    return pre*mp.e**ps*gmoment(0,p2,om)*b
def I1(k,a):
    # g_1 integral modulus: |I_+^{(1)}|
    N=k+1; E=mp.e**(1j*a); gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    g1=-mp.mpf(3)/4-1j*E/(4*mp.sin(a))
    xs=mp.sqrt(N*a)
    J=mp.quad(lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(5)/2)),[0,xs/8,xs/2,xs,2*xs,6*xs,mp.inf])
    I1v=gp*g1*(1j/mp.pi)*mp.gamma(mp.mpf(5)/2)*mp.gamma(N-mp.mpf(5)/2)/mp.gamma(N)*J
    return abs(I1v)
print("N : |2Re MV0 - lead|/sc : |2Re I_+^{(1)}|/sc")
for N in [100,200,400,800,1600,3200]:
    k=N-1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th); lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    mv=MV0(k,a); d=abs(2*mp.re(mv)-lead); i1=2*I1(k,a)
    print(f"N={N}: {mp.nstr(d/sc,4)} : {mp.nstr(i1/sc,4)}")
