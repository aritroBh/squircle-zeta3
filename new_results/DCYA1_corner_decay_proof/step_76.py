import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def Phider(t,N,a,n,j=0):
    s=(-1)**n*math.factorial(n-1); return -s*(N-1)/t**n + s*(N-j-mp.mpf(3)/2)*b**n/(b*t+1j*a)**n
def gmoment(n,p2,om):
    P=[mp.mpc(1)]
    for kk in range(1,n+1): P.append((om*P[kk-1]+(kk-1)*P[kk-2] if kk>=2 else om*P[kk-1])/p2)
    return (-1j)**n*mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))*P[n]
def poly_exp(c,D):
    P=[mp.mpc(0)]*(D+1); P[0]=1
    for d in range(1,D+1): P[d]=mp.fsum(k*c[k]*P[d-k] for k in range(1,d+1) if k<len(c))/d
    return P
def mainterm_I0(k,a):
    N=k+1; gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a); return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
def MVK(k,a,K,D):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    c=[mp.mpc(0)]*3+[Phider(t0,N,a,n,j)/math.factorial(n) for n in range(3,K+1)]
    cK=c+[mp.mpc(0)]*(D-K); P=poly_exp(cK,D)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    S=mp.fsum(P[d]*gmoment(d,p2,om) for d in range(D+1))
    return pre*mp.e**ps*S*b
a=2*mp.pi/3
for k in [100,400]:
    I0=mainterm_I0(k,a); N=k+1; sc=N**(-1.5)*mp.e**(-mp.sqrt(2*a*N))
    for K,D in [(8,16),(12,24),(16,32)]:
        mv=MVK(k,a,K,D)
        print(f"k={k} K={K} D={D}: |2ReI0-2ReMVK|/sc={mp.nstr(abs(2*mp.re(I0)-2*mp.re(mv))/sc,4)}")
