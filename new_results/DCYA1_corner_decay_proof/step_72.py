import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def Phider(t,N,a,n,j=0):
    s=(-1)**n*math.factorial(n-1)
    return -s*(N-1)/t**n + s*(N-j-mp.mpf(3)/2)*b**n/(b*t+1j*a)**n
def mainterm_I0(k,a):
    N=k+1; gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a); J=mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
    return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*J
def gmoment(n,p2,om):
    P=[mp.mpc(1)]
    for kk in range(1,n+1): P.append((om*P[kk-1]+(kk-1)*P[kk-2] if kk>=2 else om*P[kk-1])/p2)
    return (-1j)**n*mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))*P[n]

def poly_exp(c, D):
    """P(u)=exp(sum_{n} c[n] u^n) truncated to degree D. c indexed from 0."""
    P=[mp.mpc(0)]*(D+1); P[0]=1
    for d in range(1,D+1):
        P[d]=mp.fsum(k*c[k]*P[d-k] for k in range(1,d+1) if k<len(c))/d
    return P

def model_and_saddle(k,a,K=8,D=16):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    # R(u)=Phi(t0+u)-ps - (.5 p2 u^2 + i om u) = sum_{n>=3} c_n u^n
    c=[mp.mpc(0)]*3+[Phider(t0,N,a,n,j)/math.factorial(n) for n in range(3,K+1)]
    c+= [mp.mpc(0)]*(D-K)  # pad
    P=poly_exp(c,D)  # P(u) ~ e^{R_poly(u)} to degree D
    # model value: prefactor * e^{ps} * sum p_d moment(d) * b
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    S=mp.fsum(P[d]*gmoment(d,p2,om) for d in range(D+1))
    MVK=pre*mp.e**ps*S*b
    return MVK, dict(t0=t0,om=om,p2=p2,ps=ps,c=c,P=P,pre=pre)
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [100,400]:
    N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th); lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    I0=mainterm_I0(k,a); MVK,sd=model_and_saddle(k,a)
    print(f"k={k}: |2ReI0-2ReMVK|/sc={mp.nstr(abs(2*mp.re(I0)-2*mp.re(MVK))/sc,4)}  |2ReMVK-lead|/sc={mp.nstr(abs(2*mp.re(MVK)-lead)/sc,4)}  |2ReI0-lead|/sc={mp.nstr(abs(2*mp.re(I0)-lead)/sc,4)}")
