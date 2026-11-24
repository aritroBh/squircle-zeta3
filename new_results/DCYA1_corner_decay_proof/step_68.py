import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
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
def analyze(k,a,order2=False):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    MV=pre*mp.e**ps*gmoment(0,p2,om)*b
    NTR=10
    cn=[Phider(t0,N,a,n,j)/math.factorial(n) for n in range(NTR+1)]  # R(u)=sum cn[n] u^n (n>=3)
    F1=sum(cn[n]*gmoment(n,p2,om) for n in range(3,NTR+1))
    out=pre*mp.e**ps*F1*b
    if order2:
        # R^2 coefficients
        R2=[mp.mpc(0)]*(2*NTR+1)
        for i in range(3,NTR+1):
            for l in range(3,NTR+1):
                R2[i+l]+=cn[i]*cn[l]
        F2=sum(R2[n]*gmoment(n,p2,om) for n in range(6,2*NTR+1))
        out=pre*mp.e**ps*(F1+0.5*F2)*b
    return MV,out
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [400,1600]:
    N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th)
    I0=mainterm_I0(k,a); MV,F1=analyze(k,a,False); _,F2=analyze(k,a,True)
    r=2*mp.re(I0)-2*mp.re(MV)
    print(f"k={k}: |I0-MV|/sc={mp.nstr(abs(r)/sc,3)} |I0-MV-F1|/sc={mp.nstr(abs(r-2*mp.re(F1))/sc,3)} |I0-MV-F2|/sc={mp.nstr(abs(r-2*mp.re(F2))/sc,4)}")
