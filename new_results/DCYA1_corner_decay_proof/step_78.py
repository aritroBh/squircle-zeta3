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
def mainterm_I0(k,a):
    N=k+1; gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a); return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
def bounds(k,a,K=12):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    c2=-mp.re(p2); sig=1/mp.sqrt(c2); delta=min(6*sig,t0/3)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    MV0=pre*mp.e**ps*gmoment(0,p2,om)*b
    cn=[Phider(t0,N,a,n,j)/math.factorial(n) for n in range(3,K+1)]
    # W1 = |sum c_n moment(n)|  (first-order value)
    W1=abs(pre)*abs(mp.e**ps)*abs(mp.fsum(cn[i]*gmoment(3+i,p2,om) for i in range(len(cn))))*abs(b)
    # W3 = 1/2 int e^{RePhi_m} |R_full|^2 e^{|R_full|} du
    def Rfull(u): return Phi(t0+u,N,a,j)-ps-0.5*p2*u**2-1j*om*u
    W3=0.5*abs(pre)*mp.quad(lambda u: mp.e**mp.re(ps+0.5*p2*u**2)*abs(Rfull(u))**2*mp.e**abs(Rfull(u)),[-delta,0,delta])
    return MV0,W1,W3
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [100,200,400,800,1600]:
    N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th)
    I0=mainterm_I0(k,a); MV0,W1,W3=bounds(k,a)
    res=abs(2*mp.re(I0)-2*mp.re(MV0))
    print(f"k={k}: |I0-MV0|/sc={mp.nstr(res/sc,4)}  W1/sc={mp.nstr(W1/sc,4)} W3/sc={mp.nstr(W3/sc,4)}  W1+W3={mp.nstr((W1+W3)/sc,4)}  margin={mp.nstr((W1+W3)/res,3)}")
