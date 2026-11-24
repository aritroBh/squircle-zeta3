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
def analyze(k,a,KK):
    """F1..F_{KK-1} via Gaussian moments of e^{R_poly} (degree-limited), plus modulus remainder of order KK."""
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    c2=-mp.re(p2); sig=1/mp.sqrt(c2); delta=min(6*sig,t0/3)
    c=[mp.mpc(0)]*3+[Phider(t0,N,a,n,j)/math.factorial(n) for n in range(3,KK+1)]  # R_poly to order KK
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    # model MV0 (quadratic)
    MV0=pre*mp.e**ps*gmoment(0,p2,om)*b
    # F via e^{R_poly} truncated to order KK in the *exponential* (sum R_poly^n/n! as polynomials, keep deg<=KK)
    # Build E_poly = sum_{n>=1} R_poly^n/n! truncated to polynomial degree KK (in u)
    R=c+[mp.mpc(0)]*(KK-len(c)+1)
    E=[mp.mpc(0)]*(KK+1)
    cur=[mp.mpc(0)]*(KK+1); cur[0]=1  # powers of R_poly
    import math as m
    for n in range(1,KK+1):
        # cur *= R_poly / n ... build R_poly^n
        new=[mp.mpc(0)]*(KK+1)
        for d in range(KK+1):
            for e in range(3,min(KK-d,len(R)-1)+1):
                if e<len(R): new[d+e]+=cur[d]*R[e]
        cur=new
        for d in range(KK+1): E[d]+=cur[d]/m.factorial(n)
    # F1 = pre e^{ps} sum_{d>=1} E[d] moment(d) b  (all orders up to deg KK) -- this is the full polynomial model to order KK
    F1=pre*mp.e**ps*mp.fsum(E[d]*gmoment(d,p2,om) for d in range(1,KK+1))*b
    # modulus remainder of order KK+1: |R|^{KK+1}/(KK+1)! e^{|R|} gaussian-weighted (numeric)
    MKK=max(abs(Phider(t0-delta,N,a,KK+1,j)),abs(Phider(t0,N,a,KK+1,j)),abs(Phider(t0+delta,N,a,KK+1,j)))
    def Rfull(u): return Phi(t0+u,N,a,j)-ps-0.5*p2*u**2-1j*om*u
    W3=0.5*abs(pre)*mp.quad(lambda u: mp.e**mp.re(ps+0.5*p2*u**2)*abs(Rfull(u))**2*mp.e**abs(Rfull(u)),[-delta,0,delta])
    return MV0,F1,W3
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [100,400]:
    I0=mainterm_I0(k,a); N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th); lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    for KK in [4,8,12]:
        MV0,F1,W3=analyze(k,a,KK)
        r=2*mp.re(I0)-2*mp.re(MV0)
        print(f"k={k} KK={KK}: |I0-MV0|/sc={mp.nstr(abs(r)/sc,3)} |I0-MV0-F1|/sc={mp.nstr(abs(r-2*mp.re(F1))/sc,3)} |F1|/sc={mp.nstr(abs(2*mp.re(F1))/sc,3)} W3mod/sc={mp.nstr(W3/sc,3)}")
