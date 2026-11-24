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
def RePhi(t,N,a,j=0): return -t/mp.sqrt(2)+(N-1)*mp.log(t)-(N-j-mp.mpf(3)/2)/2*mp.log(t**2+mp.sqrt(2)*a*t+a**2)
def RePhip(t,N,a,j=0):
    q=t**2+mp.sqrt(2)*a*t+a**2; return -1/mp.sqrt(2)+(N-1)/t-(N-j-mp.mpf(3)/2)*(t+a/mp.sqrt(2))/q
def Phi3b(t,N,a):
    q=t**2+mp.sqrt(2)*a*t+a**2; return (t**3+6*a*(N-1)*t**2+6*a**2*(N-1)*t+2*a**3*(N-1))/(t**3*q**mp.mpf('1.5'))
def Phi4b(t,N,a):
    q=t**2+mp.sqrt(2)*a*t+a**2; return (3*t**4+24*a*(N-1)*t**3+36*a**2*(N-1)*t**2+24*a**3*(N-1)*t+6*a**4*(N-1))/(t**4*q**2)
def gmoment(n,p2,om):
    P=[mp.mpc(1)]
    for kk in range(1,n+1): P.append((om*P[kk-1]+(kk-1)*P[kk-2] if kk>=2 else om*P[kk-1])/p2)
    return (-1j)**n*mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))*P[n]
def mom(j,c): return mp.gamma((j+1)/2)/c**((j+1)/2)
def mainterm_I0(k,a):
    N=k+1; gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a); return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*mp.quad(f,[0,xs/8,xs/2,xs,2*xs,6*xs,mp.inf])
def full_bound(k,a,J=48):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    c2=-RePhipp(t0,N,a,j); sig=1/mp.sqrt(c2); delta=min(sig*mp.sqrt(mp.log(N)+1.0),t0/3)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    P0=abs(pre)*mp.e**mp.re(ps); MV0=pre*mp.e**ps*gmoment(0,p2,om)*b
    W1=P0*abs(om)*(Phi3b(t0,N,a)/6*mom(4,c2)); W2=mp.mpf(0); W3=mp.mpf(0)
    for side in [1,-1]:
        for jj in range(J):
            u0=delta*jj/J; u1=delta*(jj+1)/J; xsegL=t0-u1 if side==-1 else t0
            M3=Phi3b(xsegL,N,a); M4=Phi4b(xsegL,N,a)
            R4b=M4*u1**4/24; V=M3*u1**3/6; gw=mp.e**(-c2*u0**2/2); du=u1-u0
            W2+=P0*gw*R4b*du; W3+=P0*gw*(V*V/2)*du
    xr=t0+delta; xl=t0-delta
    TR=abs(pre)*mp.e**RePhi(xr,N,a,j)/abs(RePhip(xr,N,a,j))
    TLm=abs(pre)*mp.e**RePhi(xl,N,a,j)/RePhip(xl,N,a,j)
    sg=1/mp.sqrt(c2); x=delta/sig; TM=P0*2*sg*mp.e**(-x**2/2)/x
    return MV0,W1,W2,W3,TM,TR,TLm
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
print("N : |I0-MV|/sc : bound/sc : margin")
for N in [100,200,400,800,1600,3200,6400,12800,25600,10**5,10**6]:
    k=N-1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th)
    try: I0=mainterm_I0(k,a)
    except: I0=mp.nan
    MV0,W1,W2,W3,TM,TR,TLm=full_bound(k,a)
    res=abs(2*mp.re(I0)-2*mp.re(MV0)) if not mp.isnan(I0) else mp.nan
    tot=W1+W2+W3+TM+TR+TLm
    print(f"N={N}: {mp.nstr(res/sc,4) if not mp.isnan(I0) else 'NA'} : {mp.nstr(tot/sc,4)} : {mp.nstr(tot/res,3) if not mp.isnan(I0) else 'NA'}")
