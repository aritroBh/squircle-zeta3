import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=40
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def Phi3b(t,N,a): return 2*(N-1)/t**3 + 2*(N-mp.mpf('1.5'))/(t**2+mp.sqrt(2)*a*t+a**2)**mp.mpf('1.5')
def Phi4b(t,N,a): return 6*(N-1)/t**4 + 6*(N-mp.mpf('1.5'))/(t**2+mp.sqrt(2)*a*t+a**2)**2
def mom(j,c): return mp.gamma((j+1)/2)/c**((j+1)/2)
def window_bound(k,a,J=40,cw=mp.mpf('0.5')):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    c2=-mp.re(p2); sig=1/mp.sqrt(2*c2); w=mp.sqrt(mp.log(N)+cw)
    delta=min(w*sig, t0/3)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    P0=abs(pre)*mp.e**mp.re(ps)
    W1=P0*abs(om)*(Phi3b(t0,N,a)/6*mom(4,c2))
    W2=mp.mpf(0); W3=mp.mpf(0)
    for side in [1,-1]:
        for jj in range(J):
            u0=delta*jj/J; u1=delta*(jj+1)/J
            xsegL=t0-u1 if side==-1 else t0
            M3=Phi3b(xsegL,N,a); M4=Phi4b(xsegL,N,a)
            R4b=M4*u1**4/24; V=M3*u1**3/6
            gw=mp.e**(-c2*u0**2); du=u1-u0
            W2+=P0*gw*R4b*du; W3+=P0*gw*(V*V/2)*du
    return W1,W2,W3
a=2*mp.pi/3
for k in [100,400,1600,6400,25600,102400]:
    N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th)
    W1,W2,W3=window_bound(k,a)
    print(f"k={k}: W1/sc={mp.nstr(W1/sc,3)} W2/sc={mp.nstr(W2/sc,3)} W3/sc={mp.nstr(W3/sc,4)}")
