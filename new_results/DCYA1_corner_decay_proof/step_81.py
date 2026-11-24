import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def gmoment(n,p2,om):
    P=[mp.mpc(1)]
    for kk in range(1,n+1): P.append((om*P[kk-1]+(kk-1)*P[kk-2] if kk>=2 else om*P[kk-1])/p2)
    return (-1j)**n*mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))*P[n]
def Phi3b(t,N,a): return 2*(N-1)/t**3 + 2*(N-mp.mpf('1.5'))/(t**2+mp.sqrt(2)*a*t+a**2)**mp.mpf('1.5')
def Phi4b(t,N,a): return 6*(N-1)/t**4 + 6*(N-mp.mpf('1.5'))/(t**2+mp.sqrt(2)*a*t+a**2)**2
def mom(j,c): return mp.gamma((j+1)/2)/c**((j+1)/2)
def mainterm_I0(k,a):
    N=k+1; gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a); return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
def window_bound(k,a,J=40,cw=mp.mpf('0.5')):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    c2=-mp.re(p2); sig=1/mp.sqrt(2*c2); w=mp.sqrt(mp.log(N)+cw)
    delta=min(w*sig, t0/3)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    P0=abs(pre)*mp.e**mp.re(ps)
    MV0=pre*mp.e**ps*gmoment(0,p2,om)*b
    # W1 (odd cubic, omega trick)
    W1=P0*abs(om)*(Phi3b(t0,N,a)/6*mom(4,c2))
    W2=mp.mpf(0); W3=mp.mpf(0)
    for side in [1,-1]:
        for jj in range(J):
            u0=delta*jj/J; u1=delta*(jj+1)/J
            xsegL=t0-u1 if side==-1 else t0  # extreme (smallest t) point on this subinterval
            if side==1: xsegL=t0  # on right, min t is t0
            M3=Phi3b(xsegL,N,a); M4=Phi4b(xsegL,N,a)
            R4b=M4*u1**4/24; Rb=M3*u1**3/6; V=Rb
            gw=mp.e**(-c2*u0**2); du=u1-u0
            q=mp.e**max(mp.mpf(0),mp.mpf('0'))  # curvature variation ~0 for now
            W2+=P0*gw*R4b*du
            W3+=P0*gw*(V*V/2)*q*du
    # model tail (TM) and actual tail (Tact)
    p2mag=abs(p2); omh=abs(om); dlo=delta
    den=p2mag*dlo-omh
    ecd=mp.e**(-c2*dlo**2)
    TM=P0*2*(ecd/den + p2mag*ecd/(den**2*2*c2*dlo)) if den>0 else mp.mpf('1e99')
    Tact=abs(pre)*(mp.quad(lambda t: mp.e**mp.re(Phi(t,N,a,j)),[t0+delta,3*t0,10*t0,mp.inf])+mp.quad(lambda t: mp.e**mp.re(Phi(t,N,a,j)),[mp.mpf('0.3'),t0-delta]))
    return MV0,W1,W2,W3,TM,Tact
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [100,400,1600]:
    N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th)
    I0=mainterm_I0(k,a); MV0,W1,W2,W3,TM,Tact=window_bound(k,a)
    res=abs(2*mp.re(I0)-2*mp.re(MV0)); tot=W1+W2+W3+TM+Tact
    print(f"k={k}: |I0-MV0|/sc={mp.nstr(res/sc,4)} | W1={mp.nstr(W1/sc,3)} W2={mp.nstr(W2/sc,3)} W3={mp.nstr(W3/sc,3)} TM={mp.nstr(TM/sc,3)} Ta={mp.nstr(Tact/sc,3)} tot={mp.nstr(tot/sc,4)} margin={mp.nstr(tot/res,3)}")
