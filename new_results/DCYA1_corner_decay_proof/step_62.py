import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0): return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0): return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0): return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def Phippp(t,N,a,j=0): return 2*(N-1)/t**3-2*(N-j-mp.mpf(3)/2)*b**3/(b*t+1j*a)**3
def Phipppp(t,N,a,j=0): return -6*(N-1)/t**4+6*(N-j-mp.mpf(3)/2)*b**4/(b*t+1j*a)**4
def Phip5(t,N,a,j=0): return 24*(N-1)/t**5-24*(N-j-mp.mpf(3)/2)*b**5/(b*t+1j*a)**5
def mainterm_I0(k,a):
    N=k+1; gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a); J=mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
    return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*J

def gauss_moment(n, p2, om):
    # int_RR u^n e^{.5 p2 u^2 + i om u} du  (Re p2<0)
    Z=mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))
    # moments via derivatives w.r.t (i om): M_n = d^n Z/d(i om)^n
    # use recurrence: represent Z as function of w=i*om; M_n=Z^{(n)}(w)
    # do by series: Z(w)=Z0 e^{w^2/(2 p2)} * ... actually Z depends on om via om^2/(2 p2).
    # Let s=om; M_n = d^n/d(is)^n [sqrt(2pi/(-p2)) e^{s^2/(2 p2)}]
    # Use explicit: int u^n e^{.5 p2 u^2 + i s u}du. Set c=-p2>0 complex... use hermite-like
    # Compute by numerical formula via Faa: just use sympy-free recurrence with bell polynomials
    # Simpler: M_n = (-i)^n d^n/ds^n Z(s), Z(s)=Zc e^{s^2/(2 p2)}, Zc=sqrt(2pi/(-p2))
    Zc=mp.sqrt(2*mp.pi/(-p2))
    # Taylor of e^{s^2/(2p2)} around s=om up to order n, get derivative
    # d^n/ds^n e^{s^2/(2p2)} = e^{s^2/(2p2)} * H_n(s/p2) / p2^{n}... use recursion
    # P_0=1, P_{k+1}= (s P_k + k P_{k-1})/p2  (for d^n e^{s^2/(2p2)}/ds^n = e^{...} P_n(s))
    P=[mp.mpc(1)]
    for k in range(1,n+1):
        P.append((om*P[k-1]+(k-1)*P[k-2] if k>=2 else om*P[k-1])/p2)
    return (-1j)**n*Zc*mp.e**(om**2/(2*p2))*P[n]

def window_bound(k,a):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j); Reps=mp.re(ps)
    c2=-mp.re(p2); sig=1/mp.sqrt(c2); delta=sig*mp.sqrt(2*mp.log(N))
    lo,hi=t0-delta,t0+delta
    M3=max(abs(Phippp(x,N,a,j)) for x in [lo,t0,hi])
    M4=max(abs(Phipppp(x,N,a,j)) for x in [lo,t0,hi])
    M5=max(abs(Phip5(x,N,a,j)) for x in [lo,t0,hi])
    # Model value
    MV=mp.e**ps*gauss_moment(0,p2,om)*b
    gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    pre=gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    MV_full=pre*MV
    # W1: first-order remainder |pre * int e^{Phi_m} R_{34} du|, R_{34}=p3 u^3/6 + p4 u^4/24
    p3=Phippp(t0,N,a,j); p4=Phipppp(t0,N,a,j)
    R34=abs(pre)*(abs(p3)/6*abs(gauss_moment(3,p2,om))+abs(p4)/24*abs(gauss_moment(4,p2,om)))
    # W3: second-order: .5 e^{Reps} e^{Rmax} int e^{.5 c2 u^2... } (M3|u|^3/6+M4 u^4/24+M5|u|^5/120)^2 * e^{.5 Re p2 u^2}du
    Rmax=M3*delta**3/6+M4*delta**4/24+M5*delta**5/120
    def w3int(u):
        Rmag=M3*abs(u)**3/6+M4*u**4/24+M5*abs(u)**5/120
        return mp.e**(0.5*mp.re(p2)*u**2)*Rmag**2
    W3=0.5*abs(pre)*mp.e**Reps*mp.e**Rmax*mp.quad(w3int,[-delta,0,delta])
    # model tail: |pre| int_{|u|>delta} |e^{Phi_m}| du
    Tmodel=2*abs(pre)*mp.e**Reps*mp.e**(-0.5*c2*delta**2)/(c2*delta)
    # actual tail: bound Re Phi beyond delta via concavity (numeric estimate of e^{RePhi} integral)
    def RePhi(t): return mp.re(Phi(t,N,a,j))
    # find Re Phi at delta points, approximate tail integral assuming exp decay from delta value
    # crude: integrate e^{RePhi} on [t0+delta, 3*t0] and [t0/2, t0-delta] numerically + far tail
    Tact=abs(pre)*(mp.quad(lambda t: mp.e**RePhi(t),[t0+delta, 3*t0])+mp.quad(lambda t: mp.e**RePhi(t),[mp.mpf('0.5'),t0-delta]))
    return dict(MV_full=MV_full,W1=R34,W3=W3,Tmodel=Tmodel,Tact=Tact,delta=delta,sig=sig,om=om,Rmax=Rmax,pre=pre)

a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [400]:
    N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th); lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    wb=window_bound(k,a)
    I0=mainterm_I0(k,a)
    tot=wb['W1']+wb['W3']+wb['Tmodel']+wb['Tact']
    print(f"k={k}: |2ReI0-2ReMV|/sc={mp.nstr(abs(2*mp.re(I0)-2*mp.re(wb['MV_full']))/sc,4)}")
    print(f"  W1/sc={mp.nstr(wb['W1']/sc,3)} W3/sc={mp.nstr(wb['W3']/sc,3)} Tmodel/sc={mp.nstr(wb['Tmodel']/sc,3)} Tact/sc={mp.nstr(wb['Tact']/sc,3)}")
    print(f"  total bound/sc={mp.nstr(tot/sc,4)}  (should exceed |2ReI0-2ReMV|/sc={mp.nstr(abs(2*mp.re(I0)-2*mp.re(wb['MV_full']))/sc,3)})")
    print(f"  |2ReMV-lead|/sc={mp.nstr(abs(2*mp.re(wb['MV_full'])-lead)/sc,3)}  Rmax={mp.nstr(wb['Rmax'],3)} delta/sig={mp.nstr(wb['delta']/wb['sig'],3)}")
