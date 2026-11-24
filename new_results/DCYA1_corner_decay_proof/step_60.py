import sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
b=mp.e**(1j*mp.pi/4)
def Phi(t,N,a,j=0):
    return -b*t+(N-1)*mp.log(b*t)-(N-j-mp.mpf(3)/2)*mp.log(b*t+1j*a)
def Phip(t,N,a,j=0):
    return -b+(N-1)/t-(N-j-mp.mpf(3)/2)*b/(b*t+1j*a)
def Phipp(t,N,a,j=0):
    return -(N-1)/t**2+(N-j-mp.mpf(3)/2)*b**2/(b*t+1j*a)**2
def Phippp(t,N,a,j=0):
    return 2*(N-1)/t**3-2*(N-j-mp.mpf(3)/2)*b**3/(b*t+1j*a)**3
def Phipppp(t,N,a,j=0):
    return -6*(N-1)/t**4+6*(N-j-mp.mpf(3)/2)*b**4/(b*t+1j*a)**4

def saddle(N,a,j=0):
    ts=mp.findroot(lambda t: Phip(t,N,a,j), mp.sqrt(N*a))
    return ts

def mainterm_true(k,a):
    # 2 Re I_+^{(0)} by direct quad
    N=k+1; gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a)
    J=mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])
    I0=gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*J
    return 2*mp.re(I0)

def saddle_bound(k,a):
    """Compute a candidate E(N) bounding |2 Re I_+^{(0)} - MV|, MV=model value."""
    N=k+1; j=0
    ts=saddle(N,a,j)
    ps=Phi(ts,N,a,j); p2=Phipp(ts,N,a,j)
    c2=-mp.re(p2); sig=1/mp.sqrt(c2)
    delta=sig*mp.sqrt(2*mp.log(N))
    # M3, M4 on window (bound by endpoint evals, they decrease away from ts; take max)
    lo,hi=ts-delta,ts+delta
    M3=max(abs(Phippp(lo,N,a,j)),abs(Phippp(ts,N,a,j)),abs(Phippp(hi,N,a,j)))
    M4=max(abs(Phipppp(lo,N,a,j)),abs(Phipppp(ts,N,a,j)),abs(Phipppp(hi,N,a,j)))
    Reps=mp.re(ps)
    # model value
    MV=mp.e**ps*mp.sqrt(2*mp.pi/(-p2))*b
    # window remainder: int_W e^{Re ps + .5 Re p2 u^2} |e^R-1| du, |R|<=M3|u|^3/6+M4 u^4/24
    def integrand(u):
        Rmag=M3*abs(u)**3/6+M4*u**4/24
        return mp.e**(Reps+0.5*mp.re(p2)*u**2)*Rmag*mp.e**Rmag
    W_in=mp.quad(integrand,[-delta,0,delta])
    # model tail (Gaussian): 2 * int_delta^inf e^{Re ps - .5 c2 u^2} du <= 2 e^{Reps} e^{-.5 c2 delta^2}/(c2 delta)
    T_model=2*mp.e**Reps*mp.e**(-0.5*c2*delta**2)/(c2*delta)
    # actual tail: bound Re Phi beyond delta. Use secant: at delta, RePhi decays; bound integral
    # numerically compute the actual tail integral of e^{Re Phi} for the bound design
    return dict(ts=ts,ps=ps,p2=p2,c2=c2,sig=sig,delta=delta,M3=M3,M4=M4,MV=MV,
                W_in=W_in,T_model=T_model,Reps=Reps)
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [400]:
    sb=saddle_bound(k,a)
    true=mainterm_true(k,a); N=k+1; th=mp.sqrt(2*a*N); lead=-(A/N)*mp.e**(-th)*mp.cos(th); sc=N**(-1.5)*mp.e**(-th)
    print(f"k={k}: true 2ReI0={mp.nstr(true,8)}")
    print(f"  MV model (2 Re) = {mp.nstr(2*mp.re(sb['MV']),8)}")
    print(f"  |true-MV|/sc = {mp.nstr(abs(true-2*mp.re(sb['MV']))/sc,4)}")
    print(f"  W_in/sc={mp.nstr(sb['W_in']/sc,4)}  T_model/sc={mp.nstr(sb['T_model']/sc,4)}")
    print(f"  M3 delta^3/6={mp.nstr(sb['M3']*sb['delta']**3/6,4)} M4 delta^4/24={mp.nstr(sb['M4']*sb['delta']**4/24,4)}")
