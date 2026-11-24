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
    P=[mp.mpc(1)]
    for k in range(1,n+1):
        P.append((om*P[k-1]+(k-1)*P[k-2] if k>=2 else om*P[k-1])/p2)
    return (-1j)**n*mp.sqrt(2*mp.pi/(-p2))*mp.e**(om**2/(2*p2))*P[n]
def bound(k,a,verbose=False):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j); Reps=mp.re(ps)
    c2=-mp.re(p2); sig=1/mp.sqrt(c2)
    delta=min(sig*mp.sqrt(2*mp.log(N)), t0/3)
    lo,hi=t0-delta,t0+delta
    def M(fn): return max(abs(fn(lo,N,a,j)),abs(fn(t0,N,a,j)),abs(fn(hi,N,a,j)))
    M3,M4,M5=M(Phippp),M(Phipppp),M(Phip5)
    p3,p4=Phippp(t0,N,a,j),Phipppp(t0,N,a,j)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    eReps=mp.e**Reps
    # model value
    MV=pre*mp.e**ps*gauss_moment(0,p2,om)*b
    # W1 first-order (full-line gauss moments)
    W1=abs(pre)*eReps*(abs(p3)/6*abs(gauss_moment(3,p2,om))+abs(p4)/24*abs(gauss_moment(4,p2,om)))
    # W3 second-order, gaussian-weighted, M's on window
    def S(u): return M3*abs(u)**3/6+M4*u**4/24+M5*abs(u)**5/120
    W3=0.5*abs(pre)*eReps*mp.e**max(S(-delta),S(delta))*mp.quad(lambda u: mp.e**(0.5*mp.re(p2)*u**2)*S(u)**2,[-delta,0,delta])
    # model tail
    Tmodel=2*abs(pre)*eReps*mp.e**(-0.5*c2*delta**2)/(c2*delta)
    # actual tail via RePhi concavity (numeric)
    Tact=abs(pre)*(mp.quad(lambda t: mp.e**mp.re(Phi(t,N,a,j)),[t0+delta,3*t0,10*t0,mp.inf])+mp.quad(lambda t: mp.e**mp.re(Phi(t,N,a,j)),[mp.mpf('0.3'),t0-delta]))
    return dict(MV=MV,W1=W1,W3=W3,Tmodel=Tmodel,Tact=Tact,delta=delta,t0=t0)
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [400,1600]:
    N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th); lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    bd=bound(k,a); I0=mainterm_I0(k,a)
    true_res=abs(2*mp.re(I0)-2*mp.re(bd['MV']))
    tot=bd['W1']+bd['W3']+bd['Tmodel']+bd['Tact']
    print(f"k={k}: |2ReI0-2ReMV|/sc={mp.nstr(true_res/sc,4)}")
    print(f"   W1/sc={mp.nstr(bd['W1']/sc,3)} W3/sc={mp.nstr(bd['W3']/sc,3)} Tm/sc={mp.nstr(bd['Tmodel']/sc,3)} Ta/sc={mp.nstr(bd['Tact']/sc,3)} tot/sc={mp.nstr(tot/sc,4)}")
    print(f"   |2ReMV-lead|/sc={mp.nstr(abs(2*mp.re(bd['MV'])-lead)/sc,3)}")
