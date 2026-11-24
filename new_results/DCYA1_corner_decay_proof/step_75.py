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
def poly_exp(c,D):
    P=[mp.mpc(0)]*(D+1); P[0]=1
    for d in range(1,D+1): P[d]=mp.fsum(k*c[k]*P[d-k] for k in range(1,d+1) if k<len(c))/d
    return P
def mainterm_I0(k,a):
    N=k+1; gp=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))
    f=lambda x: mp.e**(-x)*x**(N-1)*(x+1j*a)**(-(N-mp.mpf(3)/2))
    xs=mp.sqrt(N*a); return gp*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)*mp.quad(f,[0,xs/4,xs,4*xs,mp.inf])

def saddle_full(k,a,K,D):
    N=k+1; j=0
    t0=mp.findroot(lambda t: mp.re(Phip(t,N,a,j)), mp.sqrt(N*a))
    om=mp.im(Phip(t0,N,a,j)); p2=Phipp(t0,N,a,j); ps=Phi(t0,N,a,j)
    c=[mp.mpc(0)]*3+[Phider(t0,N,a,n,j)/math.factorial(n) for n in range(3,K+1)]
    cK=c+[mp.mpc(0)]*(D-K)
    P=poly_exp(cK,D)
    pre=mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4-a/2))*(1j/mp.pi)*mp.gamma(mp.mpf(3)/2)*mp.gamma(N-mp.mpf(3)/2)/mp.gamma(N)
    S=mp.fsum(P[d]*gmoment(d,p2,om) for d in range(D+1))
    MVK=pre*mp.e**ps*S*b
    # R_tail bound: R(u)-R_poly(u), i.e. terms n>=K+1 of phase Taylor. bound by M_{K+1}/(K+1)! |u|^{K+1} * 1/(1-r) 
    c2=-mp.re(p2); sig=1/mp.sqrt(c2); delta=min(6*sig, t0/3)
    # M_{K+1}: max |Phi^{K+1}| on window; and assume geometric after
    lo=t0-delta
    MK1=max(abs(Phider(lo,N,a,K+1,j)),abs(Phider(t0,N,a,K+1,j)),abs(Phider(t0+delta,N,a,K+1,j)))
    # R_poly(u) max modulus on window
    def Rpoly(u): return mp.fsum(c[n]*u**n for n in range(3,K+1))
    # eps1 = max |e^{R_poly}-P| on window (sample)
    us=[-delta+i*2*delta/40 for i in range(41)]
    eps1=max(abs(mp.e**Rpoly(u)-mp.fsum(P[d]*u**d for d in range(D+1))) for u in us)
    # eps2 = max |e^{R_poly}| |R_tail| e^{|R_tail|} on window; R_tail ~ MK1/(K+1)! |u|^{K+1}/(1-q), q~ratio
    def Rtail_bound(u):
        return MK1/math.factorial(K+1)*abs(u)**(K+1)/(1-mp.mpf('0.5'))
    eps2=max(abs(mp.e**Rpoly(u))*Rtail_bound(u)*mp.e**Rtail_bound(u) for u in us)
    # remainder E0 = (eps1+eps2) * |pre| e^{Re ps} sqrt(2 pi/c2)  + tails
    E0win=(eps1+eps2)*abs(pre)*mp.e**mp.re(ps)*mp.sqrt(2*mp.pi/c2)
    # tails: actual modulus tail + model tail
    Tact=abs(pre)*(mp.quad(lambda t: mp.e**mp.re(Phi(t,N,a,j)),[t0+delta,3*t0,10*t0,mp.inf])+mp.quad(lambda t: mp.e**mp.re(Phi(t,N,a,j)),[mp.mpf('0.3'),t0-delta]))
    return MVK, E0win+Tact, dict(eps1=eps1,eps2=eps2,delta=delta,sig=sig)
a=2*mp.pi/3; A=mp.sqrt(2*a*mp.sin(a))
for k in [100,400]:
    N=k+1; th=mp.sqrt(2*a*N); sc=N**(-1.5)*mp.e**(-th); lead=-(A/N)*mp.e**(-th)*mp.cos(th)
    I0=mainterm_I0(k,a); MVK,E0,info=saddle_full(k,a,8,16)
    print(f"k={k}: |2ReI0-2ReMVK|/sc={mp.nstr(abs(2*mp.re(I0)-2*mp.re(MVK))/sc,4)}  E0/sc={mp.nstr(E0/sc,4)}  (eps1={mp.nstr(info['eps1'],3)},eps2={mp.nstr(info['eps2'],3)})")
