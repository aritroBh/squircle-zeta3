# Cross-check (numerical, mpmath 60dps): the certified W1,W2,W3 from step_34 must
# upper-bound the EXACT pieces they were designed to bound:
#  W1t = |P0 int e^{p2 u^2/2} (Ph3c u^3/6 + gap u) i sin(om u) du|   (odd part)
#  W2t = |P0 int e^{p2 u^2/2} (R4 + ga2) e^{i om u} du|             (even 2nd order)
#  W3t = |P0 int e^{p2 u^2/2} e^{i om u} (e^{R+ga}-1-(R+ga)) du|    (exp remainder)
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, gamma, findroot, quad, sin, cos
mp.dps = 60
b = exp(mpc(0,pi/4)); a = pi/(2*sqrt(2)); rt2 = sqrt(2); sqpi = sqrt(pi)
def H(s): return sqpi/4*gamma(s/2)/gamma(s/2+mpf("1.5"))
def Phi(x,N):
    s=b*x; return mpc(0,1)*pi*s/2 + N*log(s/(1+s))
def RePhi(x,N): return -a*x + N*log(x) - N*log(1+rt2*x+x*x)/2
def ReP1(x,N):  return -a + N/x - N*(x+1/rt2)/(1+rt2*x+x*x)
def ReP2(x,N):  return -N/x**2 + N*x*(x+rt2)/(1+rt2*x+x*x)**2
def Phi2(x,N):
    s=b*x; return -N*(1+2*s)/(x*x*(1+s)**2)
def Phi3(x,N):
    s=b*x; return 2*N/x**3 - 2*N*b**3/(1+s)**3
def g(x): return H(b*x)/x
def gap(x):
    s=b*x
    return b*(mp.polygamma(0,s/2)-mp.polygamma(0,s/2+mpf("1.5")))/2 - 1/x
def x0_of(N): return findroot(lambda x: ReP1(x,N), sqrt(2*N/pi))

def pieces(Nv):
    N = mpf(Nv); th = sqrt(pi*N)
    x0 = x0_of(N); y0 = 1/x0; ys = sqrt(pi/(2*N))
    om = (N/rt2)*(ys**2 - y0**2/(1+rt2*y0+y0**2))
    c2 = abs(ReP2(x0,N))/2
    sig = 1/sqrt(2*c2); w = sqrt(log(N)+1)
    d = min(w*sig, x0/3)
    Ph0 = Phi(x0,N); g0 = g(x0); P0 = abs(g0)*exp(RePhi(x0,N))
    p2 = Phi2(x0,N); Ph3c = Phi3(x0,N); gapc = gap(x0); lg0 = log(g0)
    def R(u):  return Phi(x0+u,N) - Ph0 - mpc(0,1)*om*u - p2*u*u/2
    def R4(u): return R(u) - Ph3c*u**3/6
    def ga(u): return log(g(x0+u)) - lg0
    def ga2(u): return ga(u) - gapc*u
    sc = N**mpf("-1.5")*exp(-th)
    f1 = lambda u: abs(exp(p2*u*u/2))*(Ph3c*u**3/6 + gapc*u)*mpc(0,1)*sin(om*u)
    W1t = abs(P0*quad(f1, [-d,0,d]))
    f2 = lambda u: exp(p2*u*u/2)*(R4(u)+ga2(u))*exp(mpc(0,1)*om*u)
    W2t = abs(P0*quad(f2, [-d,0,d]))
    f3 = lambda u: exp(p2*u*u/2)*exp(mpc(0,1)*om*u)*(exp(R(u)+ga(u))-1-(R(u)+ga(u)))
    W3t = abs(P0*quad(f3, [-d,0,d]))
    # full I - MV truth
    MV = g0*exp(Ph0)*sqrt(-2*pi/p2)*exp(om**2/(2*p2))
    If = lambda u: exp(Phi(x0+u,N)-Ph0)*g(x0+u)/g0
    Iwin = P0*quad(If, [-d, 0, d])
    print(f"N={Nv}: d={mp.nstr(d,6)}")
    print(f"  W1 truth/sc = {mp.nstr(W1t/sc,6)}   (cert W1/sc = 10.7454/11.9079)")
    print(f"  W2 truth/sc = {mp.nstr(W2t/sc,6)}   (cert W2/sc = 1.8813/4.4679)")
    print(f"  W3 truth/sc = {mp.nstr(W3t/sc,6)}   (cert W3/sc = 11.0078/28.0595)")
    print(f"  W1+W2+W3 truth/sc = {mp.nstr((W1t+W2t+W3t)/sc,6)} vs cert sum = {mp.nstr(mpf('10.7454')+mpf('1.8813')+mpf('11.0078') if Nv==100 else mpf('11.9079')+mpf('4.4679')+mpf('28.0595'),6)}")
for Nv in [100, 400]:
    pieces(Nv)
