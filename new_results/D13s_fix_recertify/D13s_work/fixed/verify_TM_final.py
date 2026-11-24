
import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi1(x, N): return (-1+1j)*a + N*(1/x - bb/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x

# certified pieces from fixed final cell (arb): TMreg, TRreg, TLmreg
TMreg = mp.mpf('0.1223723258162209663165629878873688709761810860711856105338671421354680')
TRreg = mp.mpf('0.280825611949804834267966285156637300776651197108579577869439550913293')
TLmreg = mp.mpf('1.535002056440645689329067312875911123804586144112063536292603212900568')
cv = (2/pi)*(pi/2)**mp.mpf('1.5')
print(f"certified TM C-part = {mp.nstr(cv*TMreg,12)}  TR = {mp.nstr(cv*TRreg,12)}  TLm = {mp.nstr(cv*TLmreg,12)}")

for Nv in [10**8, 10**10, 10**12, 10**14]:
    N = mp.mpf(Nv); th = mp.sqrt(pi*N)
    y = mp.sqrt(pi/(2*N))
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    x0 = 1/(rho*y); t=rho*y
    Phi0 = Phi(x0,N); p2 = Phi2(x0,N); g0 = g(x0)
    om = mp.im(Phi1(x0,N))
    c2 = abs(mp.re(p2))/2
    w = mp.sqrt(mp.log(N)+1.0)
    d = min(w/mp.sqrt(2*c2), x0/3)
    sc = N**mp.mpf('-1.5')*mp.e**(-th)
    rt = mp.sqrt(-p2/2)
    tailR = mp.e**(om**2/(2*p2))*mp.sqrt(-2/p2)*(sqpi/2)*mp.erfc(rt*(d + 1j*om/p2))
    tailL = mp.e**(om**2/(2*p2))*mp.sqrt(-2/p2)*(sqpi/2)*mp.erfc(rt*(d - 1j*om/p2))
    TM_exact = abs(g0*mp.e**Phi0*(tailR+tailL))
    # TR and TLm truths: direct integration
    x2 = x0/2; xl = x0-d; xr = x0+d
    TRt = mp.quad(lambda x: abs(mp.e**Phi(x,N)*g(x)), [xr, x0+6*d, mp.inf])
    TLmt = mp.quad(lambda x: abs(mp.e**Phi(x,N)*g(x)), [x2, xl])
    print(f"N={Nv:g}: TM truth/sc = {mp.nstr(TM_exact/sc,12)} (bound {mp.nstr(cv*TMreg,10)}: margin {mp.nstr(cv*TMreg/(TM_exact/sc),6)})")
    print(f"          TR truth/sc = {mp.nstr(TRt/sc,12)} (bound {mp.nstr(cv*TRreg,10)}: margin {mp.nstr(cv*TRreg/(TRt/sc),6)})")
    print(f"          TLm truth/sc = {mp.nstr(TLmt/sc,12)} (bound {mp.nstr(cv*TLmreg,10)}: margin {mp.nstr(cv*TLmreg/(TLmt/sc),6)})")
