
import mpmath as mp
mp.mp.dps = 50
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi1(x, N): return (-1+1j)*a + N*(1/x - bb/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x

def truth(Nv):
    N = mp.mpf(Nv); th = mp.sqrt(pi*N)
    y = mp.sqrt(pi/(2*N))
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-45'))
    x0 = 1/(rho*y)
    Phi0 = Phi(x0,N); p2 = Phi2(x0,N); g0 = g(x0)
    om = mp.im(Phi1(x0,N))
    c2 = abs(mp.re(p2))/2
    w = mp.sqrt(mp.log(N)+1.0)
    d = min(w/mp.sqrt(2*c2), x0/3)
    xl, xr = x0-d, x0+d; x2 = x0/2
    sc = N**mp.mpf('-1.5')*mp.e**(-th)
    integrand = lambda x: mp.e**Phi(x,N)*g(x)
    model = lambda u: g0*mp.e**Phi0*mp.e**(1j*om*u + p2*u*u/2)
    MV = g0*mp.e**Phi0*mp.sqrt(-2*pi/p2)*mp.e**(om**2/(2*p2))
    Wt = mp.quad(lambda u: abs(integrand(x0+u)-model(u)), [0, d]) + mp.quad(lambda u: abs(integrand(x0+u)-model(u)), [-d, 0])
    Min = mp.quad(lambda u: model(u), [-d, d])
    TMt = abs(MV - Min)
    TRt = mp.quad(lambda x: abs(integrand(x)), [xr, x0+6*d, mp.inf])
    TLmt = mp.quad(lambda x: abs(integrand(x)), [x2, xl])
    TLft = mp.quad(lambda x: abs(integrand(x)), [1, x2])
    TL0t = mp.quad(lambda x: abs(integrand(x)), [0, mp.mpf('0.5'), 1])
    Lt = -1j*(pi**mp.mpf('1.5')/(2*N))*mp.e**(-th)*mp.e**(1j*th)
    Dt = (2/pi)*abs(MV-Lt)
    lhs = abs((2/pi)*mp.im(mp.quad(integrand, [0, 1, x2, xl, x0, xr, mp.inf])) - (2/pi)*mp.im(Lt))
    tot = (Wt+TMt+TRt+TLmt+TLft+TL0t)*(2/pi)/sc + Dt/sc
    print(f"N={Nv}: |Mmain-lead|/sc = {mp.nstr(lhs/sc,10)}  ; piecewise-sum/sc = {mp.nstr(tot,10)}  (d/x0={mp.nstr(d/x0,5)})")
for Nv in [1891, 10000, 100000]:
    truth(Nv)
