import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi1(x, N): return (-1+1j)*a + N*(1/x - bb/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x

def check(Nv):
    N = mp.mpf(Nv); th = mp.sqrt(pi*N)
    y = mp.sqrt(pi/(2*N))
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    x0 = 1/(rho*y)
    Phi0 = Phi(x0,N); p2 = Phi2(x0,N); g0 = g(x0)
    om = mp.im(Phi1(x0,N))
    c2 = abs(mp.re(p2))/2
    sg = 1/mp.sqrt(2*c2)
    w = mp.sqrt(mp.log(N)+1.0)
    d = min(w*sg, x0/3)
    xl, xr = x0-d, x0+d; x2 = x0/2
    P0 = abs(g0)*mp.e**mp.re(Phi0)
    def integrand(x): return mp.e**Phi(x,N)*g(x)
    def model(u): return g0*mp.e**Phi0*mp.e**(1j*om*u + p2*u*u/2)
    sc = N**mp.mpf('-1.5')*mp.e**(-th)
    # SIGNED window difference (oscillatory): split into subintervals for quad stability
    Js = 24
    W_signed = mp.mpc(0)
    for j in range(Js):
        u0, u1 = -d + 2*d*j/Js, -d + 2*d*(j+1)/Js
        W_signed += mp.quad(lambda u: integrand(x0+u)-model(u), [u0, u1])
    W_abs = mp.quad(lambda u: abs(integrand(x0+u)-model(u)), [0, d]) + mp.quad(lambda u: abs(integrand(x0+u)-model(u)), [-d, 0])
    # model tail exact via erfc: tail_R = e^{om^2/(2p2)} sqrt(-2/p2) (sqrtpi/2) erfc(sqrt(-p2/2)(d + i om/p2))
    rt = mp.sqrt(-p2/2)
    tailR = mp.e**(om**2/(2*p2))*mp.sqrt(-2/p2)*(sqpi/2)*mp.erfc(rt*(d + 1j*om/p2))
    tailL = mp.e**(om**2/(2*p2))*mp.sqrt(-2/p2)*(sqpi/2)*mp.erfc(rt*(d - 1j*om/p2))
    TM_exact = abs(g0*mp.e**Phi0*(tailR+tailL))
    # crude IBP estimate for sanity: 2 P0 e^{-c2 d^2}/(|p2| d - om)
    TM_crude = 2*P0*mp.e**(-c2*d*d)/(abs(p2)*d - om)
    TRt = mp.quad(lambda x: abs(integrand(x)), [xr, x0+8*d, mp.inf])
    TLmt = mp.quad(lambda x: abs(integrand(x)), [x2, xl])
    TLft = mp.quad(lambda x: abs(integrand(x)), [1, x2])
    print(f"N={Nv}: d/x0={mp.nstr(d/x0,5)}")
    print(f"  |signed window diff|/sc = {mp.nstr(abs(W_signed)/sc,10)}   (abs-integral/sc = {mp.nstr(W_abs/sc,10)})")
    print(f"  TM exact(erfc)/sc = {mp.nstr(TM_exact/sc,10)}   crude-IBP/sc = {mp.nstr(TM_crude/sc,10)}")
    print(f"  TR/sc = {mp.nstr(TRt/sc,10)}  TLm/sc = {mp.nstr(TLmt/sc,10)}  TLf/sc = {mp.nstr(TLft/sc,10)}")
for Nv in [10**4, 10**6, 10**8]:
    check(Nv)
