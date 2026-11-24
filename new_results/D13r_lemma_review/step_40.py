import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi1(x, N): return (-1+1j)*a + N*(1/x - bb/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x

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
    TM_crude = 2*abs(g0)*mp.e**mp.re(Phi0)*mp.e**(-c2*d*d)/(abs(p2)*d - om)
    print(f"N={Nv:g}: TM_exact/sc = {mp.nstr(TM_exact/sc,12)}   TM_crude/sc = {mp.nstr(TM_crude/sc,12)}   d/x0={mp.nstr(d/x0,6)}")
print("certified final-cell bound for TM_true/sc (from TMreg=0.10290245): (2/pi)(pi/2)^1.5*TMreg =", 
      mp.nstr((2/pi)*(pi/2)**mp.mpf('1.5')*mp.mpf('0.1029024501051197981475520730109980944078443916754430317318924352606686'),12))
