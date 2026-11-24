import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = pi/(2*sqrt2) if (pi:=mp.pi) else None
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
random.seed(45)
mn = mp.mpf('9'); mx = mp.mpf('-9')
for _ in range(300):
    y = mp.mpf(10)**random.uniform(-8,-0.4)
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    rho = 1+y*y*random.random()
    sig = (rho-1)/y**2
    x0 = 1/(rho*y); t = rho*y
    phi_true = mp.im(Phi(x0,N)) - th
    zz = sqrt2*t/(2+sqrt2*t); z2=zz*zz
    Ap = 1 - mp.atan(zz)/zz
    breg_h = -sig/(2*rho) + (-2*y*sig + 2*sqrt2*sig + sqrt2*rho - 2*sqrt2*rho*Ap)/(2*sqrt2*(2+sqrt2*t))
    Btrue = (phi_true + pi/4)/((pi/sqrt2)*y)
    d = Btrue - breg_h + mp.mpf(1)/12
    mn = min(mn, d); mx = max(mx, d)
print("Btrue - breg_h + 1/12 over 300 random (arbitrary rho): min =", mp.nstr(mn,8), " max =", mp.nstr(mx,8))

sqpi = mp.sqrt(pi)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
wz = mp.mpf('0')
for _ in range(200):
    y = mp.mpf(10)**random.uniform(-8,-0.4)
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    x0 = 1/(rho*y)
    Phi0 = Phi(x0,N); p2 = Phi2(x0,N)
    om = mp.im((-1+1j)*a + N*(1/x0 - bb/(1+bb*x0)))
    g0 = Hc(bb*x0)/x0
    MV = g0*mp.e**Phi0*mp.sqrt(-2*pi/p2)*mp.e**(om**2/(2*p2))
    Lt = -1j*(pi**mp.mpf('1.5')/(2*N))*mp.e**(-th)*mp.e**(1j*th)
    z_true = MV/Lt
    z_fact = (2*N*1j/pi**mp.mpf('1.5'))*g0*mp.e**(mp.re(Phi0)+th)*mp.sqrt(-2*pi/p2)*mp.e**(om**2/(2*p2))*mp.e**(1j*(mp.im(Phi0)-th))
    wz = max(wz, abs(z_true-z_fact)/abs(z_true))
print("z_true = factored form: worst rel residual =", mp.nstr(wz,10))
N = mp.mpf('400'); th = mp.sqrt(pi*N)
Lt = -1j*(pi**mp.mpf('1.5')/(2*N))*mp.e**(-th)*mp.e**(1j*th)
print("leading check: (2/pi)Im Lt =", mp.nstr((2/pi)*mp.im(Lt),12),
      " ; -(sqpi/N)e^{-th}cos th =", mp.nstr(-sqpi/N*mp.e**(-th)*mp.cos(th),12))
