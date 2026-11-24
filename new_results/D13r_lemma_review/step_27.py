import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x

def get(N):
    y = mp.sqrt(pi/(2*N)); th = pi/(sqrt2*y)
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    sig = (rho-1)/y**2; x0 = 1/(rho*y); t = rho*y
    Phi0 = Phi(x0,N); p2 = Phi2(x0,N)
    om = mp.im((-1+1j)*a + N*(1/x0 - bb/(1+bb*x0)))
    MV = g(x0)*mp.e**Phi0*mp.sqrt(-2*pi/p2)*mp.e**(om**2/(2*p2))
    Lt = -1j*(pi**mp.mpf('1.5')/(2*N))*mp.e**(-th)*mp.e**(1j*th)
    z_true = MV/Lt
    phi_true = mp.im(Phi0) - th
    argz0 = mp.arg(z_true*mp.e**(-1j*phi_true))   # arg of z0 = all other factors
    return y, rho, sig, t, phi_true, argz0, z_true

print(" N        phi_true     arg(z0)      phi*=-arg(z0)   arg(z_true)   safe?(phi*<(phi_true+phi_hi)/2)")
for Nv in [100, 150, 225, 340, 400, 700, 2000, 1e4, 1e6, 1e8, 1e10, 1e12]:
    N = mp.mpf(Nv)
    y, rho, sig, t, phi_true, argz0, z_true = get(N)
    zz = sqrt2*t/(2+sqrt2*t); z2 = zz*zz
    def Bh(A1):
        T2n = 2*y*y*sig - sqrt2*rho*y - 2*rho*A1*y
        return -y*sig/(2*rho) + T2n/(2*y*(2+sqrt2*t))
    phi_lo = (pi/sqrt2)*Bh(z2/(1+z2)); phi_hi = (pi/sqrt2)*Bh(0)
    phistar = -argz0
    # condition for box max |z-1| >= |z_true-1|: NOT(phistar > (phi_true+phi_hi)/2 and phi_lo closer)
    cond = (phistar > (phi_true+phi_hi)/2)
    print(f"{Nv:>9g} {mp.nstr(phi_true,8):>13} {mp.nstr(argz0,8):>12} {mp.nstr(phistar,8):>13} {mp.nstr(mp.arg(z_true),10):>12}   DANGER={cond}")
# quantitatively: compute |z(phi)-1| at phi_true vs max over box, with z0 exact
def mag_at(phi, z0): return abs(z0*mp.e**(1j*phi)-1)
worst_ratio = mp.mpf('0'); warg=None
for k in range(200):
    Nv = 100*mp.mpf(1.05)**k
    if Nv > 1e14: break
    y, rho, sig, t, phi_true, argz0, z_true = get(Nv)
    z0 = abs(z_true)*mp.e**(1j*(mp.arg(z_true)-phi_true))
    zz = sqrt2*t/(2+sqrt2*t); z2 = zz*zz
    def Bh(A1):
        T2n = 2*y*y*sig - sqrt2*rho*y - 2*rho*A1*y
        return -y*sig/(2*rho) + T2n/(2*y*(2+sqrt2*t))
    phi_lo = (pi/sqrt2)*Bh(z2/(1+z2)); phi_hi = (pi/sqrt2)*Bh(0)
    mbox = max(mag_at(phi_lo,z0), mag_at(phi_hi,z0))
    mtrue = abs(z_true-1)
    r = mbox/mtrue
    if r < worst_ratio or worst_ratio==0: pass
    if r < 1: print("UNDERESTIMATE at N=", mp.nstr(Nv,6), "ratio box/true =", mp.nstr(r,8))
    worst_ratio = r if r<1 else worst_ratio
print("scan complete: any ratio<1 printed above (none => box dominates everywhere, but by luck of phase alignment)")
