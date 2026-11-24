import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x

def analyze(N):
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
    # code phireg box: Bh decreasing in A1, A1m = [0, z2/(1+z2)]
    zz = sqrt2*t/(2+sqrt2*t); z2 = zz*zz
    def Bh(A1):
        T2n = 2*y*y*sig - sqrt2*rho*y - 2*rho*A1*y
        return -y*sig/(2*rho) + T2n/(2*y*(2+sqrt2*t))
    Bh_lo, Bh_hi = Bh(z2/(1+z2)), Bh(0)
    phi_box = ((pi/sqrt2)*Bh_lo, (pi/sqrt2)*Bh_hi)
    A1_code_needed = (1 - mp.atan(zz)/zz)/y
    in_box = phi_box[0] <= phi_true <= phi_box[1]
    # |z-1| over the box arc vs true: z = z0 e^{i phi}, z0 = z_true e^{-i phi_true}
    z0 = z_true*mp.e**(-1j*phi_true)
    def mag(phi): return abs(z0*mp.e**(1j*phi) - 1)
    # max of mag over [phi_lo, phi_hi]: endpoints (mag is monotonic in |phi| near here; check both + interior min at 0 if contained)
    mbox = max(mag(phi_box[0]), mag(phi_box[1]))
    mtrue = abs(z_true-1)
    D_true = sqpi*mtrue*mp.sqrt(N)
    D_box = sqpi*mbox*mp.sqrt(N)
    print(f"N={N}: phi_true={mp.nstr(phi_true,8)}  phi_box=[{mp.nstr(phi_box[0],8)},{mp.nstr(phi_box[1],8)}]  contains: {in_box}")
    print(f"   A1_code={mp.nstr(A1_code_needed,8)} vs A1_up={mp.nstr(z2/(1+z2),8)};  |z_true-1|={mp.nstr(mtrue,8)}  max|z_box-1|={mp.nstr(mbox,8)}")
    print(f"   D_true={mp.nstr(D_true,8)}  D_boxupper={mp.nstr(D_box,8)}  D_box >= D_true: {mbox >= mtrue}")
for Nv in [100, 157, 400, 1000, 10000, 10**6, 10**10]:
    analyze(mp.mpf(Nv))
