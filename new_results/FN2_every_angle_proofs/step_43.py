import mpmath as mp, time, json
mp.mp.dps = 50
t0 = time.time()

def hc(s, c):
    if s <= 0: return mp.mpf(0)
    if s < 1: return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c) - (1-c**2)/(r + s + c)
def S12(z):
    return mp.zeta(3) - mp.polylog(3, 1-z) + mp.log(1-z)*mp.polylog(2, 1-z) \
           + mp.mpf('0.5')*mp.log(1-z)**2*mp.log(z)
def cumulants(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,z)-u**2, -12*w*S12(z)+2*(1+u**3))

# f = sin(pi x/2): Phi(tau) = -log cos(pi tau/2), beta = pi^2/8, m=2
beta = mp.pi**2/8; m = mp.mpf(2)
def f_sin(x): return mp.sin(mp.pi*x/2)
def fp_sin(x): return mp.pi/2*mp.cos(mp.pi*x/2)

def Dn_sin(th, n):
    nn = mp.mpf(n); c = mp.cos(th)
    Nn = mp.sqrt(beta*nn)
    mid = mp.log(nn/Nn)
    def integrand(v):
        tt = mp.exp(v)/nn
        if tt >= 1: return mp.mpf(0)
        fx = f_sin(1-tt)
        if fx <= 0: return mp.mpf(0)
        yp = nn*mp.exp((nn-1)*mp.log(fx))*fp_sin(1-tt)
        return hc(yp, c)*mp.exp(v)/nn
    lnn = mp.log(nn)
    pts = [-200, mid-12, mid-5, mid, mid+4, min(mid+9, lnn-mp.mpf('0.5')), lnn]
    return mp.quad(integrand, pts, maxdegree=10)

# Lambert-W model pieces (m=2): c_m = 1/2; U solves U - (1/2) log U = y
cm = 1 - 1/m
ym = cm*(1 - mp.log(cm))
def Uroot(y):
    if y < ym: return None
    return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, -1)
def Vroot(y):
    if y < ym: return None
    return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, 0)
def Ffun(y):  # U^{1/m}
    U = Uroot(y)
    return None if U is None else U**(1/m)
def ellfun(y):
    U = Uroot(y)
    if U is None: return mp.mpf(0)
    return U**(1/m) - Vroot(y)**(1/m)

def full_model(th, n):
    """(1-c)/N * E[ell(Lambda + xi_theta) ; > y_m]  via integration over s."""
    c = mp.cos(th); Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn)
    def phi_st(s):
        r = mp.sqrt(1 + 2*c*s + s**2)
        return (1-c**2)/(r*(r+s+c))
    # integrate over u = log s in (-inf, inf): ell(Lam - u) phi(e^u) e^u du
    def integrand(u):
        return ellfun(Lam - u)*phi_st(mp.exp(u))*mp.exp(u)
    I = mp.quad(integrand, [-60, Lam-12, Lam-4, Lam, Lam+4, 60], maxdegree=10)
    return I/Nn   # note: includes (1-c) since phi_st has mass (1-c)

def laws(th, n):
    c = mp.cos(th); k1,k2,k3 = cumulants(th)
    Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn)
    L0 = (1-c)*(mp.log(n)/(m*beta*n))**(1/m)
    ustar = Ffun(Lam)                       # = u*
    L1 = (1-c)*ustar/Nn
    L2 = (1-c)*Ffun(Lam+k1)/Nn
    # k2 correction: F'' via finite differences (or analytic); F'' = F/( (2U-1)^2 ) * (1-2U+... ) do analytic:
    # F = U^{1/2}; U' = U/(U-1/2); F' = U^{1/2}/(2(U-1/2)); F'' = d/dy ...
    def Fp(y):
        U = Uroot(y); return U**(1/m-1)/(m*(U-cm)) * U**(1/m)*0 + U**(1/m-1)*U/(m*(U-cm))*U**0  # placeholder
    # do it cleanly: F(y) = U(y)^{1/m}, U' = U/(U-cm)
    def F1(y):
        U = Uroot(y); Up = U/(U-cm); return (1/m)*U**(1/m-1)*Up
    def F2(y):
        U = Uroot(y); Up = U/(U-cm); Upp = Up**2*(cm)/(U*(U-cm))*0  # derive: U'' = d/dy[U/(U-cm)] = Up*( (U-cm)-U )/(U-cm)^2 = -cm Up/(U-cm)^2
        Upp = -cm*Up/(U-cm)**2
        return (1/m)*((1/m-1)*U**(1/m-2)*Up**2 + U**(1/m-1)*Upp)
    L3 = (1-c)*(Ffun(Lam+k1) + (k2/2)*F2(Lam+k1))/Nn
    return L0, L1, L2, L3, Lam, ustar

angles = {'pi/3': mp.pi/3, 'pi/2': mp.pi/2, '2pi/3': 2*mp.pi/3}
out = {}
for thn, th in angles.items():
    for n in [10**4, 10**5, 10**6, 10**7, 10**8]:
        D = Dn_sin(th, n)
        L0,L1,L2,L3,Lam,us = laws(th, n)
        Fm = full_model(th, n)
        out[(thn,n)] = dict(D=D, L0=L0,L1=L1,L2=L3 if False else L2, L3=L3, Fm=Fm)
        print(f"th={thn:6s} n={n:8d}: D/L0={mp.nstr(D/L0,7)} D/L1={mp.nstr(D/L1,7)} "
              f"D/L2={mp.nstr(D/L2,7)} D/L3={mp.nstr(D/L3,7)} D/Fm={mp.nstr(D/Fm,8)}")
json.dump({str(k): {kk: mp.nstr(vv,30) for kk,vv in v.items()} for k,v in out.items()},
          open('contact_results.json','w'), indent=1)
print("elapsed", time.time()-t0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} contact numerics: sin, 3 angles, n=1e4..1e8, D vs L0/L1/L2/L3/full-model  [NUMERICAL]\n")
