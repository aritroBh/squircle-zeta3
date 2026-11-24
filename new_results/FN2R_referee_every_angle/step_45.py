import mpmath as mp, time
mp.mp.dps = 60
t0=time.time()
def hc(s, c):
    if s <= 0: return mp.mpf(0)
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c) - (1-c**2)/(r + s + c)
def S12(zv):
    return mp.zeta(3) - mp.polylog(3, 1-zv) + mp.log(1-zv)*mp.polylog(2, 1-zv) \
           + mp.mpf('0.5')*mp.log(1-zv)**2*mp.log(zv)
def kforms(th):
    zv = mp.sin(th/2)**2; w2 = 1/mp.tan(th/2)**2; u = w2*mp.log(mp.cos(th/2)**2)
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,zv)-u**2, 2+2*u**3-12*w2*S12(zv))
def phi_s(s,c):
    r = mp.sqrt(1+2*c*s+s**2)
    return (1-c**2)/(r*(r+s+c))

def setup(mv, betav):
    m = mp.mpf(mv); beta = mp.mpf(betav); cm = 1-1/m
    if m > 1:
        ym = cm*(1-mp.log(cm))
        U = lambda y: -cm*mp.lambertw(-mp.exp(-y/cm)/cm, -1)
        V = lambda y: -cm*mp.lambertw(-mp.exp(-y/cm)/cm, 0)
        def ell(y):
            if y <= ym: return mp.mpf(0)
            return (U(y)**(1/m) - V(y)**(1/m)).real
    else:
        ym = -mp.inf
        U = lambda y: abs(cm)*mp.lambertw(mp.exp(y/abs(cm))/abs(cm), 0)
        ell = lambda y: (U(y)**(1/m)).real
    F = lambda y: U(y)**(1/m)
    def F2(y):
        Uy = U(y); Up = Uy/(Uy-cm); Upp = -cm*Up/(Uy-cm)**2
        return (1/m)*((1/m-1)*Uy**(1/m-2)*Up**2 + Uy**(1/m-1)*Upp)
    return m, beta, cm, ym, U, ell, F, F2

def Dn_exact(th, n, m, beta):
    """D_n = int h_c(y') dx for f = exp(-beta (1-x)^m); v = log(n t)."""
    nn = mp.mpf(n); c = mp.cos(th); lnn = mp.log(nn)
    def bump(v):
        tt = mp.exp(v)/nn
        if tt > 1: return mp.mpf(0)
        lyp = mp.log(nn) + mp.log(beta*m*tt**(m-1)) - nn*beta*tt**m
        return hc(mp.exp(lyp), c)*mp.exp(v)/nn
    # window center: t ~ (log n/(beta n))^{1/m}: v ~ log n + (1/m)(log log n - log(beta n)) = (1-1/m)lnn + (1/m)log(lnn/beta)
    mid = (1-1/m)*lnn + mp.log(lnn/beta)/m
    return mp.quad(bump, sorted(set([-80, mid-10, mid-3, mid+3, mid+10, lnn])), maxdegree=10)

def Fm_model(th, n, m, beta, ell, ym):
    c = mp.cos(th); Nn = (beta*mp.mpf(n))**(1/m); Lam = mp.log(m*Nn)
    if ym == -mp.inf:
        pts = [-80, 0, Lam-4, Lam, Lam+4, 60]
    else:
        us = Lam - ym
        pts = sorted(set([-80, 0, us-4, us, us+2, us+8, 60]))
    I = mp.quad(lambda u: ell(Lam-u)*phi_s(mp.exp(u),c)*mp.exp(u), pts, maxdegree=10)
    return I/Nn

for mv, betav in [(mp.mpf('1.5'), 1), (mp.mpf('0.5'), 1)]:
    m, beta, cm, ym, U, ell, F, F2 = setup(mv, betav)
    print("="*20, "m =", mv, " beta =", betav, " cm =", mp.nstr(cm,6), " ym =", mp.nstr(ym if ym>-mp.inf else mp.mpf('-9e9'),6))
    for thn, th in [('pi/3', mp.pi/3), ('pi/2', mp.pi/2)]:
        c = mp.cos(th); k1,k2,k3 = kforms(th)
        for n in [10**4, 10**6]:
            D = Dn_exact(th, n, m, beta)
            Nn = (beta*mp.mpf(n))**(1/m); Lam = mp.log(m*Nn)
            L0 = (1-c)*(mp.log(n)/(m*beta*n))**(1/m)
            L3 = (1-c)*(F(Lam+k1) + (k2/2)*F2(Lam+k1))/Nn
            Fm = Fm_model(th, n, m, beta, ell, ym)
            predlogn = 1 + ((m-1)*mp.log(mp.log(n)) + m*k1 + mp.log(m*beta))/(m*mp.log(n))
            sclogn = (D/L0 - predlogn)*(mp.log(n)/mp.log(mp.log(n)))**2
            print(" th=%5s n=%8d: D=%s D/Fm=%s D/L0=%s predln=%s scln=%s D/L3=%s" % (
                thn,n,mp.nstr(D,6),mp.nstr(D/Fm,10),mp.nstr(D/L0,7),mp.nstr(predlogn,7),mp.nstr(sclogn,5),mp.nstr(D/L3,7)))
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f2:
    f2.write(f"{time.strftime('%H:%M')} step_10: contact regimes m=3/2 (E[S^{1/(m-1)}]=infty gap regime) and m=1/2 (infinite slope), exact Phi models, theta=pi/3,pi/2: full model and refined law tested [referee, NUMERICAL]\n")
