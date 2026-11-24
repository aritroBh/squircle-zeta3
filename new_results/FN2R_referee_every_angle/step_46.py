import mpmath as mp, time
mp.mp.dps = 60
t0=time.time()
def hc(s, c):
    if s <= 0: return mp.mpf(0)
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c) - (1-c**2)/(r + s + c)
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
    return m, beta, cm, ym, U, ell
def Dn_exact(th, n, m, beta):
    nn = mp.mpf(n); c = mp.cos(th); lnn = mp.log(nn)
    def bump(v):
        tt = mp.exp(v)/nn
        if tt > 1: return mp.mpf(0)
        lyp = mp.log(nn) + mp.log(beta*m*tt**(m-1)) - nn*beta*tt**m
        return hc(mp.exp(lyp), c)*mp.exp(v)/nn
    mid = (1-1/m)*lnn + mp.log(lnn/beta)/m
    return mp.quad(bump, sorted(set([-80, mid-10, mid-3, mid+3, mid+10, lnn])), maxdegree=10)
def Fm_model(th, n, m, beta, ell, ym):
    c = mp.cos(th); Nn = (beta*mp.mpf(n))**(1/m); Lam = mp.log(m*Nn)
    if ym == -mp.inf: pts = [-80, 0, Lam-4, Lam, Lam+4, 60]
    else:
        us = Lam - ym
        pts = sorted(set([-80, 0, us-4, us, us+2, us+8, 60]))
    I = mp.quad(lambda u: ell(Lam-u)*phi_s(mp.exp(u),c)*mp.exp(u), pts, maxdegree=10)
    return I/Nn
th = mp.pi/3
for mv in [mp.mpf('1.5'), mp.mpf('0.5'), mp.mpf(2), mp.mpf(3)]:
    m, beta, cm, ym, U, ell = setup(mv, 1)
    for n in [10**4, 10**5, 10**6]:
        D = Dn_exact(th, n, m, beta); Fm = Fm_model(th, n, m, beta, ell, ym)
        print("m=%s n=%8d: 1-D/Fm = %s" % (mp.nstr(m,4), n, mp.nstr(1-D/Fm,3)))
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f2:
    f2.write(f"{time.strftime('%H:%M')} step_10b: exact-Phi models: |1-D/Fm| at theta=pi/3 for m=1.5,0.5,2,3 quantified (see result.md table); Lambert-W reduction (a) essentially exact [referee]\n")
