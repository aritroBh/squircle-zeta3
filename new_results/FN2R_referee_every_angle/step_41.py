import mpmath as mp, time, json
mp.mp.dps = 60
t0 = time.time()
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

beta = mp.pi**2/8; m = mp.mpf(2); cm = 1-1/m; ym = cm*(1-mp.log(cm))
f  = lambda x: mp.sin(mp.pi*x/2)
fp = lambda x: mp.pi/2*mp.cos(mp.pi*x/2)
def Dn_sin(th, n):
    nn = mp.mpf(n); c = mp.cos(th); lnn = mp.log(nn)
    def bump(v):
        tt = mp.exp(v)/nn
        if tt > 1: return mp.mpf(0)
        fx = f(1-tt)
        if fx <= 0: return mp.mpf(0)
        lyp = mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fp(1-tt))
        return hc(mp.exp(lyp), c)*mp.exp(v)/nn
    # peak at t ~ (2/pi) sqrt(log n / n): v ~ (ln n + ln ln n)/2 + const
    mid = (lnn + mp.log(lnn))/2
    pts = sorted(set([-60, mid-8, mid-2, mid+2, mid+8, lnn]))
    return mp.quad(bump, pts, maxdegree=10)

def Uroot(y):  # outer root of u - cm log u = y, m=2
    return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, -1)
def Vroot(y):
    return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, 0)
def Ffun(y): return Uroot(y)**(1/m)
def F2(y):
    U = Uroot(y); Up = U/(U-cm); Upp = -cm*Up/(U-cm)**2
    return (1/m)*((1/m-1)*U**(1/m-2)*Up**2 + U**(1/m-1)*Upp)
def ell(y):
    if y <= ym: return mp.mpf(0)
    return Uroot(y)**(1/m) - Vroot(y)**(1/m)

def full_model(th, n):
    c = mp.cos(th); Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn)
    # E[ell(Lam+xi)] = int_0^inf phi(s)/(1-c) ell(Lam + log(1/s)) ds ; s = e^{-u}
    def phi_s(s):
        r = mp.sqrt(1+2*c*s+s**2)
        return (1-c**2)/(r*(r+s+c))
    i1 = mp.quad(lambda u: phi_s(mp.exp(-u))*mp.exp(-u)*ell(Lam+u), [0, mp.inf])
    i2 = mp.quad(lambda u: phi_s(mp.exp(u))*mp.exp(u)*ell(Lam-u), [0, mp.inf])
    return (i1+i2)/Nn

res = {}
for thn, th in [('pi/3', mp.pi/3), ('2pi/3', 2*mp.pi/3), ('pi/2', mp.pi/2)]:
    c = mp.cos(th); k1,k2,k3 = kforms(th)
    for n in [10**4, 10**5, 10**6, 10**7]:
        D = Dn_sin(th, n)
        Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn)
        L0 = (1-c)*(mp.log(n)/(m*beta*n))**(1/m)
        L1 = (1-c)*Ffun(Lam)/Nn
        L3 = (1-c)*(Ffun(Lam+k1) + (k2/2)*F2(Lam+k1))/Nn
        Fm = full_model(th, n)
        predLam = 1 + ((m-1)*mp.log(Lam) + m*k1)/(m**2*Lam)
        predlogn = 1 + ((m-1)*mp.log(mp.log(n)) + m*k1 + mp.log(m*beta))/(m*mp.log(n))
        sc = (D/L0 - predLam)*(Lam/mp.log(Lam))**2
        sclogn = (D/L0 - predlogn)*(mp.log(n)/mp.log(mp.log(n)))**2
        res[(thn,n)] = [mp.nstr(x,30) for x in [D,L0,L1,L3,Fm]]
        print("th=%5s n=%8d: D/L0=%s predLam=%s predln=%s sc=%s scln=%s D/L3=%s D/Fm=%s" % (
            thn, n, mp.nstr(D/L0,7), mp.nstr(predLam,7), mp.nstr(predlogn,7),
            mp.nstr(sc,5), mp.nstr(sclogn,5), mp.nstr(D/L3,7), mp.nstr(D/Fm,9)))
json.dump({str(k): v for k,v in res.items()}, open("fn2r_contact.json","w"), indent=1)
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f2:
    f2.write(f"{time.strftime('%H:%M')} step_09: contact sin(pi x/2) m=2, 3 angles, n=1e4..1e7: D/L0 vs refined predictions, scaled residuals, D/L3, D/full-model [referee, NUMERICAL]\n")
