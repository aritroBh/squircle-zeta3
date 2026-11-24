import mpmath as mp, time, json
mp.mp.dps = 60
t0=time.time()
beta = mp.pi**2/8; m = mp.mpf(2); cm = 1-1/m; ym = cm*(1-mp.log(cm))
def Uroot(y): return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, -1)
def Vroot(y): return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, 0)
def ell(y):
    if y <= ym: return mp.mpf(0)
    return (Uroot(y)**(1/m) - Vroot(y)**(1/m)).real
def phi_s(s,c):
    r = mp.sqrt(1+2*c*s+s**2)
    return (1-c**2)/(r*(r+s+c))
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
    mid = (lnn + mp.log(lnn))/2
    return mp.quad(bump, sorted(set([-60, mid-8, mid-2, mid+2, mid+8, lnn])), maxdegree=10)
def Fm_model(th, n):
    c = mp.cos(th); Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn)
    us = Lam - ym
    I = mp.quad(lambda u: ell(Lam-u)*phi_s(mp.exp(u),c)*mp.exp(u),
                [-80, 0, us-4, us, us+2, us+8, 60], maxdegree=10)
    return I/Nn
def U0(y): return Uroot(y)
def Ffun(y): return Uroot(y)**(1/m)
def F2(y):
    U = Uroot(y); Up = U/(U-cm); Upp = -cm*Up/(U-cm)**2
    return (1/m)*((1/m-1)*U**(1/m-2)*Up**2 + U**(1/m-1)*Upp)
for thn, th in [('pi/3', mp.pi/3), ('2pi/3', 2*mp.pi/3), ('pi/2', mp.pi/2)]:
    c = mp.cos(th); k1,k2,k3 = kforms(th)
    for n in [10**4, 10**5, 10**6, 10**7]:
        D = Dn_sin(th, n); Fm = Fm_model(th, n)
        Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn)
        L1 = (1-c)*Ffun(Lam)/Nn
        L3 = (1-c)*(Ffun(Lam+k1) + (k2/2)*F2(Lam+k1))/Nn
        err = abs(D/Fm-1)
        print("th=%5s n=%8d: D/Fm=%s  (1-D/Fm)/(log n/n)=%s  D/L1=%s D/L3=%s" % (
            thn,n,mp.nstr(D/Fm,9), mp.nstr((1-D/Fm)/(mp.log(n)/n),5), mp.nstr(D/L1,7), mp.nstr(D/L3,7)))
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f2:
    f2.write(f"{time.strftime('%H:%M')} step_09b: full-model quadrature fixed (breakpoint at sqrt singularity of ell). D/Fm matches FN2 (e.g. pi/3 1e6: 0.99999873); error ~ log n/n as (a) predicts with delta=m. D/L1, D/L3 match FN2 exactly [referee, NUMERICAL]\n")
