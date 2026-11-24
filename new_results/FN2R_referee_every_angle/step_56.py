import mpmath as mp, time
mp.mp.dps = 60
def hc(s, c_):
    if s <= 0: return mp.mpf(0)
    r = mp.sqrt(1 + 2*c_*s + s**2)
    return (1-c_) - (1-c_**2)/(r + s + c_)
def S12(zv):
    return mp.zeta(3) - mp.polylog(3, 1-zv) + mp.log(1-zv)*mp.polylog(2, 1-zv) \
           + mp.mpf('0.5')*mp.log(1-zv)**2*mp.log(zv)
def kforms(thv):
    zv = mp.sin(thv/2)**2; w2 = 1/mp.tan(thv/2)**2; u = w2*mp.log(mp.cos(thv/2)**2)
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,zv)-u**2, 2+2*u**3-12*w2*S12(zv))
def Jpaper(thv, a_, A_, B_, n):
    k1,k2,k3 = kforms(thv)
    zv = mp.sin(thv/2)**2; cv = mp.cos(thv)
    Xi = mp.log(n*a_) + k1
    mu1 = Xi; mu2 = Xi**2+k2; mu3 = Xi**3+3*k2*Xi+k3
    return (2*zv*mu1/a_,
            zv*(A_*(mu2-2*mu1) - a_**2*mu2 + 2*a_*(1-zv))/a_**3,
            ( zv*(a_**4*mu3 - (3*A_*a_**2+B_*a_)*(mu3-3*mu2) + 3*A_**2*(mu3-4*mu2+2*mu1))/(3*a_**5)
              - 4*A_*zv*(1-zv)/a_**4 - cv*zv*(1-zv)/a_**3 ))
def make(which):
    if which=='exp':
        return (lambda x: (mp.exp(x)-1)/(mp.e-1), lambda x: mp.exp(x)/(mp.e-1), (mp.e/(mp.e-1),)*3)
    return (lambda x: (1+x)/2*mp.exp(x-1), lambda x: (2+x)/2*mp.exp(x-1), (mp.mpf('1.5'), mp.mpf(2), mp.mpf('2.5')))
def Dn_robust(which, thv, n):
    f, fp, _ = make(which)
    nn = mp.mpf(n); cv = mp.cos(thv); lnn = mp.log(nn)
    def bump(v):
        tt = mp.exp(v)/nn
        if tt > 1: return mp.mpf(0)
        fx = f(1-tt)
        if fx <= 0: return mp.mpf(0)
        lyp = mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fp(1-tt))
        return hc(mp.exp(lyp), cv)*mp.exp(v)/nn
    pts = sorted(set([-90, -30, lnn-8, lnn-2, lnn+2, lnn]))
    return mp.quad(bump, pts, maxdegree=10)
for which in ['exp','linexp']:
    a_,A_,B_ = make(which)[2]
    thv = 2*mp.pi/3
    n = 10**6
    D = Dn_robust(which, thv, n)
    J1,J2,J3 = Jpaper(thv, a_, A_, B_, mp.mpf(n))
    eps_ = 1/mp.mpf(n-1)
    sc = (D - (J1*eps_+J2*eps_**2+J3*eps_**3))/(eps_**4*mp.log(n)**4)
    print("robust: f=%6s th=2pi/3 n=1e6: scaled resid = %s (D=%s)" % (which, mp.nstr(sc,7), mp.nstr(D,14)))
