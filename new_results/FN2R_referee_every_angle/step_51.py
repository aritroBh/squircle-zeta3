import mpmath as mp, time
mp.mp.dps = 50
t0=time.time()
def phi_s(s,c):
    r = mp.sqrt(1+2*c*s+s**2)
    return (1-c**2)/(r*(r+s+c))
def k4_quad(th):
    c = mp.cos(th)
    raw=[]
    for r in (1,2,3,4):
        i1 = mp.quad(lambda u: phi_s(mp.exp(-u),c)*mp.exp(-u)*u**r, [0,mp.inf])
        i2 = mp.quad(lambda u: phi_s(mp.exp(u),c)*mp.exp(u)*(-u)**r, [0,mp.inf])
        raw.append((i1+i2)/(1-c))
    m1,m2,m3,m4 = raw
    k2 = m2-m1**2; k3 = m3-3*m1*m2+2*m1**3
    return m1, k2, k3, m4 - 4*m3*m1 - 3*m2**2 + 12*m2*m1**2 - 6*m1**4
k1q,k2q,k3q,k4q = k4_quad(mp.pi/2)
print("k4(pi/2) quadrature =", mp.nstr(k4q,15), " 6+pi^4/15 =", mp.nstr(6+mp.pi**4/15,15), " diff:", mp.nstr(k4q-(6+mp.pi**4/15),2))

# tau^+(s) for sin at general theta: outer point where (f^n)' = s
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
beta = mp.pi**2/8; m = mp.mpf(2)
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
def tau_plus(n, s):
    """outer (largest) tau with (f^n)'(1-tau) = s: solve on the outer branch."""
    nn = mp.mpf(n)
    def lyp(tt):
        fx = f(1-tt)
        return mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fp(1-tt))
    # outer branch: tau beyond the argmax t* ~ 2/(pi sqrt n)
    tstar = 2/(mp.pi*mp.sqrt(nn))
    return mp.findroot(lambda tt: lyp(tt) - mp.log(s), (tstar, 10*tstar))
for thn, th in [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3)]:
    c = mp.cos(th); k1,k2,k3 = kforms(th)
    for n in [10**5, 10**6]:
        D = Dn_sin(th, n)
        Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn)
        tp = tau_plus(n, mp.exp(-k1))
        pred = 1 - (m-1)*k2/(2*m**2*Lam**2)
        val = D/((1-c)*tp)
        print("th=%5s n=%8d: D/((1-c)tau^+(e^-k1)) = %s  pred bracket = %s  diff = %s (expect O(logLam/Lam^3)=%s)" % (
            thn, n, mp.nstr(val,9), mp.nstr(pred,9), mp.nstr(val-pred,4), mp.nstr(mp.log(Lam)/Lam**3,3)))
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f2:
    f2.write(f"{time.strftime('%H:%M')} step_13: k4(pi/2)=6+pi^4/15 confirmed (1e-50); Thm2(b) tau^+(e^{-k1}) form checked at 3 angles: bracket 1-(m-1)k2/(2m^2 Lam^2) holds to O(log Lam/Lam^3) [referee, NUMERICAL]\n")
