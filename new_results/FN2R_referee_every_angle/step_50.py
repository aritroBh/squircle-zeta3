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
def S12(zv):
    return mp.zeta(3) - mp.polylog(3, 1-zv) + mp.log(1-zv)*mp.polylog(2, 1-zv) \
           + mp.mpf('0.5')*mp.log(1-zv)**2*mp.log(zv)
def kforms(th):
    zv = mp.sin(th/2)**2; w2 = 1/mp.tan(th/2)**2; u = w2*mp.log(mp.cos(th/2)**2)
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,zv)-u**2, 2+2*u**3-12*w2*S12(zv))
def k4_quad(th):
    c = mp.cos(th)
    raw=[]
    for r in (1,2,3,4):
        i1 = mp.quad(lambda u: phi_s(mp.exp(-u),c)*mp.exp(-u)*u**r, [0,mp.inf])
        i2 = mp.quad(lambda u: phi_s(mp.exp(u),c)*mp.exp(u)*(-u)**r, [0,mp.inf])
        raw.append((i1+i2)/(1-c))
    m1,m2,m3,m4 = raw
    k2 = m2-m1**2; k3 = m3-3*m1*m2+2*m1**3
    k4 = m4 - 4*m3*m1 - 3*m2**2 + 12*m2*m1**2 - 6*m1**4
    return m1,k2,k3,k4
# J's including J4 (master form)
def Jvals(th, a, A, B, n, k4):
    k1,k2,k3 = kforms(th)
    z = mp.sin(th/2)**2; c = mp.cos(th)
    Xi = mp.log(n*a) + k1
    mu1 = Xi; mu2 = Xi**2+k2; mu3 = Xi**3+3*k2*Xi+k3
    mu4 = Xi**4 + 6*k2*Xi**2 + 4*k3*Xi + (k4+3*k2**2)
    M1,M2,M3,M4 = 2*z*mu1, 2*z*mu2, 2*z*mu3, 2*z*mu4
    J1 = 2*z*mu1/a
    J2 = z*(A*(mu2-2*mu1) - a**2*mu2 + 2*a*(1-z))/a**3
    J3 = ( z*(a**4*mu3 - (3*A*a**2+B*a)*(mu3-3*mu2) + 3*A**2*(mu3-4*mu2+2*mu1))/(3*a**5)
           - 4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3 )
    th40 = (15*A**3 - 10*A*B*a + 7*A*a**4 + 6*B*a**3 - a**6 + a**2*(-18*A**2 + (a*0) ))/(24*a**7)
    return J1,J2,J3
def J4val(th, a, A, B, D4, n, k4):
    k1,k2,k3 = kforms(th)
    z = mp.sin(th/2)**2; c = mp.cos(th)
    Xi = mp.log(n*a) + k1
    mu1 = Xi; mu2 = Xi**2+k2; mu3 = Xi**3+3*k2*Xi+k3
    mu4 = Xi**4 + 6*k2*Xi**2 + 4*k3*Xi + (k4+3*k2**2)
    M1,M2,M3,M4 = 2*z*mu1, 2*z*mu2, 2*z*mu3, 2*z*mu4
    th40 = (15*A**3 - 10*A*B*a + 7*A*a**4 + 6*B*a**3 - a**6 + a**2*(-18*A**2 + D4))/(24*a**7)
    th31 = -(7*a**4*A+6*a**3*B-24*a**2*A**2+a**2*D4-13*a*A*B+23*A**3)/(6*a**7)
    th22 = -3*A*(2*a**2*A+a*B-3*A**2)/(2*a**7)
    th13 = -A**3/a**7
    T4 = -(c-1)*(c+1)*(96*A**2 - 48*A*a**2 + 9*A*a*c - 24*B*a + 12*a**4 + 9*a**3*c + 5*a**2*c**2 - a**2)/(24*a**6)
    return th40*M4 + th31*M3 + th22*M2 + th13*M1 + T4

def make(which):
    if which=='exp':
        f  = lambda x: (mp.exp(x)-1)/(mp.e-1)
        fp = lambda x: mp.exp(x)/(mp.e-1)
        aa = mp.e/(mp.e-1)
        return f, fp, (aa, aa, aa, aa)
    if which=='linexp':
        f  = lambda x: (1+x)/2*mp.exp(x-1)
        fp = lambda x: (2+x)/2*mp.exp(x-1)
        return f, fp, (mp.mpf('1.5'), mp.mpf(2), mp.mpf('2.5'), mp.mpf(3))
def Dn(which, th, n):
    f, fp, _ = make(which)
    nn = mp.mpf(n); c = mp.cos(th); lnn = mp.log(nn)
    def bump(v):
        tt = mp.exp(v)/nn
        if tt > 1: return mp.mpf(0)
        fx = f(1-tt)
        if fx <= 0: return mp.mpf(0)
        lyp = mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fp(1-tt))
        return hc(mp.exp(lyp), c)*mp.exp(v)/nn
    pts = sorted(set([-90, -30, lnn-8, lnn-2, lnn+2, lnn]))
    return mp.quad(bump, pts, maxdegree=10)
for which in ['exp','linexp']:
    a,A,B,D4 = make(which)[2]
    for thn, th in [('pi/3', mp.pi/3), ('2pi/3', 2*mp.pi/3)]:
        k1q,k2q,k3q,k4 = k4_quad(th)
        # sanity: k4 vs pi/2 value 6+pi^4/15 at pi/2 only; here just use
        for n in [10**4, 10**5, 10**6]:
            D = Dn(which, th, n)
            J1,J2,J3 = Jvals(th, a, A, B, mp.mpf(n), k4)
            J4 = J4val(th, a, A, B, D4, mp.mpf(n), k4)
            eps = 1/mp.mpf(n-1)
            resid = D - (J1*eps+J2*eps**2+J3*eps**3)
            print("f=%6s th=%5s n=%7d: resid/(eps^4 J4)=%s  J4=%s" % (
                which, thn, n, mp.nstr(resid/(eps**4*J4),8), mp.nstr(J4,6)))
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f2:
    f2.write(f"{time.strftime('%H:%M')} step_12b: resid/(eps^4 J4) with my own theta_4j, T4(c), k4-quadrature: tested exp & linexp at pi/3, 2pi/3, n=1e4..1e6 [referee, NUMERICAL]\n")
