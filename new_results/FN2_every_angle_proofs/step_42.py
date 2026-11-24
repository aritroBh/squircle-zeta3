import mpmath as mp, time, json
mp.mp.dps = 60
t0 = time.time()
def hc(s, c):
    if s <= 0: return mp.mpf(0)
    if s < 1: return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c) - (1-c**2)/(r + s + c)
def phi_st(s, c):
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c**2)/(r*(r+s+c))
def S12(z):
    return mp.zeta(3) - mp.polylog(3, 1-z) + mp.log(1-z)*mp.polylog(2, 1-z) \
           + mp.mpf('0.5')*mp.log(1-z)**2*mp.log(z)
def cumulants(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    k1 = 1+u; k2 = 1+mp.pi**2/3-2*mp.polylog(2,z)-u**2; k3 = -12*w*S12(z)+2*(1+u**3)
    return k1,k2,k3
def k4_numeric(th):
    c = mp.cos(th)
    m = []
    for r in range(1,5):
        i1 = mp.quad(lambda u: phi_st(mp.exp(-u),c)*mp.exp(-u)*u**r, [0, mp.inf])
        i2 = mp.quad(lambda u: phi_st(mp.exp(u),c)*mp.exp(u)*(-u)**r, [0, mp.inf])
        m.append((i1+i2)/(1-c))
    m1,m2,m3,m4 = m
    return m4-4*m3*m1-3*m2**2+12*m2*m1**2-6*m1**4

def make_f(which):
    if which=='exp':
        e = mp.e; a = e/(e-1)
        return (lambda x: (mp.exp(x)-1)/(e-1), lambda x: mp.exp(x)/(e-1), (a,a,a,a))
    if which=='quad':
        return (lambda x: (x+x**2)/2, lambda x: (1+2*x)/2, (mp.mpf('1.5'), mp.mpf(1), mp.mpf(0), mp.mpf(0)))

def theta_ij(a,A,B,D4):
    return {
     (1,0): 1/a, (2,0): (A-a**2)/(2*a**3), (1,1): -A/a**3,
     (3,0): (3*A**2-3*A*a**2-B*a+a**4)/(6*a**5),
     (2,1): (-4*A**2+3*A*a**2+B*a)/(2*a**5), (1,2): A**2/a**5,
     (4,0): -(a**6-7*a**4*A-6*a**3*B+18*a**2*A**2-a**2*D4+10*a*A*B-15*A**3)/(24*a**7),
     (3,1): -(7*a**4*A+6*a**3*B-24*a**2*A**2+a**2*D4-13*a*A*B+23*A**3)/(6*a**7),
     (2,2): -3*A*(2*a**2*A+a*B-3*A**2)/(2*a**7), (1,3): -A**3/a**7 }

def Jvals_all(th, prms, n, k4):
    a,A,B,D4 = prms
    k1,k2,k3 = cumulants(th)
    Xi = mp.log(n*a)+k1
    z = mp.sin(th/2)**2; c = mp.cos(th)
    mu1 = Xi; mu2 = Xi**2+k2; mu3 = Xi**3+3*k2*Xi+k3; mu4 = Xi**4+6*k2*Xi**2+4*k3*Xi+k4+3*k2**2
    J1 = 2*z*mu1/a
    J2 = z*(A*(mu2-2*mu1) - a**2*mu2 + 2*a*(1-z))/a**3
    J3 = ( z*(a**4*mu3 - (3*A*a**2+B*a)*(mu3-3*mu2) + 3*A**2*(mu3-4*mu2+2*mu1))/(3*a**5)
           - 4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3 )
    # J4 via master formula
    t = theta_ij(a,A,B,D4)
    M1,M2,M3,M4 = 2*z*mu1, 2*z*mu2, 2*z*mu3, 2*z*mu4
    d2 = 2*z*(1-z); d3 = -4*c*z*(1-z); d4c = 3*(1-c**2)*(5*c**2-1)/8
    T4 = ( d2/a*(6*t[(3,0)]-2*t[(2,1)]+t[(1,2)]+2*t[(2,0)]-t[(1,1)]+t[(1,0)])
         + d3/a**2*(-t[(2,0)]/4+t[(1,1)]/4-t[(1,0)]/2)
         + t[(1,0)]*d4c/(9*a**3) )
    J4 = t[(4,0)]*M4 + t[(3,1)]*M3 + t[(2,2)]*M2 + t[(1,3)]*M1 + T4
    return J1,J2,J3,J4

def Dn(which, th, n):
    f, fp, _ = make_f(which)[:2], make_f(which)[1], None
    f = make_f(which)[0]; fp = make_f(which)[1]
    nn = mp.mpf(n); c = mp.cos(th)
    def integrand(v):
        tt = mp.exp(v)/nn
        if tt >= 1: return mp.mpf(0)
        fx = f(1-tt)
        if fx <= 0: return mp.mpf(0)
        yp = nn*mp.exp((nn-1)*mp.log(fx))*fp(1-tt)
        return hc(yp, c)*mp.exp(v)/nn
    mid = mp.log(mp.log(nn)); lnn = mp.log(nn)
    pts = [-160, mid-6, mid-2, mid+2, mid+6, min(mid+12, lnn-mp.mpf('0.5')), lnn]
    return mp.quad(integrand, pts, maxdegree=10)

angles = {'pi/3': mp.pi/3, 'pi/2': mp.pi/2, '2pi/3': 2*mp.pi/3}
print("K=4 check: resid3 vs eps^4*J4 (ratio ~ 1+O(eps log n)); resid4 scaled by eps^5 log^5 n")
k4s = {name: k4_numeric(th) for name, th in angles.items()}
for name, k4v in k4s.items(): print(f"  k4({name}) = {mp.nstr(k4v,12)}")
rows = []
for which in ['exp','quad']:
    prms = make_f(which)[2]
    for thn in ['pi/3','2pi/3','pi/2']:
        th = angles[thn]; k4 = k4s[thn]
        for n in [10**4, 10**5, 10**6]:
            num = Dn(which, th, n)
            J1,J2,J3,J4 = Jvals_all(th, prms, mp.mpf(n), k4)
            eps = 1/mp.mpf(n-1)
            resid3 = num - (J1*eps+J2*eps**2+J3*eps**3)
            resid4 = resid3 - J4*eps**4
            ratio = resid3/(J4*eps**4)
            sc4 = resid4/(eps**5*mp.log(n)**5)
            rows.append((which,thn,n,mp.nstr(ratio,8),mp.nstr(sc4,8)))
            print(f"f={which:4s} th={thn:6s} n={n:7d}: resid3/(eps^4 J4)={mp.nstr(ratio,7)}  resid4/(eps^5 log^5 n)={mp.nstr(sc4,5)}")
json.dump(rows, open('residuals_K4.json','w'), indent=1)
print("elapsed", time.time()-t0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} K=4 check: resid after J3 = eps^4*J4*(1+o(1)) (ratios ~1), resid after J4 scaled by eps^5 log^5 n bounded  [NUMERICAL]\n")
