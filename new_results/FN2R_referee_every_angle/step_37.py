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
def Jpaper(th, a, A, B, n):
    k1,k2,k3 = kforms(th)
    z = mp.sin(th/2)**2; c = mp.cos(th)
    Xi = mp.log(n*a) + k1
    mu1 = Xi; mu2 = Xi**2+k2; mu3 = Xi**3+3*k2*Xi+k3
    J1 = 2*z*mu1/a
    J2 = z*(A*(mu2-2*mu1) - a**2*mu2 + 2*a*(1-z))/a**3
    J3 = ( z*(a**4*mu3 - (3*A*a**2+B*a)*(mu3-3*mu2) + 3*A**2*(mu3-4*mu2+2*mu1))/(3*a**5)
           - 4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3 )
    return J1,J2,J3
# linexp with CORRECT derivatives: a=3/2, A=2, B=5/2
f  = lambda x: (1+x)/2*mp.exp(x-1)
fp = lambda x: (2+x)/2*mp.exp(x-1)
aa,AA,BB = mp.mpf('1.5'), mp.mpf(2), mp.mpf('2.5')
# verify derivatives numerically
for val, d, name in [(aa,1,'a'),(AA,2,'A'),(BB,3,'B')]:
    print(name, mp.nstr(mp.diff(f, mp.mpf(1), n=int(d)) - val, 2))

def Dn_linexp(th, n):
    nn = mp.mpf(n); c = mp.cos(th); lnn = mp.log(nn)
    def bump(v):
        tt = mp.exp(v)/nn
        fx = f(1-tt)
        if fx <= 0: return mp.mpf(0)
        lyp = mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fp(1-tt))
        yp = mp.exp(lyp)
        return hc(yp, c)*mp.exp(v)/nn
    pts = sorted(set([-90, -30, lnn-8, lnn-2, lnn+2, lnn+8, lnn+14, lnn]))
    return mp.quad(bump, pts, maxdegree=10)
def Ln_direct(th, n):
    """length = int sqrt(1+2c yp+yp^2) dx, then 2 - L - f(0)^n; independent route."""
    nn = mp.mpf(n); c = mp.cos(th); lnn = mp.log(nn)
    def integrand(v):
        tt = mp.exp(v)/nn
        fx = f(1-tt)
        if fx <= 0: return 1*mp.exp(v)/nn
        lyp = mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fp(1-tt))
        yp = mp.exp(lyp)
        if yp > mp.mpf('1e200'): return (yp + c)*mp.exp(v)/nn
        return mp.sqrt(1+2*c*yp+yp**2)*mp.exp(v)/nn
    pts = sorted(set([-90, -30, lnn-8, lnn-2, lnn+2, lnn+8, lnn+14, lnn]))
    L = mp.quad(integrand, pts, maxdegree=10)
    return 2 - L - mp.mpf('0.5')**0*f(mp.mpf(0))**n
angles = [('pi/3', mp.pi/3), ('2pi/3', 2*mp.pi/3), ('pi/2', mp.pi/2)]
rows=[]
for thn, th in angles:
    for n in [10**3, 10**4, 10**5, 10**6]:
        D = Dn_linexp(th, n)
        D2 = Ln_direct(th, n)
        J1,J2,J3 = Jpaper(th, aa, AA, BB, mp.mpf(n))
        eps = 1/mp.mpf(n-1)
        pre = J1*eps + J2*eps**2 + J3*eps**3
        resid = D - pre
        scaled = resid/(eps**4*mp.log(n)**4)
        rows.append((thn,n,mp.nstr(D,40),mp.nstr(resid,10),mp.nstr(scaled,10),mp.nstr(abs(D-D2),3)))
        print("linexp th=%5s n=%7d: D=%s |D-Dlen|=%s resid=%s scaled=%s" %
              (thn, n, mp.nstr(D,10), mp.nstr(abs(D-D2),2), mp.nstr(resid,3), mp.nstr(scaled,6)))
json.dump(rows, open("fn2r_residuals_linexp.json","w"), indent=1)
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f2:
    f2.write(f"{time.strftime('%H:%M')} step_08b: linexp corrected (a=3/2,A=2,B=5/2 — my slip fixed); cross-check vs direct length integral; scaled residuals [referee, NUMERICAL]\n")
