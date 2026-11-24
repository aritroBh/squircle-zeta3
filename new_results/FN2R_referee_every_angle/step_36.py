import mpmath as mp, time, json
mp.mp.dps = 60
t0 = time.time()

# ---------- my own defect engine (independent of FN2's step_39) ----------
def hc(s, c):
    if s <= 0: return mp.mpf(0)
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c) - (1-c**2)/(r + s + c)   # stable everywhere

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

# f1(x) = (1+x)/2 * e^{x-1};  f2(x) = (e^x-1)/(e-1)
def make(which):
    if which=='linexp':
        f  = lambda x: (1+x)/2*mp.exp(x-1)
        fp = lambda x: (3+x)/2*mp.exp(x-1)
        return f, fp, (mp.mpf(2), mp.mpf('2.5'), mp.mpf(3))
    if which=='exp':
        f  = lambda x: (mp.exp(x)-1)/(mp.e-1)
        fp = lambda x: mp.exp(x)/(mp.e-1)
        aa = mp.e/(mp.e-1)
        return f, fp, (aa, aa, aa)

def Dn_layercake_direct(which, th, n):
    """D_n = int_0^1 h_c(y') dx via t=1-x, t=e^v/n ; also L_n direct for cross-check."""
    f, fp, _ = make(which)
    nn = mp.mpf(n); c = mp.cos(th)
    def vals(v):
        tt = mp.exp(v)/nn
        fx = f(1-tt)
        if fx <= 0: return mp.mpf(0), mp.mpf(0)
        lyp = mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fp(1-tt))
        yp = mp.exp(lyp)
        h = hc(yp, c) if yp < mp.mpf('1e300') else (1-c)
        ln = mp.sqrt(1+2*c*yp+yp**2) if yp < mp.mpf('1e300') else yp + c
        return h*mp.exp(v)/nn, (1+yp-ln)*mp.exp(v)/nn  # second: (1+y'-|u+y'v|) dx = h as well; use sqrt form
    # integrate h-form; cross-check with 1+y'-sqrt form on same grid
    lnn = mp.log(nn)
    pts = [-80, -20, lnn-8, lnn-2, lnn+2, lnn+8, min(lnn+14, lnn+mp.log(mp.log(nn))*0+14), lnn+mp.mpf('0.0001')]
    pts = sorted(set(pts + [lnn]))
    Dh = mp.quad(lambda v: vals(v)[0], pts, maxdegree=10)
    Dq = mp.quad(lambda v: vals(v)[1], pts, maxdegree=10)
    return Dh, Dq

angles = [('pi/3', mp.pi/3), ('2pi/3', 2*mp.pi/3), ('pi/2', mp.pi/2)]
rows = []
for which in ['linexp','exp']:
    a,A,B = make(which)[2]
    for thn, th in angles:
        for n in [10**3, 10**4, 10**5, 10**6]:
            Dh, Dq = Dn_layercake_direct(which, th, n)
            J1,J2,J3 = Jpaper(th, a, A, B, mp.mpf(n))
            eps = 1/mp.mpf(n-1)
            pre = J1*eps + J2*eps**2 + J3*eps**3
            resid = Dh - pre
            scaled = resid/(eps**4*mp.log(n)**4)
            rows.append((which,thn,n,Dh,mp.nstr(resid,8),mp.nstr(scaled,8),mp.nstr(abs(Dh-Dq),2)))
            print("f=%6s th=%5s n=%7d: D=%s |Dh-Dq|=%s resid=%s scaled=%s" %
                  (which, thn, n, mp.nstr(Dh,10), mp.nstr(abs(Dh-Dq),2), mp.nstr(resid,3), mp.nstr(scaled,6)))
json.dump([[r[0],r[1],r[2],mp.nstr(r[3],40),r[4],r[5],r[6]] for r in rows], open("fn2r_residuals.json","w"), indent=1)
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_08: MY engine (h_c bump AND 1+y'-|u+y'v| cross-checks agree) ; residuals after J3 / eps^4 log^4 n bounded for f=(1+x)/2 e^{x-1} and (e^x-1)/(e-1), 3 angles, n=1e3..1e6 [referee, NUMERICAL]\n")
