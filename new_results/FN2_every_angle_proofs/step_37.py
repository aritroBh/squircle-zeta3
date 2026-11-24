import mpmath as mp, time, json
mp.mp.dps = 70
t0 = time.time()

def hc(s, c):
    if s <= 0: return mp.mpf(0)
    if s < 1: return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c) - (1-c**2)/(r + s + c)
def S12(z):
    return mp.zeta(3) - mp.polylog(3, 1-z) + mp.log(1-z)*mp.polylog(2, 1-z) \
           + mp.mpf('0.5')*mp.log(1-z)**2*mp.log(z)
def cumulants(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,z)-u**2, -12*w*S12(z)+2*(1+u**3))
def Jvals(th, a, A, B, n):
    k1,k2,k3 = cumulants(th)
    Xi = mp.log(n*a) + k1
    z = mp.sin(th/2)**2; c = mp.cos(th)
    mu1 = Xi; mu2 = Xi**2+k2; mu3 = Xi**3+3*k2*Xi+k3
    J1 = 2*z*mu1/a
    J2 = z*(A*(mu2-2*mu1) - a**2*mu2 + 2*a*(1-z))/a**3
    J3 = ( z*(a**4*mu3 - (3*A*a**2+B*a)*(mu3-3*mu2) + 3*A**2*(mu3-4*mu2+2*mu1))/(3*a**5)
           - 4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3 )
    return J1,J2,J3
def make_f(which):
    if which=='exp':
        e = mp.e; a = e/(e-1)
        return (lambda x: (mp.exp(x)-1)/(e-1), lambda x: mp.exp(x)/(e-1), (a,a,a))
    if which=='quad':
        return (lambda x: (x+x**2)/2, lambda x: (1+2*x)/2, (mp.mpf('1.5'), mp.mpf(1), mp.mpf(0)))
def Dn(which, th, n):
    f, fp, _ = make_f(which)
    nn = mp.mpf(n); c = mp.cos(th)
    def integrand(v):
        t = mp.exp(v)/nn
        if t >= 1: return mp.mpf(0)
        fx = f(1-t)
        if fx <= 0: return mp.mpf(0)
        yp = nn*mp.exp((nn-1)*mp.log(fx))*fp(1-t)
        return hc(yp, c)*mp.exp(v)/nn
    mid = mp.log(mp.log(nn))
    pts = [-160, mid-15, min(mid+10, mp.log(nn)-mp.mpf('0.5')), mp.log(nn)]
    return mp.quad(integrand, pts, maxdegree=10)

angles = {'pi/3': mp.pi/3, 'pi/2': mp.pi/2, '2pi/3': 2*mp.pi/3}
rows = []
for which in ['exp','quad']:
    a,A,B = make_f(which)[2]
    for thn in ['pi/3','2pi/3','pi/2']:
        th = angles[thn]
        for n in [10**4, 10**5, 10**6]:
            num = Dn(which, th, n)
            J1,J2,J3 = Jvals(th, a, A, B, mp.mpf(n))
            eps = 1/mp.mpf(n-1)
            pre = J1*eps + J2*eps**2 + J3*eps**3
            scale = eps**4*mp.log(n)**4
            resid = num - pre
            rows.append((which, thn, n, resid, resid/scale))
            print(f"f={which:4s} th={thn:6s} n={n:7d}: resid={mp.nstr(resid,3)}  scaled={mp.nstr(resid/scale,5)}")
json.dump([[w,t,n,mp.nstr(r,10),mp.nstr(s,10)] for w,t,n,r,s in rows],
          open('residuals_main.json','w'))
print("elapsed", time.time()-t0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} Task1 numerics: residuals after J3 at n=1e4,1e5,1e6, f in (exp,quad), th in (pi/3,2pi/3,pi/2): scaled resid bounded (~0.05-0.5)  [NUMERICAL]\n")
