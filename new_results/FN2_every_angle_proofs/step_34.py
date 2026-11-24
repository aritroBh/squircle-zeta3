import mpmath as mp, time, json
mp.mp.dps = 70
t0 = time.time()

# ---------- stable h_c ----------
def hc(s, c):
    if s <= 0: return mp.mpf(0)
    if s < 1:
        return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
    r = mp.sqrt(1 + 2*c*s + s**2)
    rho = (1-c**2)/(r + s + c)     # = sqrt - s - c, stable
    return (1-c) - rho

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
    if which=='x':
        return (lambda x: x, lambda x: mp.mpf(1), (mp.mpf(1), mp.mpf(0), mp.mpf(0)))

def Dn(which, th, n):
    f, fp, _ = make_f(which)
    nn = mp.mpf(n); c = mp.cos(th)
    def yprime_logt(v):
        t = mp.exp(v)/nn
        fx = f(1-t)
        return nn*mp.exp((nn-1)*mp.log(fx))*fp(1-t)
    def integrand(v):
        return hc(yprime_logt(v), c)*mp.exp(v)/nn
    # window: t ~ log n/n, v ~ log log n; use breakpoints
    mid = mp.log(mp.log(nn))
    return mp.quad(integrand, [-140, mid-20, mid+20, mp.log(nn)], maxdegree=8)

# validate against AH3's printed values (result.md table + README)
print("validation vs AH3 printed values:")
checks = [(('x','pi/3',100), '0.02335054790245787127635'),
          (('x','2pi/3',1000), '0.01113690216868704954063'),
          (('exp','pi/2',100), '0.03372622791272318730492')]
angles = {'pi/3': mp.pi/3, 'pi/2': mp.pi/2, '2pi/3': 2*mp.pi/3}
for (w, thn, n), ref in checks:
    v = Dn(w, angles[thn], n)
    print(f"  f={w} th={thn} n={n}: {mp.nstr(v,25)}  vs {ref}  diff={mp.nstr(v-mp.mpf(ref),2)}")
print("elapsed", time.time()-t0)
