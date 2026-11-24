
import mpmath as mp, json, time
mp.mp.dps = 55

def S12(z):
    ell = mp.log(1-z)
    return mp.zeta(3) - mp.polylog(3, 1-z) + ell*mp.polylog(2, 1-z) + mp.mpf('0.5')*ell**2*mp.log(z)
def cumulants(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    k1 = 1 + u
    k2 = 1 + mp.pi**2/3 - 2*mp.polylog(2, z) - u**2
    k3 = -12*w*S12(z) + 2*(1 + u**3)
    return k1, k2, k3
def hc(s, th):
    c = mp.cos(th)
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)

def Jvals(th, a, A, B, n):
    k1, k2, k3 = cumulants(th)
    Xi = mp.log(n*a) + k1
    z = mp.sin(th/2)**2; c = mp.cos(th)
    mu1 = Xi; mu2 = Xi**2 + k2; mu3 = Xi**3 + 3*k2*Xi + k3
    J1 = 2*z*mu1/a
    J2 = z*(A*(mu2 - 2*mu1) - a**2*mu2 + 2*a*(1-z))/a**3
    J3 = ( z*( a**4*mu3 - (3*A*a**2 + B*a)*(mu3 - 3*mu2) + 3*A**2*(mu3 - 4*mu2 + 2*mu1) )/(3*a**5)
           - 4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3 )
    return J1, J2, J3
def predict(th, a, A, B, n):
    eps = 1/(n-1); J1,J2,J3 = Jvals(th, a, A, B, mp.mpf(n))
    return J1*eps + J2*eps**2 + J3*eps**3

def make_f(which):
    if which=='x':
        return (lambda x: x, lambda x: mp.mpf(1), (mp.mpf(1), mp.mpf(0), mp.mpf(0)))
    if which=='exp':
        e = mp.e; a = e/(e-1)
        return (lambda x: (mp.exp(x)-1)/(e-1), lambda x: mp.exp(x)/(e-1), (a, a, a))
    if which=='sin':
        return (lambda x: mp.sin(mp.pi*x/2), lambda x: mp.pi/2*mp.cos(mp.pi*x/2), (mp.mpf(0), -mp.pi**2/4, mp.mpf(0)))

def quad_deficit(which, th, n):
    f, fp, _ = make_f(which)
    nn = mp.mpf(n)
    def yprime(x):
        fx = f(x)
        if fx <= 0: return mp.mpf(0)
        return nn*mp.exp((nn-1)*mp.log(fx))*fp(x)
    def integrand(v):
        t = mp.exp(v)/nn
        return hc(yprime(1-t), th)*mp.exp(v)/nn
    if which == 'sin':
        mid = 0.5*mp.log(nn)
        return mp.quad(integrand, [-80, mid-15, mid+5, mp.log(nn)], maxdegree=8)
    return mp.quad(integrand, [-130, 0, 3, mp.log(nn)], maxdegree=8)

t0=time.time()
results = {}
angles = [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3)]
for which in ['x','exp']:
    a,A,B = make_f(which)[2]
    for name, th in angles:
        for n in [100, 1000, 10000, 100000]:
            num = quad_deficit(which, th, n)
            pre = predict(th, a, A, B, n)
            eps = 1/mp.mpf(n-1)
            scale = eps**4*mp.log(n)**4
            resid = num - pre
            results[(which,name,n)] = (num, pre, resid, resid/scale)
            print(f"f={which:3s} th={name:6s} n={n:6d}: 2-L={mp.nstr(num,17)} resid={mp.nstr(resid,3)} scaled={mp.nstr(resid/scale,5)}")
print('elapsed', time.time()-t0)
json.dump({str(k): [mp.nstr(x,40) for x in v] for k,v in results.items()}, open('results_main.json','w'))
