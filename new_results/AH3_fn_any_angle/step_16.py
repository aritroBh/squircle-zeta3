
import mpmath as mp
mp.mp.dps = 60

# ---- cumulants at angle theta (closed forms) ----
def S12(z):
    ell = mp.log(1-z)
    return mp.zeta(3) - mp.polylog(3, 1-z) + ell*mp.polylog(2, 1-z) + mp.mpf('0.5')*ell**2*mp.log(z)
def cumulants(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    k1 = 1 + u
    k2 = 1 + mp.pi**2/3 - 2*mp.polylog(2, z) - u**2
    k3 = -12*w*S12(z) + 2*(1 + u**3)
    return k1, k2, k3

# ---- J polynomials (from symbolic, hard-coded for speed) ----
def Jvals(th, a, A, B, n):
    k1, k2, k3 = cumulants(th)
    Xi = mp.log(n*a) + k1
    z = mp.sin(th/2)**2; c = mp.cos(th)
    mu1 = Xi; mu2 = Xi**2 + k2; mu3 = Xi**3 + 3*k2*Xi + k3
    J1 = 2*z*mu1/a
    # from sympy output:
    J2 = (A*mu2*z - 2*A*mu1*z - mu2*a**2*z + a*2*z*(1-z))/a**3
    # recompute J2 from the general formula to be safe:
    J2 = 2*z*(e1v(a,A)*mu1/a**2 - d2v(a,A)*mu2/a**3) + 2*z*(1-z)/a**2
    J3 = (2*z*( e1v(a,A)**2*mu1/a**3 + (e2v(a,A,B)/a**3 - 3*e1v(a,A)*d2v(a,A)/a**4)*mu2
              + (2*d2v(a,A)**2/a**5 - d3v(a,A,B)/a**4)*mu3 )
          + 2*z*(1-z)/a*(1/a - e1v(a,A)/a**2 - 2*d2v(a,A)/a**3) - c*z*(1-z)/a**3 )
    return J1, J2, J3
def d2v(a,A): return (a**2-A)/2
def d3v(a,A,B): return a**3/3 - a*A/2 + B/6
def e1v(a,A): return -A/a
def e2v(a,A,B): return B/(2*a) - A**2/(2*a**2)

def predict(th, a, A, B, n):
    eps = 1/(n-1)
    J1,J2,J3 = Jvals(th, a, A, B, n)
    return J1*eps + J2*eps**2 + J3*eps**3

# ---- direct quadrature of 2 - L_n = int h_c(y') dx, x = 1 - exp(v)/n ----
def hc(s, th):
    c = mp.cos(th)
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)

def make_f(which):
    if which=='x':
        f = lambda x: x; fp = lambda x: mp.mpf(1)
        data = (mp.mpf(1), mp.mpf(0), mp.mpf(0))
    elif which=='exp':
        e = mp.e
        f = lambda x: (mp.exp(x)-1)/(e-1); fp = lambda x: mp.exp(x)/(e-1)
        a = e/(e-1); data = (a, a, a)  # f' = f'' = f''' = e^x/(e-1) at 1
    elif which=='sin':
        f = lambda x: mp.sin(mp.pi*x/2); fp = lambda x: mp.pi/2*mp.cos(mp.pi*x/2)
        data = (mp.mpf(0), -mp.pi**2/4, mp.mpf(0))  # f'(1)=0! degenerate
    return f, fp, data

def Ln_deficit(which, th, n):
    f, fp, _ = make_f(which)
    nn = mp.mpf(n)
    def yprime(x):
        if x <= 0: return fp(0)*0 if which!='x' else mp.mpf(0)
        fx = f(x)
        if fx <= 0: return mp.mpf(0)
        return nn*mp.exp((nn-1)*mp.log(fx))*fp(x)
    def integrand(v):
        t = mp.exp(v)/nn
        x = 1 - t
        return hc(yprime(x), th)*mp.exp(v)/nn
    # v in (-inf, log n]; integrand ~ e^v as v-> -inf
    val = mp.quad(integrand, [-mp.inf, 0, mp.log(nn)], maxdegree=8)
    return val

# quick validation at n=1000
for which in ['x','exp']:
    for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
        a,A,B = make_f(which)[2]
        num = Ln_deficit(which, th, 1000)
        pre = predict(th, a, A, B, mp.mpf(1000))
        print(f" f={which} th={mp.nstr(th,6)}: 2-L={mp.nstr(num,16)} pred={mp.nstr(pre,16)} resid={mp.nstr(num-pre,4)}")
