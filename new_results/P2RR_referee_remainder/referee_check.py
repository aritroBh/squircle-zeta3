# Independent referee numerics for thm:P2asym remainder. Own code; does not import P2R files.
from mpmath import mp, mpf, pi, sin, cos, log, exp, hyp2f1, legendre, taylor, quad, sqrt, polylog, zeta, nstr
mp.dps = 80

def g(e, z):          # R(eps)/(2z) = (pi e/sin pi e) 2F1(-e,1+e;2;z)
    if e == 0: return mpf(1)
    return pi*e/sin(pi*e)*hyp2f1(-e, 1+e, 2, z)

def Lam_m(m, c):      # int_c^1 P_m
    return (legendre(m-1, c) - legendre(m+1, c))/(2*m+1)

def L_exact(th, n):   # thm:P2
    n = mpf(n); e = 1/(n-1); c = cos(th); z = sin(th/2)**2
    s = mpf(0)
    for m in range(1, 400):
        t = (-1)**m*e/(m-e)*Lam_m(m, c)*n**(-m)
        s += t
        if m > 4 and abs(t) < mpf(10)**(-mp.dps-10): break
    return 1 + c + 2*z*g(e, z)*n**(-e) + s

def L_quad(th, n):    # literal arc length of x u + x^n v
    c = cos(th); n = mpf(n)
    f = lambda x: sqrt(1 + 2*c*n*x**(n-1) + n**2*x**(2*n-2))
    pts = [mpf(0)] + [1 - mpf(10)**(-k) for k in range(1, int(log(n,10))+3)] + [mpf(1)]
    pts = sorted(set(pts))
    return quad(f, pts, maxdegree=10)

for th in (mpf(1), mpf('2.5')):
    z = sin(th/2)**2
    # a1, a2 numerically (Taylor of log g), independent of the closed forms
    co = taylor(lambda e: log(g(e, z)), mpf('0'), 3, method='quad', radius=mpf(1)/8)
    a1, a2, a3 = co[1], co[2], co[3]
    w = (1-z)/z; l = log(1-z); u = w*l
    a2cf = zeta(2) + mpf(1)/2 - polylog(2, z) - u**2/2
    print("theta=", nstr(th,6), "z=", nstr(z,12), "a1=", nstr(a1,20), "(-1-u=", nstr(-1-u,20), ") a2=", nstr(a2,20), "closed form diff", nstr(a2-a2cf,3))
    for n in (50, 400):
        print("  formula - quadrature at n=%d:" % n, nstr(L_exact(th, n) - L_quad(th, n), 3))
    mu = exp(a1); c2 = 2*z*(a2 - 1 + z); K = -2*z*a2
    print("  claimed log coefficient -2 z a2 =", nstr(K, 15))
    prev = None
    for k in range(2, 13):
        n = mpf(10)**k; e = 1/(n-1)
        E = L_exact(th, n) - (1 + cos(th) + 2*z*(mu/n)**e + c2/n**2)
        D = E*n**3/log(n)
        slope = '' if prev is None else nstr((E*n**3 - prev)/log(10), 15)
        print("  n=1e%-2d  E*n^3/log n = %s   slope of E*n^3 in log n = %s" % (k, nstr(D, 15), slope))
        prev = E*n**3
