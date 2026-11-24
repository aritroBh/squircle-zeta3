"""
contact_all_m.py -- Theorem 4 for every contact exponent m > 0, m != 1 (including INFINITE slope at x=1, m<1,
and non-integer m), checked on
   f = 1 - sqrt(1-x)        (m = 1/2, beta = 1:  f'(1) = +infinity)
   f = 1 - (1-x)^{3/2}      (m = 3/2, beta = 1:  f'(1) = 0, non-integer contact)
Model:  D_n = N^{-1} E[ (U^{1/m} - V^{1/m})(Lam + xi) ]  (V absent when m < 1),
        N = (beta n)^{1/m}, Lam = log(m N), U - c ln U = y, c = 1 - 1/m.
   m > 1:  U = -c W_{-1}(-e^{-y/c}/c),  V = -c W_0(-e^{-y/c}/c)
   m < 1:  U = |c| W_0(e^{y/|c|}/|c|)   (unique root; slope of f^n is infinite at x = 1)
Universal law: D_n ~ (log n/(m beta n))^{1/m}; slope-2/e: D_n = tau^+(2/e) [1 - (m-1) k2/(2 m^2 Lam^2) + ...].
"""
import mpmath as mp
from fnlib import g, gprime

K1 = 1 - mp.log(2)
K2 = 1 + mp.pi**2/6

def U_of(y, m):
    c = 1 - mp.mpf(1)/m
    if m > 1:
        return -c*mp.lambertw(-mp.exp(-y/c)/c, -1).real
    cc = -c
    return cc*mp.lambertw(mp.exp(y/cc)/cc, 0).real

def V_of(y, m):
    c = 1 - mp.mpf(1)/m
    return -c*mp.lambertw(-mp.exp(-y/c)/c, 0).real

def model(Lam, m):
    w = lambda xi: gprime(mp.exp(-xi))*mp.exp(-xi)
    if m > 1:
        c = 1 - mp.mpf(1)/m
        ymin = c*(1 - mp.log(c))
        lo = ymin - Lam
        H = lambda xi: U_of(Lam + xi, m)**(1/mp.mpf(m)) - V_of(Lam + xi, m)**(1/mp.mpf(m))
        pts = [lo, lo + mp.mpf('0.01'), lo + 1, lo + 4]
        pts = [p for p in pts if p < 0] + [mp.mpf(0), mp.mpf(2), mp.mpf(10), mp.mpf(40), mp.inf]
    else:
        H = lambda xi: U_of(Lam + xi, m)**(1/mp.mpf(m))
        pts = [-mp.inf, -Lam - 10, -Lam, -Lam/2, mp.mpf(-4), mp.mpf(0), mp.mpf(2), mp.mpf(10), mp.mpf(40), mp.inf]
    return mp.quad(lambda xi: w(xi)*H(xi), sorted(set(pts)))

def defect(f, fp, n, m, beta):
    n = mp.mpf(n); N = (beta*n)**(1/mp.mpf(m)); nm1 = n - 1
    Lam = mp.log(m*N)
    def integrand(v):
        t = v/N
        if t <= 0:
            return mp.mpf(1) if m < 1 else mp.mpf(0)
        fx = f(t)
        if fx <= 0:
            return mp.mpf(0)
        return g(n*mp.exp(nm1*mp.log(fx))*fp(t))
    vmax = min((Lam + (mp.mp.dps + 20)*mp.log(10))**(1/mp.mpf(m))*mp.mpf('1.2'), N)
    c0 = Lam**(1/mp.mpf(m))
    pts = [mp.mpf(0)]
    if m > 1:
        inner = N**(-1/(mp.mpf(m) - 1))
        p = inner*mp.mpf(10)**-3
        while p < mp.mpf('0.5'):
            pts.append(p); p *= 4
    else:
        p = mp.mpf(10)**-12
        while p < mp.mpf('0.5'):
            pts.append(p); p *= 8
    pts += [mp.mpf('0.5'), mp.mpf(1)]
    for k in [-0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 2]:
        q = c0*(1 + mp.mpf(k))
        if 1 < q < vmax:
            pts.append(q)
    pts.append(vmax)
    return mp.quad(integrand, sorted(set(pts)))/N, N, Lam

def tau_plus(f, fp, n, s, m, beta):
    """outer root tau of (f^n)'(1-tau) = s, by bisection in log tau (slope decreasing beyond its maximum)"""
    n = mp.mpf(n)
    F = lambda t: (n - 1)*mp.log(f(t)) + mp.log(fp(t)) + mp.log(n) - mp.log(s)   # >0 inside, <0 outside
    N = (beta*n)**(1/mp.mpf(m))
    t0 = (mp.log(m*N) + mp.log(1/s))**(1/mp.mpf(m))/N
    lo, hi = t0/50, min(t0*50, mp.mpf('0.999'))
    assert F(lo) > 0 and F(hi) < 0
    a_, b_ = mp.log(lo), mp.log(hi)
    for _ in range(mp.mp.prec + 20):
        mid = (a_ + b_)/2
        if F(mp.exp(mid)) > 0:
            a_ = mid
        else:
            b_ = mid
    return mp.exp((a_ + b_)/2)

CASES = [
    # functions of t = 1 - x:  f(1-t), f'(1-t)
    ('1-sqrt(1-x)', lambda t: 1 - mp.sqrt(t), lambda t: 1/(2*mp.sqrt(t)), mp.mpf(1)/2, 1),
    ('1-(1-x)^(3/2)', lambda t: 1 - t**(mp.mpf(3)/2), lambda t: mp.mpf(3)/2*mp.sqrt(t), mp.mpf(3)/2, 1),
]
if __name__ == '__main__':
    for name, f, fp, m, beta in CASES:
        print('=' * 104)
        print('f = %s  (m=%s, beta=%s)' % (name, mp.nstr(m, 3), beta))
        print('%-7s %-26s %-12s %-12s %-12s %-12s' % ('n', 'D_n', 'model 4(a)', 'tau+(2/e)', 'pred. corr', 'lead law'))
        for e in [4, 8, 12, 20, 30, 40]:
            mp.mp.dps = 50 + 3*e
            D, N, Lam = defect(f, fp, mp.mpf(10)**e, m, mp.mpf(beta))
            mod = model(Lam, m)/N
            tp = tau_plus(f, fp, mp.mpf(10)**e, 2/mp.e, m, mp.mpf(beta))
            corr = -(m - 1)*K2/(2*m**2*Lam**2)
            lead = (mp.log(mp.mpf(10)**e)/(m*beta*mp.mpf(10)**e))**(1/mp.mpf(m))
            print('1e%-5d %-26s %-12s %-12s %-12s %-12s' % (e, mp.nstr(D, 18), mp.nstr(mod/D - 1, 3),
                  mp.nstr(D/tp - 1, 3), mp.nstr(corr, 3), mp.nstr(lead/D - 1, 3)))
        print('  model 4(a), tau+(2/e), lead law: relative errors; "D/tau+ - 1" should approach "pred. corr" = -(m-1)k2/(2m^2 Lam^2)')
