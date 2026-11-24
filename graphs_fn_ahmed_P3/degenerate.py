"""
degenerate.py -- the case f'(1) = 0 (e.g. f = sin(pi x/2)), not covered by Ahmed's hypothesis.

If f(1-t) = 1 - beta t^m + O(t^{m+1}) with m >= 2, beta > 0, then (Theorem 4 in result.md)
  D_n = N^{-1} E[ F(Lam + xi) ] * (1 + O(power-small)),     N = (beta n)^{1/m},  Lam = log(m N),
  F(y) = U(y)^{1/m},  U = larger root of  U - (1 - 1/m) log U = y   (U = -c W_{-1}(-e^{-y/c}/c), c = 1-1/m)
  minus the inner-root term N^{-1} E[ V(Lam+xi)^{1/m} ]   (V = smaller root; power-small).
Consequences:  D_n = N^{-1}[F(Lam+k1) + (k2/2) F''(Lam+k1) + O(Lam^{1/m-3})]
               D_n = tau_n^+(2/e) (1 + O(1/log^2 n)),   tau_n^+(s) = outer point (distance from x=1) where (f^n)' = s
               D_n ~ (log n/(m beta n))^{1/m}.
For f = sin(pi x/2): D_n ~ (2/pi) sqrt(log n / n).

Test functions: 2x - x^2 (m=2, beta=1), sin(pi x/2) (m=2, beta=pi^2/8), 1-(1-x)^3 (m=3, beta=1).
"""
import mpmath as mp
import sympy as sp
from fnlib import g, gprime

xs = sp.Symbol('x', positive=True)
CASES = {
    '2x-x^2':     (2*xs - xs**2, 2, sp.Integer(1)),
    'sin(pi x/2)': (sp.sin(sp.pi*xs/2), 2, sp.pi**2/8),
    '1-(1-x)^3':  (1 - (1 - xs)**3, 3, sp.Integer(1)),
}
K1 = 1 - mp.log(2)
K2 = 1 + mp.pi**2/6
K3 = 2 - mp.mpf(3)/2*mp.zeta(3)

def Uroot(y, c):
    """larger root of U - c log U = y (y >= ymin):  U = -c W_{-1}(-e^{-y/c}/c)"""
    return -c*mp.lambertw(-mp.exp(-y/c)/c, -1).real

def Vroot(y, c):
    """smaller root of V - c log V = y :  V = -c W_0(-e^{-y/c}/c)"""
    return -c*mp.lambertw(-mp.exp(-y/c)/c, 0).real

def Fderivs(y, m):
    """F = U^{1/m} and its derivatives up to order 3 at y, via U' = U/(U - c)."""
    c = 1 - mp.mpf(1)/m
    U = Uroot(y, c)
    # use mpmath diff on the implicit function
    F = lambda yy: Uroot(yy, c)**(mp.mpf(1)/m)
    return [F(y)] + [mp.diff(F, y, k) for k in (1, 2, 3)]

def defect_deg(f, fp, n, m, beta):
    n = mp.mpf(n)
    N = (beta*n)**(mp.mpf(1)/m)
    nm1 = n - 1
    def integrand(v):
        t = v/N
        xx = 1 - t
        fx = f(xx)
        if fx <= 0:
            return mp.mpf(0)
        s = n*mp.exp(nm1*mp.log(fx))*fp(xx)
        return g(s)
    Lam = mp.log(m*N)
    vmax = (Lam + (mp.mp.dps + 20)*mp.log(10))**(mp.mpf(1)/m)*mp.mpf('1.2')
    vmax = min(vmax, N)
    inner = N**(-mp.mpf(1)/(m - 1))          # inner boundary layer where slope ~ 1 near x = 1
    pts = [mp.mpf(0)]
    p = inner*mp.mpf(10)**-3
    while p < mp.mpf('0.5'):
        pts.append(p); p *= 4
    pts += [mp.mpf('0.5'), mp.mpf(1)]
    c0 = Lam**(mp.mpf(1)/m)
    for k in [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 4, 8]:
        q = c0 + mp.mpf(k)/max(c0**(m - 1), 1)
        if 1 < q < vmax:
            pts.append(q)
    pts.append(vmax)
    pts = sorted(set(pts))
    val = mp.quad(integrand, pts)/N
    return val, N, Lam

def tau_plus(f, fp, n, s, m, beta):
    """outer solution tau of (f^n)'(1-tau) = s"""
    n = mp.mpf(n)
    F = lambda t: (n - 1)*mp.log(f(1 - t)) + mp.log(fp(1 - t)) + mp.log(n) - mp.log(s)
    N = (beta*n)**(mp.mpf(1)/m)
    t0 = (mp.log(m*N) + mp.log(1/s))**(mp.mpf(1)/m)/N
    return mp.findroot(F, t0)

def EH(Hfun, Lam, ymin):
    """E[H(Lam + xi)], xi = log(1/S), restricted to Lam + xi >= ymin"""
    lo = ymin - Lam
    w = lambda xi: gprime(mp.exp(-xi))*mp.exp(-xi)
    pts = [lo, lo + mp.mpf('0.01'), lo + 1, lo + 4]
    pts = [p for p in pts if p < 0] + [mp.mpf(0), mp.mpf(2), mp.mpf(10), mp.mpf(40), mp.inf]
    pts = sorted(set(pts))
    return mp.quad(lambda xi: w(xi)*Hfun(Lam + xi), pts)

if __name__ == '__main__':
    for name, (expr, m, beta_s) in CASES.items():
        mp.mp.dps = 60
        f = sp.lambdify(xs, expr, 'mpmath'); fp = sp.lambdify(xs, sp.diff(expr, xs), 'mpmath')
        beta = mp.mpf(sp.N(beta_s, 80))
        c = 1 - mp.mpf(1)/m
        ymin = c*(1 - mp.log(c))
        print('=' * 110)
        print('f = %s   (m=%d, beta=%s)   leading law D_n ~ (log n/(m beta n))^(1/m)' % (name, m, mp.nstr(beta, 12)))
        print('%-7s %-24s %-11s %-11s %-11s %-11s %-11s' % ('n', 'D_n', 'lead', 'tau+(2/e)', 'F+k2F"/2', 'E[F]', 'E[F]-E[V]'))
        for e in [4, 6, 10, 15, 20, 30, 40]:
            n = mp.mpf(10)**e
            mp.mp.dps = 50 + 2*e
            f = sp.lambdify(xs, expr, 'mpmath'); fp = sp.lambdify(xs, sp.diff(expr, xs), 'mpmath')
            beta = mp.mpf(sp.N(beta_s, mp.mp.dps + 10))
            c = 1 - mp.mpf(1)/m; ymin = c*(1 - mp.log(c))
            D, N, Lam = defect_deg(f, fp, n, m, beta)
            lead = (mp.log(n)/(m*beta*n))**(mp.mpf(1)/m)
            tp = tau_plus(f, fp, n, 2/mp.e, m, beta)
            Fd = Fderivs(Lam + K1, m)
            two = (Fd[0] + K2/2*Fd[2])/N
            EF = EH(lambda y: Uroot(y, c)**(mp.mpf(1)/m), Lam, ymin)/N
            EV = EH(lambda y: Vroot(y, c)**(mp.mpf(1)/m), Lam, ymin)/N
            full = EF - EV
            rel = lambda P: mp.nstr(P/D - 1, 3)
            print('1e%-5d %-24s %-11s %-11s %-11s %-11s %-11s' % (e, mp.nstr(D, 18), rel(lead), rel(tp), rel(two), rel(EF), rel(full)))
        print('  columns: relative error (prediction/D_n - 1).  last column = full pure-power model; expected')
        print('  residual O(log n / n) (Phi even in t for these three cases).')
