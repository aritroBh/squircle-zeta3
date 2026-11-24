"""
verify_lemmas.py -- 40-digit checks of the exact ingredients.

 (L1) layer-cake / random-slope identity  2 - length(y) = int_0^inf g'(s) |{x: y'(x) > s}| ds
      (convex cases f^n, n in {3, 7.5, 20}, for all test functions; plus a NON-monotone slope case
       y = (2x - x^2)^2, whose superlevel sets {y' > s} are intervals (x_-(s), x_+(s)))
 (L2) Mellin transform  int_0^inf g'(s) s^w ds = R(w)  at real and complex w in (-1,1)
 (L3) moments/cumulants of xi = log(1/S): E xi = 1 - log 2,  Var xi = 1 + pi^2/6,
      third cumulant 2 - 3 zeta(3)/2,  geometric mean of S = 2/e
 (L4) Mellin-Barnes formula  D_n = (1/2 pi i) int_{(c)} R(-z) Z_n(z) dz/z,  Z_n(z) = int_0^1 (y_n')^z dx
"""
import mpmath as mp
import sympy as sp
from fnlib import TESTS, make_callables, g, gprime, R, logR_coeffs, moments_xi

mp.mp.dps = 50
TOL = mp.mpf(10)**-40
worst = mp.mpf(0)

def report(label, val, ref):
    global worst
    err = abs(val - ref)
    rel = err/abs(ref) if ref != 0 else err
    worst = max(worst, rel)
    print('%-58s %s   rel.resid = %s' % (label, mp.nstr(val, 25), mp.nstr(rel, 3)))

print('=== (L2) Mellin transform of g\' ===')
for w in [mp.mpf('-0.7'), mp.mpf('-0.3'), mp.mpf('0.2'), mp.mpf('0.5'), mp.mpf('0.9'), mp.mpc('0.3', '2.0')]:
    # substitute s = e^u ; integrand g'(e^u) e^{u(w+1)}
    val = mp.quad(lambda u: gprime(mp.exp(u))*mp.exp(u*(w + 1)), [-mp.inf, -5, 0, 5, mp.inf])
    report('int g\'(s) s^w ds, w=%s' % mp.nstr(w, 3), val, R(w))

print('=== (L3) cumulants of xi = log(1/S) ===')
def Exi(k):
    return mp.quad(lambda u: gprime(mp.exp(u))*mp.exp(u)*(-u)**k, [-mp.inf, -5, 0, 5, mp.inf])
m1, m2, m3 = Exi(1), Exi(2), Exi(3)
report('E[xi] = 1 - log 2', m1, 1 - mp.log(2))
report('Var[xi] = 1 + pi^2/6', m2 - m1**2, 1 + mp.pi**2/6)
report('k3[xi] = 2 - 3 zeta(3)/2', m3 - 3*m2*m1 + 2*m1**3, 2 - 3*mp.zeta(3)/2)
report('geometric mean of S = 2/e', mp.exp(-m1), 2/mp.e)
rho = moments_xi(4)
report('rho_2 from exact log R series', rho[2], m2)
report('rho_3 from exact log R series', rho[3], m3)
report('int_0^inf g\'(s) ds = 1 (probability density)', mp.quad(gprime, [0, 1, mp.inf]), mp.mpf(1))

print('=== (L1) layer-cake identity ===')
xs = sp.Symbol('x', positive=True)
for name, expr in TESTS.items():
    f, fp = make_callables(expr)
    fpp = sp.lambdify(xs, sp.diff(expr, xs, 2), 'mpmath')
    for n in [mp.mpf(3), mp.mpf('7.5'), mp.mpf(20)]:
        yp = lambda t: n*f(t)**(n - 1)*fp(t)
        ypp = lambda t: n*((n - 1)*f(t)**(n - 2)*fp(t)**2 + f(t)**(n - 1)*fpp(t))
        Ddirect = mp.quad(lambda t: g(yp(t)), mp.linspace(0, 1, 9))
        smax = yp(mp.mpf(1))
        def m(s):
            if s >= smax:
                return mp.mpf(0)
            # y' increasing (checked: f^n convex for these n); solve y'(X)=s by bisection+Newton
            lo, hi = mp.mpf(0), mp.mpf(1)
            for _ in range(60):
                mid = (lo + hi)/2
                if yp(mid) > s:
                    hi = mid
                else:
                    lo = mid
            X = mp.findroot(lambda t: yp(t) - s, (lo + hi)/2)
            return 1 - X
        Dlayer = mp.quad(lambda s: gprime(s)*m(s), [0, smax/8, smax/2, smax])
        report('%-7s n=%-4s  layer-cake vs direct' % (name, mp.nstr(n, 3)), Dlayer, Ddirect)
# non-monotone slope: y=(2x-x^2)^2, y'=4x(2-x)(1-x)*... compute generally
y = (2*xs - xs**2)**2
ypf = sp.lambdify(xs, sp.diff(y, xs), 'mpmath')
Ddirect = mp.quad(lambda t: g(ypf(t)), [0, mp.mpf(1)/3, mp.mpf(2)/3, 1])
xstar = mp.findroot(lambda t: sp.lambdify(xs, sp.diff(y, xs, 2), 'mpmath')(t), 0.4)
smax = ypf(xstar)
def bis(F, lo, hi, it=180):
    # F(lo) < 0 < F(hi) or reverse; plain bisection to full precision
    flo = F(lo)
    for _ in range(it):
        mid = (lo + hi)/2
        fm = F(mid)
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi)/2
def m2fun(s):
    if s >= smax:
        return mp.mpf(0)
    xm = bis(lambda t: ypf(t) - s, mp.mpf(0), xstar)
    xp = bis(lambda t: ypf(t) - s, xstar, mp.mpf(1))
    return xp - xm
Dlayer = mp.quad(lambda s: gprime(s)*m2fun(s), [0, smax/2, smax])
report('y=(2x-x^2)^2 (slope not monotone) layer-cake vs direct', Dlayer, Ddirect)

print('=== (L4) Mellin inversion g(s) = (1/2 pi i) int_(c) R(-z) s^z dz/z  (=> Mellin-Barnes formula for D_n by Fubini) ===')
for sv in [mp.mpf('0.3'), mp.mpf(1), mp.mpf(7)]:
    c = mp.mpf('0.5')
    val = mp.quad(lambda tt: (R(-mp.mpc(c, tt))*sv**mp.mpc(c, tt)/mp.mpc(c, tt)).real, mp.linspace(-90, 90, 31))/(2*mp.pi)
    report('Mellin inversion at s=%s' % mp.nstr(sv, 3), val, g(sv))
print()
print('WORST relative residual over all checks:', mp.nstr(worst, 3))
