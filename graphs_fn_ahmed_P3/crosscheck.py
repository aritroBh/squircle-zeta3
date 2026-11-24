"""
crosscheck.py
 (a) independent check of the defect routine: direct quadrature of 2 - int_0^1 sqrt(1 + y'^2) dx
     (high working precision to absorb cancellation) vs fnlib.defect, at n = 10^3 and 10^5, all test functions;
 (b) locality: f2 = x + 0.02 sin^2(pi (x-0.1)/0.3) on [0.1,0.4] (f2 = x elsewhere) has the same germ at x=1
     as f1 = x; Proposition 2 predicts |D_n(f2) - D_n(f1)| <= f2(1/2)^n + (1/2)^n = 2^{1-n}.
"""
import mpmath as mp
import sympy as sp
from fnlib import TESTS, make_callables, defect, g

xs = sp.Symbol('x', positive=True)
print('(a) direct sqrt-quadrature vs defect routine')
worst = 0
for name, expr in TESTS.items():
    f, fp = make_callables(expr)
    a = mp.mpf(sp.N(sp.diff(expr, xs).subs(xs, 1), 120))
    for n in [10**3, 10**5]:
        mp.mp.dps = 110
        nn = mp.mpf(n)
        pts = sorted(set([mp.mpf(0), mp.mpf(1)] + [1 - mp.mpf(k)/(nn*a) for k in [0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256]
                                                   if 1 - mp.mpf(k)/(nn*a) > 0]))
        Ldirect = mp.quad(lambda t: mp.sqrt(1 + (nn*f(t)**(nn - 1)*fp(t))**2), pts)
        Dd = 2 - Ldirect
        mp.mp.dps = 80
        Dr, _ = defect(f, fp, n, a, dps_target=70)
        rel = abs(Dd - Dr)/Dr
        worst = max(worst, rel)
        print('   %-7s n=%-7d D_direct=%s   rel.diff=%s' % (name, n, mp.nstr(Dd, 30), mp.nstr(rel, 3)))
print('   worst relative difference:', mp.nstr(worst, 3))

print('(b) locality: D_n(f2) - D_n(x), f2 = x + bump on [0.1,0.4]')
mp.mp.dps = 80
def f2(t):
    if mp.mpf('0.1') < t < mp.mpf('0.4'):
        return t + mp.mpf('0.02')*mp.sin(mp.pi*(t - mp.mpf('0.1'))/mp.mpf('0.3'))**2
    return t
def f2p(t):
    if mp.mpf('0.1') < t < mp.mpf('0.4'):
        u = mp.pi*(t - mp.mpf('0.1'))/mp.mpf('0.3')
        return 1 + mp.mpf('0.02')*2*mp.sin(u)*mp.cos(u)*mp.pi/mp.mpf('0.3')
    return mp.mpf(1)
for n in [10, 20, 40, 80, 160]:
    nn = mp.mpf(n)
    pts = [mp.mpf(0), mp.mpf('0.1'), mp.mpf('0.25'), mp.mpf('0.4'), mp.mpf('0.7'), mp.mpf('0.9'), 1 - 1/nn, mp.mpf(1)]
    pts = sorted(set(pts))
    D2 = mp.quad(lambda t: g(nn*f2(t)**(nn - 1)*f2p(t)), pts)
    D1 = mp.quad(lambda t: g(nn*t**(nn - 1)), pts)
    print('   n=%-4d  D_n(x)=%s  |D_n(f2)-D_n(x)| = %-10s  bound 2^(1-n) = %-10s  0.42^n = %s' % (
        n, mp.nstr(D1, 15), mp.nstr(abs(D2 - D1), 4), mp.nstr(mp.mpf(2)**(1 - n), 4), mp.nstr(mp.mpf('0.42')**n, 4)))
