# Independent numeric cross-check: the REAL S_m(t) polynomial (genintegrand.S_poly) vs the proved N_l/Pi_l
import sympy as sp, time
from genintegrand import S_poly
import tow2lib as T
d, m = T.d, T.m; i = sp.I
t = sp.Symbol('t'); ts = (1+i)/2
t0 = time.time(); bad = 0; n = 0
for dd in range(3, 10):
    for mm in range(0, dd+1):
        P = S_poly(dd, mm)
        poly = sum(sp.Rational(c.numerator, c.denominator)*t**k for k, c in enumerate(P))
        S0 = sp.expand(poly.subs(t, ts))
        for l in range(1, 5):
            lhs = sp.expand(sp.diff(poly, t, l).subs(t, ts))/S0
            rhs = (T.N[l]/T.Pi(l)).subs({d: dd, m: mm})
            diff = sp.simplify(sp.expand(lhs - rhs))
            n += 1
            if diff != 0: bad += 1; print('MISMATCH', dd, mm, l, diff)
print(f'S^(l)(t*)/S(t*) vs N_l/Pi_l using real S_poly: {n} cells (d=3..9, all m, l=1..4), {bad} bad  ({time.time()-t0:.1f}s)')
