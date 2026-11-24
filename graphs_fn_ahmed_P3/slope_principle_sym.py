"""
slope_principle_sym.py -- the "slope 2/e" form of the expansion (generalizes Ahmed's Cor. 2(b) and Prop. A).

  x_n^* := the point where (f^n)'(x) = 2/e,  tau_* = 1 - x_n^* = Theta(mu Xi, mu)   (exactly).
  D_n - tau_* = c2 mu^2 + mu^3 (c31 Xi + c30) + mu^4 (c42 Xi^2 + c41 Xi + c40) + ...
  i.e.  L_n = 1 + x_n^* - (D_n - tau_*).
Checks: f = x reproduces Ahmed's Proposition A:  L(n) = 1 + (2/(en))^{1/(n-1)} + pi^2/(12 n^2) + (d - a2 log n)/n^3 + ...
"""
import sympy as sp
import symbolic_Qk as S

mu, Xi, W, h = S.mu, S.Xi, S.W, S.h
a, b, c, d = S.a, S.b, S.c, S.d
K = S.K
tau_star = sp.expand(S.T.subs(h, 1).subs(W, mu*Xi))
tau_star = sum(tau_star.coeff(mu, k)*mu**k for k in range(0, K + 1))
diff = sp.expand(S.D - tau_star)
out = {}
for k in range(1, K + 1):
    out[k] = sp.collect(sp.expand(diff.coeff(mu, k)), Xi)
    print('[mu^%d] (D - tau_*) =' % k)
    for (m,), coeff in sorted(sp.Poly(out[k], Xi).terms(), reverse=True):
        print('     Xi^%d :' % m, sp.factor(coeff))

# numeric constants
lnR_vals = {S.k2: 1 + sp.pi**2/6, S.k3: 2 - sp.Rational(3, 2)*sp.zeta(3), S.k4: 6 + sp.pi**4/15}
# ---- f = x, compare with Ahmed Prop. A (in powers of 1/n): L = 1 + x* + pi^2/(12 n^2) + (dA - a2 log n)/n^3
L, nn = sp.symbols('L n', positive=True)
kap = 1 - sp.log(2)
a2 = sp.pi**2/12 + sp.Rational(1, 2)
dA = a2*(1 + sp.log(2)) + sp.zeta(3)/4 - sp.Rational(4, 3)
corr = -(out[2]*mu**2 + out[3]*mu**3)      # L - 1 - x*  through mu^3
corr = corr.subs({a: 1, b: 0, c: 0, d: 0}).subs(lnR_vals).subs(Xi, L + kap).subs(mu, 1/(nn - 1))
ser = sp.series(corr, nn, sp.oo, 4).removeO()
ahmed = sp.pi**2/(12*nn**2) + (dA - a2*L)/nn**3
print()
print('f = x:  (ours) - (Ahmed Prop. A) through 1/n^3 =', sp.simplify(sp.expand(ser - ahmed)))
print('c2(f) closed form check:', sp.simplify(-out[2] - ((a**2 - b)*sp.pi**2/12 + (a**2 - a - b)/2)/a**3).subs(lnR_vals) if False else
      sp.simplify(sp.expand(-out[2].subs(lnR_vals) - ((a**2 - b)*sp.pi**2/12 + (a**2 - a - b)/2)/a**3)))
print()
for k in (2, 3):
    print('LaTeX [mu^%d] (L_n - 1 - x_n^*) =' % k, sp.latex(sp.collect(sp.expand(-out[k]), Xi)))
