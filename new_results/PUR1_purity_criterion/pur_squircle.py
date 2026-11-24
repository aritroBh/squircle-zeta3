"""(b)+(c): squircle via the purity criterion; check (b) against (a) with numeric gamma_ab."""
import sympy as sp, mpmath as mp
from pur_core import *

mp.mp.dps = 70
N = 7
print("S(nu,y) series up to weight", N)
Sser = S_series(N)
# check S numerically against Gamma definition at a point
nu0, y0 = mp.mpf('0.031'), mp.mpf('0.017')
lam, mu = 2 * nu0 + 2 * y0, 2 * nu0
Sex = (mp.gamma(-lam / 2) * mp.gamma(lam / 2 - mu) / mp.gamma(-mu) + 2 / lam + 2 / (2 * mu - lam)) / 2
Sap = num(Sser.subs({nu: nu0, y: y0}).subs({}))  # polynomial in z's
print("S series vs Gamma form at (nu,y)=(0.031,0.017): diff = %s (truncation O(0.05^8) expected)" % mp.nstr(Sex - Sap, 5))

# ---- squircle: phi(v) = e^v/(sqrt2 cosh^{3/2} 2v), psi = 1/(sqrt2 cosh^{3/2}(2v))
psi = 1 / (sp.sqrt(2) * sp.cosh(2 * v)**sp.Rational(3, 2))
phi = sp.exp(v) * psi
pt = taylor(psi, v, N)
ph = taylor(phi, v, N)
print("psi Taylor:", [sp.simplify(pt[k]) for k in range(N + 1)])
assert all(pt[k] == 0 for k in range(1, N + 1, 2))
psi_even = {k: pt[2 * k] for k in range(0, N // 2 + 1)}
Fb = even_part_coeffs(psi_even, N, Sser)
Fa = direct_coeffs(ph, N)
known = {
    0: 4 * sp.sqrt(2), 1: 0, 2: sp.sqrt(2) * Z[2], 3: sp.sqrt(2) * Z[3], 4: -sp.Rational(33, 16) * sp.sqrt(2) * Z[4],
    5: -sp.sqrt(2) * (sp.Rational(3, 4) * Z[5] + sp.Rational(7, 4) * Z[2] * Z[3]),
    6: sp.sqrt(2) * (sp.Rational(4039, 256) * Z[6] - sp.Rational(7, 8) * Z[3]**2),
    7: sp.sqrt(2) * (sp.Rational(129, 16) * Z[7] + sp.Rational(141, 16) * Z[2] * Z[5] + sp.Rational(1071, 64) * Z[3] * Z[4]),
}
print("\nP = 4 sqrt2 + 8 F.  Coefficients c_n = 8 F_n from (b):")
for n in range(1, N + 1):
    cb = sp.expand(8 * Fb[n])
    diff = sp.simplify(to_pi(cb - known[n]))
    na = num(8 * Fa[n]) if Fa[n] != 0 else mp.mpf(0)
    nb = num(cb) if cb != 0 else mp.mpf(0)
    print("n=%d  c_n(b) = %s" % (n, sp.factor(cb)))
    print("      c_n(b) in pi: %s" % sp.simplify(to_pi(cb)))
    print("      exact minus known (symbolic, zeta(2k)->pi^2k): %s" % diff)
    print("      numeric c_n(b) = %s ;  (a) with quadrature gammas minus (b) = %s" % (mp.nstr(nb, 40), mp.nstr(na - nb, 5)))
import pickle
pickle.dump({n: sp.srepr(sp.expand(8 * Fb[n])) for n in range(N + 1)}, open('squircle_cn.pkl', 'wb'))
