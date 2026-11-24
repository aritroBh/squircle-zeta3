"""m2_structure.py -- EXACT (sympy, c symbolic) derivation of the M2(alpha) integrand and of the
weight split M2 = W3 + W2.

Corner function (corner_series.py):
  M(u;a) = int_0^1 [ (1+g)^(-1-u) sqrt(1+g^(2-2u)-2 cos(a) g^(1-u)) - 1 ] g^(u-1) dg = sum_k M_k u^k.
Steps (each printed as an exact sympy check):
 [1] u^2-coefficient of the integrand, derived independently with L = log g, m = log(1+g) as symbols:
     integrand_2 = (1/g)[ RA L^2 + RB L m + RC m^2 ]  (no other monomials).
 [2] Euler substitution t = g + sqrt(1-2ca g+g^2): g = (t^2-1)/(2(t-ca)), sqrt = P/(2(t-ca)), P = t^2-2ca t+1;
     t-integrand coefficients RA2, RB2, RC2 (rational in t over Q(c)).
 [3] KEY IDENTITY: RA2 = RA2s - s^2 d/dt[(t+1)/P], with RA2s having only simple poles at the real letters
     {1,-1,ca,r1,r2};  s^2 = 1-c^2.
 [4] Integration by parts with V = 1/2 - s^2 (t+1)/P = (t-1)(t-ca)/(2P) (V(1)=0, log g(t1)=0) and
     V * d/dt log g = 1/(2(t+1)) EXACTLY  ==>  int RA2c log^2 g dt = - int log g dt/(t+1).
 Result:  M2 = W3 + W2,  W2 = - int_1^{t1} log g dt/(t+1),
          W3 = int_1^{t1} [RA2s log^2 g + RB2 log g log(1+g) + RC2 log^2(1+g)] dt  (real letters only).
 [5] residues of RA2s, RB2, RC2 at the letters.
"""
import sympy as sp

t, c, g, L, m, u = sp.symbols('t c g L m u')
ca = 2*c**2 - 1
s2 = 1 - c**2                      # s^2, s = sin(a/2)
r1, r2 = 2*c - 1, -1 - 2*c

# ---------- [1] u-expansion with L, m symbols --------------------------------------------------
# (1+g)^(-1-u) = exp(-u m)/(1+g);  g^(-u) = exp(-u L);  g^(u-1) = exp(u L)/g
S0 = 1 + g**2 - 2*ca*g
SQ = sp.Symbol('SQ', positive=True)          # SQ = sqrt(S0)
Q = 1 + g**2*sp.exp(-2*u*L) - 2*ca*g*sp.exp(-u*L)
d = sp.series((Q - S0)/S0, u, 0, 3).removeO()
sqrtQ_ser = SQ*(1 + d/2 - d**2/8)            # (1+d)^(1/2) to O(u^3), d = O(u)
integ = sp.exp(-u*m)/(1 + g)*sqrtQ_ser*sp.exp(u*L)/g - sp.exp(u*L)/g
c2 = sp.expand(sp.series(integ, u, 0, 3).removeO().coeff(u, 2))
poly = sp.Poly(c2, L, m)
print('[1] monomials in (L,m) of the u^2 coefficient:', sorted(poly.monoms()))
RA_g = sp.simplify(poly.coeff_monomial(L**2))
RB_g = sp.simplify(poly.coeff_monomial(L*m))
RC_g = sp.simplify(poly.coeff_monomial(m**2))

# ---------- [2] Euler substitution -------------------------------------------------------------
P = t**2 - 2*ca*t + 1
g_t = (t**2 - 1)/(2*(t - ca))
sqrtS0_t = P/(2*(t - ca))
print('[2] (sqrtS0_t)^2 - S0(g_t) =', sp.simplify(sqrtS0_t**2 - S0.subs(g, g_t)))
print('    g_t + sqrtS0_t - t =', sp.simplify(g_t + sqrtS0_t - t))
dg_dt = sp.diff(g_t, t)


def to_t(expr):
    e = expr.subs(SQ, sqrtS0_t).subs(g, g_t)
    return sp.factor(sp.cancel(e*dg_dt))


RA2, RB2, RC2 = to_t(RA_g), to_t(RB_g), to_t(RC_g)
print('    RB2 =', RB2)
print('    RC2 =', RC2)
print('    1 + g_t =', sp.factor(1 + g_t))

# ---------- [3] key identity --------------------------------------------------------------------
RA2c = -s2*sp.diff((t + 1)/P, t)
RA2s = sp.factor(sp.cancel(RA2 - RA2c))
den = sp.factor(sp.denom(sp.together(RA2s)))
print('[3] RA2s = RA2 + s^2 d/dt[(t+1)/P] has denominator', den)
print('    P divides denominator of RA2s?', sp.rem(sp.Poly(sp.denom(sp.together(RA2s)), t),
                                               sp.Poly(P, t)).is_zero)

# ---------- [4] IBP identity ---------------------------------------------------------------------
V = sp.Rational(1, 2) - s2*(t + 1)/P
print('[4] V - (t-1)(t-ca)/(2P) =', sp.simplify(V - (t - 1)*(t - ca)/(2*P)))
dlogg = sp.diff(sp.log(g_t), t)
print('    V*dlog g/dt - 1/(2(t+1)) =', sp.simplify(V*dlogg - 1/(2*(t + 1))))
print('    dV/dt - RA2c =', sp.simplify(sp.diff(V, t) - RA2c))
print('    g(t1) with t1 = 1+2s (s^2=1-c^2):',
      sp.simplify((((1 + 2*sp.sqrt(s2))**2 - 1)/(2*(1 + 2*sp.sqrt(s2) - ca)))))

# ---------- [5] residues ------------------------------------------------------------------------
letters = {'1': 1, '-1': -1, 'ca': ca, 'r1': r1, 'r2': r2}
for name, R in [('RA2s', RA2s), ('RB2', RB2), ('RC2', RC2)]:
    res = {}
    for k, rho in letters.items():
        res[k] = sp.factor(sp.simplify(sp.limit((t - rho)*R, t, rho)))
    # check: R - sum res/(t-rho) == polynomial part (should be 0)
    rest = sp.simplify(R - sum(res[k]/(t - letters[k]) for k in letters))
    print(f'[5] {name}: residues', res, ' remainder =', rest)

import pickle
with open('m2_structure.pkl', 'wb') as f:
    pickle.dump({'RA2': RA2, 'RB2': RB2, 'RC2': RC2, 'RA2s': RA2s}, f)
print('saved m2_structure.pkl')
