"""PUR2 symbolic part: (A) coefficients of the family as polynomials in e = ab and zeta values,
(B) c=1 series, (C) tau^4 zero, (D) universal coefficients via R(x,tau) = (1-Q)/(tau^2-x^2)."""
import sympy as sp

t, x, e, c = sp.symbols('tau x e c')
N = 8
Z = {k: sp.Symbol('z%d' % k) for k in range(2, N + 2)}   # formal zeta(k), no relations used
zsub = {Z[k]: sp.zeta(k) for k in Z}

def p_k(k):  # a^k + b^k as polynomial in e=ab with a+b=1 (Newton recursion p_k = p_{k-1} - e p_{k-2})
    p = [sp.Integer(2), sp.Integer(1)]
    for j in range(2, k + 1):
        p.append(sp.expand(p[-1] - e * p[-2]))
    return p[k]

print("=== (A) a^k+b^k-1 in e=ab (each divisible by e) ===")
for k in range(2, N + 1):
    q = sp.expand(p_k(k) - 1)
    print("k=%d: %s ; /e = %s" % (k, q, sp.factor(q / e)))

logR = sum((p_k(k) - 1) * Z[k] * t**k / k for k in range(2, N + 1))
R = sp.expand(sp.series(sp.exp(logR), t, 0, N + 1).removeO())
perim = sp.expand(8 + sp.expand(2 * (1 - R) / e))   # perimeter minus tail
print("\n=== (A) perimeter = 8 + (2/ab)(1 - R) : coefficients (formal zeta symbols) ===")
coef = {}
for n in range(0, N + 1):
    cn = sp.expand(perim.coeff(t, n))
    coef[n] = cn
    print("tau^%d: %s" % (n, sp.collect(cn, e)))
    assert sp.simplify(cn - sp.expand(cn)) == 0
# reduce even zetas to compare tau^4 claim
c4 = sp.nsimplify(sp.expand(coef[4].subs(Z[2]**2, sp.Rational(5, 2) * Z[4])))
print("\ntau^4 after zeta(2)^2 = 5/2 zeta(4):", sp.factor(c4), "  claim (2 - 7e/2) zeta(4):",
      sp.simplify(c4 - (2 - sp.Rational(7, 2) * e) * Z[4]) == 0)
print("zero of tau^4 at e =", sp.solve(2 - sp.Rational(7, 2) * e, e), "; c^2 = 1-4e =", 1 - 4 * sp.Rational(4, 7))
c5 = sp.expand(coef[5])
print("tau^5:", sp.collect(c5, e))
print("tau^6 (zeta(2)^3, zeta(2)zeta(4) reduced to zeta(6)):",
      sp.collect(sp.expand(coef[6].subs({Z[2]**3: sp.Rational(35, 8) * Z[6]}).subs(Z[2] * Z[4], sp.Rational(7, 4) * Z[6])), e))

print("\n=== (B) e = 0 (c = 1): perimeter = 8 + 2 sum zeta(k) tau^k ===")
for n in range(0, N + 1):
    print("tau^%d: %s" % (n, sp.expand(coef[n].subs(e, 0))))
# check -tau(gamma+psi(1-tau)) series independently
ser = sp.series(-2 * t * (sp.EulerGamma + sp.polygamma(0, 1 - t)), t, 0, 8).removeO()
print("series of -2 tau (gamma + psi(1-tau)):", sp.expand(ser))
print("matches 2 sum zeta(k) tau^k through tau^7:",
      sp.simplify(sp.expand(ser) - sum(2 * sp.zeta(k) * t**k for k in range(2, 8))) == 0)

print("\n=== (D) universal coefficients: phi = e^v psi, psi even, psi(0)=1 ===")
# log Q(x,tau) = sum zeta(k)/k [((tau+x)/2)^k + ((tau-x)/2)^k - tau^k]; R = (1-Q)/(tau^2-x^2)
s = sp.Symbol('s')   # homogeneity bookkeeping: x->s x, tau->s tau
logQ = sum(Z[k] / k * (((t + x) / 2)**k + ((t - x) / 2)**k - t**k) * s**k for k in range(2, N + 1))
Q = sp.expand(sp.series(sp.exp(logQ), s, 0, N + 1).removeO())
num = sp.expand(1 - Q)
# divide by tau^2 - x^2 degree by degree (exact polynomial division in x,tau)
Rx = 0
for d in range(2, N + 1):
    part = sp.expand(num.coeff(s, d))
    qd, rd = sp.div(sp.Poly(part, x, t), sp.Poly(t**2 - x**2, x, t))
    assert rd.is_zero, "not divisible at degree %d" % d
    Rx += qd.as_expr()
Rx = sp.expand(Rx)
P = sp.symbols('P2 P4 P6')   # psi''(0), psi''''(0), psi^(6)(0)
# F ~ tau^2 sum_k psi^{(2k)}(0) tau^{2k} [x^{2k}] R ; perimeter = 8 + 8F
pd = {0: 1, 2: P[0], 4: P[1], 6: P[2]}
F = 0
for k2, pk in pd.items():
    F += t**2 * pk * t**k2 * sp.Poly(Rx, x).as_expr().coeff(x, k2) if k2 > 0 else t**2 * Rx.subs(x, 0)
F = sp.expand(F)
univ = {}
for n in range(2, N + 1):
    cn = sp.expand(8 * F.coeff(t, n))
    univ[n] = cn
    red = sp.expand(cn.subs({Z[2]**3: sp.Rational(35, 8) * Z[6]}).subs(Z[2]**2, sp.Rational(5, 2) * Z[4]).subs(Z[2] * Z[4], sp.Rational(7, 4) * Z[6]))
    print("tau^%d: %s\n        reduced: %s" % (n, sp.collect(cn, P), sp.collect(red, P)))

print("\nconsistency: psi = cosh(c v) (P2=c^2, P4=c^4, P6=c^6) vs family with e=(1-c^2)/4:")
for n in range(2, N + 1):
    d = sp.expand(univ[n].subs({P[0]: c**2, P[1]: c**4, P[2]: c**6}) - coef[n].subs(e, (1 - c**2) / 4))
    print("  tau^%d difference: %s" % (n, sp.simplify(d)))

print("\nconsistency with PUR1 squircle: psi = 1/cosh^{3/2}(2v) times sqrt2/2")
vv = sp.Symbol('v')
psq = sp.series(1 / sp.cosh(2 * vv)**sp.Rational(3, 2), vv, 0, 8).removeO()
P2s, P4s, P6s = [sp.diff(psq, vv, k).subs(vv, 0) for k in (2, 4, 6)]
print("  psi''(0), psi''''(0), psi^(6)(0) =", P2s, P4s, P6s)
known = {4: -sp.Rational(33, 16) * sp.sqrt(2) * sp.zeta(4),
         5: -sp.sqrt(2) * (sp.Rational(3, 4) * sp.zeta(5) + sp.Rational(7, 4) * sp.zeta(2) * sp.zeta(3)),
         6: sp.sqrt(2) * (sp.Rational(4039, 256) * sp.zeta(6) - sp.Rational(7, 8) * sp.zeta(3)**2),
         7: sp.sqrt(2) * (sp.Rational(129, 16) * sp.zeta(7) + sp.Rational(141, 16) * sp.zeta(2) * sp.zeta(5) + sp.Rational(1071, 64) * sp.zeta(3) * sp.zeta(4))}
for n, kv in known.items():
    mine = sp.sqrt(2) / 2 * univ[n].subs({P[0]: P2s, P[1]: P4s, P[2]: P6s}).subs(zsub)
    print("  tau^%d: (D)-formula - PUR1 value = %s" % (n, sp.N(mine - kv, 40)))
