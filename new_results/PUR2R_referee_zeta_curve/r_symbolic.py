"""Referee PUR2R: symbolic checks of (A) coefficient table, (C) values, (D) via an
independent route (moment expansion F_n = sum_{k+b=n-1} phi_k k! gamma_{kb}, gamma by quadrature)."""
import sympy as sp
import mpmath as mp

t, e = sp.symbols('tau e')
z = {k: sp.Symbol('z%d' % k) for k in range(2, 9)}
N = 8
# log R = sum (p_k - 1) zeta(k) tau^k / k with p_k = a^k + b^k, a+b=1, ab=e (Newton)
p = {0: sp.Integer(2), 1: sp.Integer(1)}
for k in range(2, N + 1):
    p[k] = sp.expand(p[k - 1] - e * p[k - 2])
print("p_k - 1:", {k: sp.factor(p[k] - 1) for k in range(2, 8)})
logR = sum((p[k] - 1) * z[k] * t**k / k for k in range(2, N + 1))
R = sp.series(sp.exp(logR), t, 0, N + 1).removeO()
P = sp.expand(8 + sp.expand((2 / e) * (1 - R)))
red = {z[2]**2: sp.Rational(5, 2) * z[4], z[2] * z[4]: sp.Rational(7, 4) * z[6], z[2]**3: sp.Rational(35, 8) * z[6]}
def reduce_(x):
    x = sp.expand(x)
    x = sp.expand(x.subs(z[2]**3, sp.Rational(35, 8) * z[6]))
    x = sp.expand(x.subs(z[2]**2, sp.Rational(5, 2) * z[4]))
    x = sp.expand(x.subs(z[2] * z[4], sp.Rational(7, 4) * z[6]))
    return x
coef = {n: reduce_(P.coeff(t, n)) for n in range(2, 8)}
for n in range(2, 8):
    print("tau^%d:" % n, sp.collect(coef[n], [z[k] for k in range(2, 8)]))
claimed = {
    4: (2 - sp.Rational(7, 2) * e) * z[4],
    5: 2 * z[5] - 2 * e * (z[5] + z[2] * z[3]),
    6: (2 - sp.Rational(13, 2) * e + sp.Rational(31, 8) * e**2) * z[6] - e * z[3]**2,
    7: 2 * z[7] - e * (2 * z[2] * z[5] + 2 * z[3] * z[4] + 4 * z[7]) + e**2 * (z[2]**2 * z[3] + 2 * z[2] * z[5] + z[3] * z[4] + 2 * z[7]),
}
for n in claimed:
    print("  claimed table tau^%d matches:" % n, sp.simplify(reduce_(claimed[n]) - coef[n]) == 0)
# (C) at e = 4/7
for n in (4, 5, 6):
    print("e=4/7 tau^%d:" % n, sp.collect(sp.expand(coef[n].subs(e, sp.Rational(4, 7))), [z[k] for k in range(2, 8)]))
mp.mp.dps = 40
zv = {z[k]: mp.zeta(k) for k in range(2, 9)}
for n in (5, 6):
    print("  numeric tau^%d at e=4/7:" % n, mp.nstr(sp.N(coef[n].subs(e, sp.Rational(4, 7)).subs(zv), 30), 25))
# e = 0 -> 2 zeta(n)
P0 = sp.limit((2 / e) * (1 - sp.exp(logR)), e, 0)
print("e->0 series:", sp.expand(sp.series(P0, t, 0, 8).removeO()))

# (D) independent route: c_n = 8 F_n, F_n = sum_{k+b=n-1, b>=1} phi_k k! gamma_{kb},
# gamma_{kb} = 1/(k! b!) int_0^inf s^k log(1+e^{-2s})^b ds, phi = e^v psi, psi = 1 + s2 v^2/2 + s4 v^4/24 + s6 v^6/720
s2, s4, s6 = sp.symbols('psi2 psi4 psi6')   # psi''(0), psi''''(0), psi^(6)(0)
v = sp.Symbol('v')
psi = 1 + s2 * v**2 / 2 + s4 * v**4 / 24 + s6 * v**6 / 720
phis = sp.expand(sp.series(sp.exp(v) * psi, v, 0, 7).removeO())
phik = [phis.coeff(v, k) for k in range(7)]
mp.mp.dps = 45
def gam(k, b):
    return mp.quad(lambda s: s**k * mp.log1p(mp.exp(-2 * s))**b, [0, 1, 4, 16, mp.inf]) / (mp.factorial(k) * mp.factorial(b))
G = {(k, b): gam(k, b) for k in range(0, 6) for b in range(1, 7) if k + b <= 7}
# sanity: gamma_{k1} = eta(k+2)/2^{k+1}
print("gamma_k1 vs eta(k+2)/2^(k+1):", max(abs(G[(k, 1)] - mp.altzeta(k + 2) / 2**(k + 1)) for k in range(6)))
X = sp.Symbol('X')
def cn(n):
    out = 0
    for k in range(0, n - 1):
        b = n - 1 - k
        out += phik[k] * sp.factorial(k) * sp.Float(mp.nstr(G[(k, b)], 44), 44)
    return sp.expand(8 * out)
Z = {k: mp.zeta(k) for k in range(2, 9)}
claimD = {
    2: 2 * Z[2] + 0 * s2,
    3: 2 * Z[3] + 0 * s2,
    4: (9 + 7 * s2) * Z[4] / 8,
    5: (3 + s2) * Z[5] / 2 + (s2 - 1) * Z[2] * Z[3] / 2,
    6: (sp.Rational(79, 128) + sp.Rational(73, 64) * s2 + sp.Rational(31, 128) * s4) * Z[6] + (s2 - 1) * Z[3]**2 / 4,
    7: (sp.Rational(9, 8) + sp.Rational(3, 4) * s2 + s4 / 8) * Z[7] + (-sp.Rational(3, 8) + s2 / 4 + s4 / 8) * Z[2] * Z[5]
       + (-sp.Rational(9, 32) + s2 / 16 + sp.Rational(7, 32) * s4) * Z[3] * Z[4],
}
for n in range(2, 8):
    d = sp.expand(cn(n) - sp.expand(claimD[n]))
    cs = [abs(sp.N(d.subs({s2: 0, s4: 0, s6: 0}))), abs(sp.N(d.coeff(s2))), abs(sp.N(d.coeff(s4))), abs(sp.N(d.coeff(s6)))]
    print("c_%d: |independent - claimed| coefficientwise [1, psi'', psi'''', psi6]:" % n, [mp.nstr(mp.mpf(str(x)), 3) for x in cs])
# coefficient of psi^(2k)(0) in c_{2k+2}: 4^{1-k} eta(2k+2)
for k in (1, 2):
    sym = [s2, s4][k - 1]
    print("coef of psi^(%d) in c_%d:" % (2 * k, 2 * k + 2), mp.nstr(mp.mpf(str(sp.N(cn(2 * k + 2).coeff(sym), 40))), 30), " 4^{1-k}eta:", mp.nstr(mp.mpf(4)**(1 - k) * mp.altzeta(2 * k + 2), 30))
# gamma_22 implied by the claimed c5
print("gamma_22 =", mp.nstr(G[(2, 2)], 30), " implied:", mp.nstr(-mp.mpf(29) / 256 * Z[5] + Z[2] * Z[3] / 16, 30))
