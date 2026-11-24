"""(d) soft-max decomposition, (e) first-failure normalisation and alternating sums,
   test example for (e) at j=0, p-norm rhombus coefficients for (g)."""
import pickle
import sympy as sp, mpmath as mp
from pur_core import *
from pur_core import _gcache

mp.mp.dps = 70
N = 7
Sser = S_series(N)


def eta_sym(n):
    return (1 - sp.Rational(1, 2**(n - 1))) * Z[n]


def reduce_g1(expr):
    """gamma_{a,1} = eta(a+2)/2^{a+1} (elementary)."""
    rep = {gam(a, 1): eta_sym(a + 2) / 2**(a + 1) for a in range(0, 12)}
    return sp.expand(expr.subs(rep))


def split(psi, N):
    pt = taylor(psi, v, N)
    ev = {k: pt[2 * k] for k in range(N // 2 + 1)}
    od = {k: pt[k] for k in range(1, N + 1, 2) if pt[k] != 0}
    return pt, ev, od


print("=== (d) soft-max: phi = 1 + v/sqrt(1-v^2), psi = e^{-v} phi ===")
phi = 1 + v / sp.sqrt(1 - v**2)
psi = sp.exp(-v) * phi
pt, ev, od = split(psi, N)
print("psi Taylor coefficients:", [pt[k] for k in range(N + 1)])
print("odd part:", od)
Fe = even_part_coeffs(ev, N, Sser)
Fo = odd_part_coeffs(od, N)
Fa = direct_coeffs(taylor(phi, v, N), N)
for n in range(2, N + 1):
    cn = reduce_g1(8 * (Fe[n] + Fo[n]))
    print("c_%d = 8F_%d = %s" % (n, n, sp.simplify(to_pi(cn))))
    print("      numeric: (b)+odd = %s ; (a) - that = %s" % (mp.nstr(num(cn), 40), mp.nstr(num(8 * Fa[n]) - num(cn), 3)))
c6 = reduce_g1(8 * (Fe[6] + Fo[6]))
sm6 = sp.Rational(11, 16) * Z[6] - Z[3]**2 / 2 + 40 * gam(3, 2)
print("c6 - SM6 form (11/16 z6 - z3^2/2 + 40 g32), symbolic in pi:", sp.simplify(to_pi(c6 - sm6)))
c7 = reduce_g1(8 * (Fe[7] + Fo[7]))
print("c7 from criterion:", c7)
sm6_c7 = 40 * gam(3, 3) + sp.Rational(3, 4) * Z[3] * Z[4] + sp.Rational(2339, 512) * Z[7]
d7 = sp.expand(c7 - sm6_c7)
print("c7 - SM6 form (40 g33 + 3/4 z3 z4 + 2339/512 z7) =", d7, " numeric:", mp.nstr(num(d7), 5))
print("  -> requires the gamma_42 term to be pure; see PSLQ below")
sm6_num = {6: mp.mpf('0.04295082788028403392943380483513526551029672769206'),
           7: mp.mpf('5.585741214319102958232826658315468061731504678')}
print("c6 numeric %s ; minus SM6 Richardson value = %s" % (mp.nstr(num(c6), 45), mp.nstr(num(c6) - sm6_num[6], 3)))
print("c7 numeric %s ; minus SM6 Richardson value = %s" % (mp.nstr(num(c7), 45), mp.nstr(num(c7) - sm6_num[7], 3)))

print("\n=== (e) gamma_{a,2} = 2^{-(a+1)} sum_{n>=2} (-1)^n H_{n-1}/n^{a+2} ===")
mp.mp.dps = 80


def altH(s):
    # sum_{n>=1} (-1)^n H_{n-1}/n^s ;  H_{n-1} = harmonic(n-1)
    return mp.nsum(lambda n: (-1)**int(n) * mp.harmonic(n - 1) / n**s, [2, mp.inf], method='alternating')


for a in range(0, 7):
    g = gamma_num(a, 2)
    s = altH(a + 2) / 2**(a + 1)
    print("a=%d  gamma_a2 = %s   series - quad = %s" % (a, mp.nstr(g, 40), mp.nstr(s - g, 3)))

print("\nClassical: A := sum_{n>=1} (-1)^{n+1} H_n/n^3 = 11/4 z4 - 7/4 z3 log2 + z2/2 log^2 2 - log^4 2/12 - 2 Li4(1/2)")
mp.mp.dps = 100
A = mp.nsum(lambda n: (-1)**(int(n) + 1) * mp.harmonic(n) / n**3, [1, mp.inf], method='alternating')
L2 = mp.log(2)
Acl = mp.mpf(11) / 4 * mp.zeta(4) - mp.mpf(7) / 4 * mp.zeta(3) * L2 + mp.zeta(2) / 2 * L2**2 - L2**4 / 12 - 2 * mp.polylog(4, mp.mpf(1) / 2)
print("A = %s ; A - closed form = %s" % (mp.nstr(A, 50), mp.nstr(A - Acl, 3)))
# gamma_12 = (1/4) sum (-1)^n H_{n-1}/n^3 = (1/4)(eta(4) - A)
_gcache.clear()
g12 = gamma_num(1, 2)
g12cl = (mp.mpf(7) / 8 * mp.zeta(4) - Acl) / 4
print("gamma_12 = %s ; minus (eta(4)-A)/4 = %s" % (mp.nstr(g12, 50), mp.nstr(g12 - g12cl, 3)))
print("gamma_12 = -15/32 z4 + 7/16 z3 log2 - 1/8 z2 log^2 2 + log^4 2/48 + 1/2 Li4(1/2):",
      mp.nstr(g12 - (-mp.mpf(15) / 32 * mp.zeta(4) + mp.mpf(7) / 16 * mp.zeta(3) * L2 - mp.zeta(2) / 8 * L2**2 + L2**4 / 48 + mp.polylog(4, 0.5) / 2), 3))

print("\nPSLQ (dps 100, maxcoeff 1e12): gamma_{a,2} against pure-zeta bases of weight a+3")
z = mp.zeta
bases = {
    1: (['z4'], [z(4)]),
    2: (['z5', 'z2z3'], [z(5), z(2) * z(3)]),
    3: (['z6', 'z3^2'], [z(6), z(3)**2]),
    4: (['z7', 'z2z5', 'z3z4'], [z(7), z(2) * z(5), z(3) * z(4)]),
    5: (['z8', 'z3z5', 'z2z3^2'], [z(8), z(3) * z(5), z(2) * z(3)**2]),
}
for a, (names, vals) in bases.items():
    _gcache.clear()
    g = gamma_num(a, 2)
    rel = mp.pslq([g] + vals, maxcoeff=10**12, maxsteps=10**6)
    print("  a=%d basis %s : %s" % (a, names, rel))
_gcache.clear()
g12 = gamma_num(1, 2)
rel = mp.pslq([g12, z(4), z(3) * L2, z(2) * L2**2, L2**4, mp.polylog(4, 0.5)], maxcoeff=10**12, maxsteps=10**6)
print("  gamma_12 vs {z4, z3 log2, z2 log^2 2, log^4 2, Li4(1/2)}:", rel)

print("\n=== test example for (e), j=0: psi = (1+v) e^{-v^2}, omega = 1 ===")
mp.mp.dps = 70
_gcache.clear()
psi = (1 + v) * sp.exp(-v**2)
pt, ev, od = split(psi, 5)
Fe = even_part_coeffs(ev, 5, S_series(5))
Fo = odd_part_coeffs(od, 5)
Fa = direct_coeffs(taylor(sp.exp(v) * psi, v, 5), 5)
known = {}
for n in range(0, 6):
    e = reduce_g1(Fe[n] + Fo[n])
    known[n] = str(num(e)) if e != 0 else '0'
    print("F_%d = %s   [(a) minus this: %s]" % (n, sp.simplify(to_pi(e)), mp.nstr((num(Fa[n]) if Fa[n] != 0 else 0) - (num(e) if e != 0 else 0), 3)))
print("  F_4 decomposes as pure + 1!*gamma_12 + 2!*gamma_21 (omega=1, j=0)")
pickle.dump(known, open('test_e_known.pkl', 'wb'))

print("\n=== (g) p-norm rhombus diag(a,b), psi = q(v)+q(-v), q = 2/(ab (e^{2v}/a^2 + e^{-2v}/b^2)^{3/2}) ===")
a_, b_ = sp.Integer(1), sp.Integer(2)
qf = lambda w: 2 / (a_ * b_ * (sp.exp(2 * w) / a_**2 + sp.exp(-2 * w) / b_**2)**sp.Rational(3, 2))
psi = qf(v) + qf(-v)
pt, ev, od = split(psi, 6)
print("odd part of psi (should be empty):", od)
Fe = even_part_coeffs(ev, 6, S_series(6))
known = {}
r = sp.sqrt(a_**2 + b_**2)
for n in range(0, 7):
    known[n] = str(num(Fe[n])) if Fe[n] != 0 else '0'
    if n >= 2:
        print("c_%d = 4F_%d = %s" % (n, n, sp.simplify(4 * Fe[n])))
print("paper thm:rhombus check: c2 = 4 a^2 b^2 zeta(2)/r^3 ->", sp.simplify(4 * Fe[2] - 4 * a_**2 * b_**2 * Z[2] / r**3),
      "; c3 = 4 a^2 b^2 zeta(3)/r^3 ->", sp.simplify(4 * Fe[3] - 4 * a_**2 * b_**2 * Z[3] / r**3))
pickle.dump(known, open('rhombus_known.pkl', 'wb'))
