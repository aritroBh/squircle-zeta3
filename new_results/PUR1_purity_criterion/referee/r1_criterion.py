# Referee check of (b), (c), (d): symbolic criterion (B3) vs route (a) with quadrature gamma_ab.
# Own code; does not import the author's scripts.
import sympy as sp
import mpmath as mp

NMAX = 7
x, mu, eps, v = sp.symbols('x mu eps v')
Z = {k: sp.Symbol('z%d' % k) for k in range(2, NMAX + 3)}   # z_k stands for zeta(k)

# log Q, Q = G(1-(mu+x)/2) G(1-(mu-x)/2) / G(1-mu), with x,mu -> eps*x, eps*mu
logQ = sum(Z[k] / k * (((mu + x) / 2) ** k + ((mu - x) / 2) ** k - mu ** k) * eps ** k
           for k in range(2, NMAX + 2))
Q = sp.series(sp.exp(logQ), eps, 0, NMAX + 2).removeO()
one_minus_Q = sp.expand(1 - Q)
# exact division by (mu^2 - x^2) degree by degree in eps
R = 0
for d in range(2, NMAX + 2):
    num = sp.expand(one_minus_Q.coeff(eps, d))
    qt, rm = sp.div(sp.Poly(num, x, mu), sp.Poly(mu ** 2 - x ** 2, x, mu))
    assert rm.is_zero, ('not divisible at degree', d)
    R += qt.as_expr() * eps ** (d - 2)
R = sp.expand(R)
# check: Q-1 vanishes on x=+-mu (independent of division)
assert sp.expand(one_minus_Q.subs(x, mu)) == 0 and sp.expand(one_minus_Q.subs(x, -mu)) == 0


def coeffs_B3(psi_even):
    """psi_even: dict k-> psi_{k} (even k). Returns dict n -> F_n via (B3):
       F ~ tau^2 sum_k psi_{2k} tau^{2k} d_x^{2k} R(x,tau)|_{x=0}."""
    tau = sp.Symbol('tau')
    tot = 0
    for k2, ps in psi_even.items():
        if ps == 0:
            continue
        dR = sp.diff(R, x, k2).subs(x, 0).subs(eps, 1)
        tot += ps * tau ** (2 + k2) * dR.subs(mu, tau)
    tot = sp.expand(tot)
    return {n: sp.simplify(tot.coeff(tau, n)) for n in range(0, NMAX + 1)}


def zeta_to_num(e):
    return e.subs({Z[k]: sp.zeta(k) for k in Z})


def weight_check(e, n):
    # every monomial in z's must have total weight n
    e = sp.expand(e)
    if e == 0:
        return True
    for term in sp.Add.make_args(e):
        w = 0
        for k, z in Z.items():
            w += k * sp.degree(term, z)
        if w != n:
            return False
    return True


out = []
P = lambda *a: (out.append(' '.join(str(t) for t in a)), print(*a))

# ---------------- squircle ----------------
psi_sq = sp.series(1 / (sp.sqrt(2) * sp.cosh(2 * v) ** sp.Rational(3, 2)), v, 0, NMAX + 1).removeO()
psi_sq_c = {k: sp.nsimplify(psi_sq.coeff(v, k)) for k in range(0, NMAX + 1)}
P('squircle psi_k:', psi_sq_c)
Fsq = coeffs_B3({k: psi_sq_c[k] for k in range(0, NMAX + 1, 2)})
claimed = {
    2: sp.sqrt(2) * Z[2] / 8,
    3: sp.sqrt(2) * Z[3] / 8,
    4: -sp.Rational(33, 16) * sp.sqrt(2) * Z[4] / 8,
    5: -sp.sqrt(2) * (sp.Rational(3, 4) * Z[5] + sp.Rational(7, 4) * Z[2] * Z[3]) / 8,
    6: sp.sqrt(2) * (sp.Rational(4039, 256) * Z[6] - sp.Rational(7, 8) * Z[3] ** 2) / 8,
    7: sp.sqrt(2) * (sp.Rational(129, 16) * Z[7] + sp.Rational(141, 16) * Z[2] * Z[5]
                     + sp.Rational(1071, 64) * Z[3] * Z[4]) / 8,
}
for n in range(0, NMAX + 1):
    cn = sp.expand(8 * Fsq[n])
    P('squircle c_%d (8F_n, from B3) =' % n, cn, '  weight-homogeneous:', weight_check(cn / sp.sqrt(2), n))
    if n in claimed:
        d = sp.expand(8 * Fsq[n] - 8 * claimed[n])
        # compare as symbols first; if not 0, compare after zeta(2k)->pi^{2k}
        dn = sp.nsimplify(sp.simplify(zeta_to_num(d)))
        P('   minus claimed (symbolic, z-basis):', d, '| after zeta(2k)->pi:', dn)
        P('   numeric value c_%d =' % n, sp.N(zeta_to_num(cn), 40))

# ---------------- route (a) with quadrature gamma_ab ----------------
mp.mp.dps = 45


def gamma_ab(a, b):
    f = lambda s: s ** a * mp.log1p(mp.exp(-2 * s)) ** b
    return mp.quad(f, [0, 1, 4, 12, mp.inf]) / (mp.factorial(a) * mp.factorial(b))


GAM = {(a, b): gamma_ab(a, b) for a in range(0, NMAX) for b in range(1, NMAX)}


def coeffs_a(phi_c):
    return {n: mp.fsum(phi_c[k] * mp.factorial(k) * GAM[(k, n - 1 - k)]
                       for k in range(0, n - 1)) for n in range(2, NMAX + 1)}


phi_sq = sp.series(sp.exp(v) / (sp.sqrt(2) * sp.cosh(2 * v) ** sp.Rational(3, 2)), v, 0, NMAX + 1).removeO()
phi_sq_c = {k: mp.mpf(sp.N(phi_sq.coeff(v, k), 60)) for k in range(0, NMAX + 1)}
Fa = coeffs_a(phi_sq_c)
for n in range(2, NMAX + 1):
    b3 = mp.mpf(sp.N(zeta_to_num(Fsq[n]), 60))
    P('squircle F_%d: route(a) - B3 = %s' % (n, mp.nstr(Fa[n] - b3, 5)))

# ---------------- soft-max ----------------
phi_sm = 1 + v / sp.sqrt(1 - v ** 2)
psi_sm = sp.series(sp.exp(-v) * phi_sm, v, 0, NMAX + 1).removeO()
psi_sm_c = {k: psi_sm.coeff(v, k) for k in range(0, NMAX + 1)}
P('soft-max psi_k:', psi_sm_c)
Fsm_even = coeffs_B3({k: psi_sm_c[k] for k in range(0, NMAX + 1, 2)})
phi_sm_s = sp.series(phi_sm, v, 0, NMAX + 1).removeO()
Fa_sm = coeffs_a({k: mp.mpf(sp.N(phi_sm_s.coeff(v, k), 60)) for k in range(0, NMAX + 1)})


def dkG_tau_tau_coeffs(k, nmax):
    # coefficients of tau^{k+1} d^k G(tau,tau) = sum gamma_ab a!/(a-k)! tau^{a-k+b+k+1}
    c = {}
    for a in range(k, nmax):
        for b in range(1, nmax):
            n = a + b + 1
            if n <= nmax:
                c[n] = c.get(n, 0) + GAM[(a, b)] * mp.factorial(a) / mp.factorial(a - k)
    return c


odd = {n: mp.mpf(0) for n in range(0, NMAX + 1)}
for k in range(1, NMAX + 1, 2):
    pk = mp.mpf(sp.N(psi_sm_c[k], 60))
    if pk == 0:
        continue
    for n, val in dkG_tau_tau_coeffs(k, NMAX).items():
        odd[n] += pk * val
z = lambda k: mp.zeta(k)
for n in range(2, NMAX + 1):
    tot = 8 * (mp.mpf(sp.N(zeta_to_num(Fsm_even[n]), 60)) + odd[n])
    P('soft-max c_%d: B1(even via B3 + odd via gammas) = %s ; route(a) = %s ; diff %s' % (
        n, mp.nstr(tot, 30), mp.nstr(8 * Fa_sm[n], 30), mp.nstr(tot - 8 * Fa_sm[n], 3)))
P('soft-max even part 8F_n symbolic:', {n: sp.expand(8 * Fsm_even[n]) for n in range(2, NMAX + 1)})
c6_claim = 40 * GAM[(3, 2)] - z(3) ** 2 / 2 + mp.mpf(11) / 16 * z(6)
P('soft-max c6 claimed 40g32 - z3^2/2 + 11/16 z6 =', mp.nstr(c6_claim, 35), ' minus route(a):',
  mp.nstr(c6_claim - 8 * Fa_sm[6], 3))
c5_claim = mp.mpf(107) / 32 * z(5) - z(2) * z(3)
P('soft-max c5 claimed minus route(a):', mp.nstr(c5_claim - 8 * Fa_sm[5], 3))
P('gamma_31 - 15 z5/256 =', mp.nstr(GAM[(3, 1)] - 15 * z(5) / 256, 3),
  '; gamma_41 - eta(6)/32 =', mp.nstr(GAM[(4, 1)] - (1 - mp.mpf(2) ** -5) * z(6) / 32, 3))
open(__file__.replace('.py', '_output.txt'), 'w').write('\n'.join(out) + '\n')
