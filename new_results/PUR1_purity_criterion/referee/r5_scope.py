# Referee check of (g): (i) a further linear image of B_p; (iii) a polygon soft-max whose gamma_12
# coefficient cancels although H_j != 1; (ii) jump sign at the split of the designed rectangle.
import sympy as sp
import mpmath as mp

out = []
P_ = lambda *a: (out.append(' '.join(str(t) for t in a)), print(*a))

NMAX = 6
x, mu, eps, v = sp.symbols('x mu eps v')
Z = {k: sp.Symbol('z%d' % k) for k in range(2, NMAX + 3)}
logQ = sum(Z[k] / k * (((mu + x) / 2) ** k + ((mu - x) / 2) ** k - mu ** k) * eps ** k for k in range(2, NMAX + 2))
Qs = sp.expand(1 - sp.series(sp.exp(logQ), eps, 0, NMAX + 2).removeO())
R = 0
for d in range(2, NMAX + 2):
    qt, rm = sp.div(sp.Poly(sp.expand(Qs.coeff(eps, d)), x, mu), sp.Poly(mu ** 2 - x ** 2, x, mu))
    assert rm.is_zero
    R += qt.as_expr()
tau = sp.Symbol('tau')


def B3(psi_even):
    tot = 0
    for k2, ps in psi_even.items():
        if ps != 0:
            tot += ps * tau ** (2 + k2) * sp.diff(R, x, k2).subs(x, 0).subs(mu, tau)
    tot = sp.expand(tot)
    return {n: sp.simplify(tot.coeff(tau, n)) for n in range(NMAX + 1)}


znum = lambda e: e.subs({Z[k]: sp.zeta(k) for k in Z})

# ---------------- (i) K = A B_p with A = [[1,1],[0,2]] ----------------
A = sp.Matrix([[1, 1], [0, 2]])
detA = A.det()
AinvT = (A.T).inv()
xi1, xi2 = sp.symbols('xi1 xi2')
Qf = sp.expand(((AinvT * sp.Matrix([xi1, xi2])).T * (AinvT * sp.Matrix([xi1, xi2])))[0])
P_('A =', A.tolist(), ' Q(xi) =', Qf)
coef = {}
psis = {}
for sgn in (1, -1):
    q = lambda vv: 2 / (detA * Qf.subs({xi1: sp.exp(vv), xi2: sgn * sp.exp(-vv)}) ** sp.Rational(3, 2))
    psi = q(v) + q(-v)
    ser = sp.series(psi, v, 0, NMAX + 1).removeO()
    pc = {k: sp.nsimplify(sp.simplify(ser.coeff(v, k))) for k in range(NMAX + 1)}
    psis[sgn] = psi
    P_('kink sign %+d: psi_k =' % sgn, pc)
    assert all(pc[k] == 0 for k in range(1, NMAX + 1, 2))
    coef[sgn] = B3({k: pc[k] for k in range(0, NMAX + 1, 2)})
cn = {n: sp.simplify(2 * coef[1][n] + 2 * coef[-1][n]) for n in range(NMAX + 1)}
for n in range(NMAX + 1):
    P_('  c_%d =' % n, cn[n])

mp.mp.dps = 60
Am = mp.matrix([[1, 1], [0, 2]])


def h_lin(th, tq):
    xx = Am[0, 0] * mp.cos(th) + Am[1, 0] * mp.sin(th)
    yy = Am[0, 1] * mp.cos(th) + Am[1, 1] * mp.sin(th)
    M = max(abs(xx), abs(yy))
    r = min(abs(xx), abs(yy)) / M
    return M * (1 + r ** (1 / tq)) ** tq


# kink angles where |x| = |y|
kinks = []
for th0 in [mp.mpf(k) * mp.pi / 8 for k in range(16)]:
    pass
fx = lambda th: (mp.cos(th) + mp.sin(th)) ** 2 - (mp.cos(th) + 2 * mp.sin(th)) ** 2 + 0  # placeholder
# x = cos, y = cos + 2 sin : |x|=|y| <=> sin = 0 or cos + sin = 0
kinks = [mp.mpf(0), mp.pi, 3 * mp.pi / 4, 7 * mp.pi / 4]
# zeros of x or y (|t| = 0: h is only C^{floor(1/tau)} there, so they must be breakpoints too)
zeros_xy = [mp.pi / 2, 3 * mp.pi / 2, mp.atan2(-1, 2) + mp.pi, mp.atan2(-1, 2) + 2 * mp.pi]


def per_theta(tq):
    pts = sorted(set([mp.mpf(0), 2 * mp.pi] + [k + s * tq * c for k in kinks + [2 * mp.pi]
                                                 for s in (-1, 1) for c in (0.5, 2, 8, 30)
                                                 if 0 < k + s * tq * c < 2 * mp.pi] + kinks + zeros_xy))
    return mp.quad(lambda th: h_lin(th, tq), pts)


fpsi = {s: sp.lambdify(v, psis[s], 'mpmath') for s in (1, -1)}


def per_v(tq):
    return 2 * mp.fsum(mp.quad(lambda vv: mp.exp(vv) * fpsi[s](vv) * (1 + mp.exp(-2 * vv / tq)) ** tq,
                               [0, tq, 4 * tq, 1, 5, 30, mp.inf]) for s in (1, -1))


for tq in [mp.mpf('0.5'), mp.mpf('0.2'), mp.mpf('0.07')]:
    a, b = per_theta(tq), per_v(tq)
    P_('  tau=%s: per by theta-quadrature %s ; v-form - theta = %s' % (tq, mp.nstr(a, 25), mp.nstr(b - a, 3)))

cnum = {n: mp.mpf(sp.N(znum(cn[n]), 70)) for n in range(NMAX + 1)}
mp.mp.dps = 70
cnum[0] = 2 * mp.fsum(mp.quad(lambda vv: mp.exp(vv) * fpsi[s](vv), [0, 1, 5, 30, mp.inf]) for s in (1, -1))
P_('  c_0 = per(A B_1) =', mp.nstr(cnum[0], 30), '; minus polygon perimeter 4+4sqrt2:',
   mp.nstr(cnum[0] - 4 - 4 * mp.sqrt(2), 3))


def neville0(xs, ys):
    ys = list(ys)
    for m in range(1, len(xs)):
        for i in range(len(xs) - m):
            ys[i] = (xs[i] * ys[i + 1] - xs[i + m] * ys[i]) / (xs[i] - xs[i + m])
    return ys[0]


mp.mp.dps = 70
taus = [mp.mpf(1) / (20 * 2 ** i) for i in range(11)]
vals = [per_theta(tq) for tq in taus]
for target in (4, 5, 6):
    g = [(Pv - mp.fsum(cnum[n] * tq ** n for n in range(target))) / tq ** target for tq, Pv in zip(taus, vals)]
    est = neville0(taus, g)
    P_('  c_%d: Neville from theta-quadrature %s ; criterion %s ; diff %s' % (
        target, mp.nstr(est, 25), mp.nstr(cnum[target], 25), mp.nstr(est - cnum[target], 3)))

# ---------------- (iii) soft-max of the rectangle [-a,a]x[-b,b], a=6/5, b=3/5 ----------------
mp.mp.dps = 60
a_, b_ = mp.mpf(6) / 5, mp.mpf(3) / 5
P_('rectangle soft-max a=6/5, b=3/5: sum_j (1-H_j)/T_j over kinks =', (1 - a_) / b_ * 2 + (1 - b_) / a_ * 2)


def h_sm(th, tq):
    n1, n2 = mp.cos(th), mp.sin(th)
    vs = [sx * a_ * n1 + sy * b_ * n2 for sx in (-1, 1) for sy in (-1, 1)]
    m = max(vs)
    return m * mp.fsum(mp.exp((w - m) / tq) for w in vs) ** tq


def per_sm(tq):
    ks = [k * mp.pi / 2 for k in range(5)]
    pts = sorted(set([k + s * tq * c for k in ks for s in (-1, 1) for c in (0.5, 2, 8, 30)
                      if 0 < k + s * tq * c < 2 * mp.pi] + ks))
    return mp.quad(lambda th: h_sm(th, tq), pts)


# predicted coefficients by route (a): kink x: phi = 2a/b + 2v/sqrt(b^2-v^2); kink y: a<->b; 2 kinks each
def gq(a, b):
    return mp.quad(lambda s: s ** a * mp.log1p(mp.exp(-2 * s)) ** b, [0, 1, 4, 12, 40, mp.inf]) / (
        mp.factorial(a) * mp.factorial(b))


GAM = {(a, b): gq(a, b) for a in range(7) for b in range(1, 7)}
vv = sp.Symbol('vv')
phis = [sp.Rational(2) * sp.Rational(6, 5) / sp.Rational(3, 5) + 2 * vv / sp.sqrt(sp.Rational(9, 25) - vv ** 2),
        sp.Rational(2) * sp.Rational(3, 5) / sp.Rational(6, 5) + 2 * vv / sp.sqrt(sp.Rational(36, 25) - vv ** 2)]
pc = []
for ph in phis:
    s_ = sp.series(ph, vv, 0, 8).removeO()
    pc.append({k: mp.mpf(sp.N(s_.coeff(vv, k), 70)) for k in range(8)})
    ps = sp.series(sp.exp(-vv) * ph, vv, 0, 6).removeO()
    P_('   psi coefficients of this kink:', [ps.coeff(vv, k) for k in range(6)])
cpred = {0: 4 * (a_ + b_), 1: mp.mpf(0)}
for n in range(2, 7):
    cpred[n] = 2 * mp.fsum(p[k] * mp.factorial(k) * GAM[(k, n - 1 - k)] for p in pc for k in range(n - 1))
P_('   predicted c_4 =', mp.nstr(cpred[4], 30), '; (5/16)zeta(4) =', mp.nstr(mp.mpf(5) / 16 * mp.zeta(4), 30))
P_('   predicted c_2 - (5/2)z2 =', mp.nstr(cpred[2] - mp.mpf(5) / 2 * mp.zeta(2), 3),
   '; c_3 - (5/2)z3 =', mp.nstr(cpred[3] - mp.mpf(5) / 2 * mp.zeta(3), 3))
mp.mp.dps = 60
taus = [mp.mpf(1) / (100 * 2 ** i) for i in range(9)]
vals = [per_sm(tq) for tq in taus]
for target in (4, 5):
    g = [(Pv - mp.fsum(cpred[n] * tq ** n for n in range(target))) / tq ** target for tq, Pv in zip(taus, vals)]
    est = neville0(taus, g)
    P_('   c_%d: Neville from direct quadrature %s ; predicted %s ; diff %s' % (
        target, mp.nstr(est, 25), mp.nstr(cpred[target], 25), mp.nstr(est - cpred[target], 3)))
# what c_4 would be with the gamma_12 term of one kink type only (to show the test is sensitive)
P_('   size of one kink-pair gamma_12 term 2*psi_1*gamma_12 =', mp.nstr(2 * mp.mpf(2) / 3 * GAM[(1, 2)], 10))

# ---------------- (ii) designed rectangle a=1, b=2: jump sign at the split ----------------
mp.mp.dps = 40
A1, B1 = mp.mpf(1), mp.mpf(2)
L1 = lambda th: mp.log(1 + (B1 / A1) * mp.sin(th) + (B1 / A1) ** 2 * (1 - mp.cos(th)))
L2 = lambda th: mp.log(1 + (A1 / B1) * mp.sin(mp.pi / 2 - th) + (A1 / B1) ** 2 * (1 - mp.cos(mp.pi / 2 - th)))
ths = mp.findroot(lambda th: L1(th) - L2(th), 0.8)
for tq in [mp.mpf('0.9'), mp.mpf('0.3'), mp.mpf('0.1')]:
    hfun = lambda th, L: (A1 * mp.cos(th) + B1 * mp.sin(th)) * (1 + mp.exp(-2 * L(th) / tq)) ** tq
    dl = mp.diff(lambda th: hfun(th, L1), ths)
    dr = mp.diff(lambda th: hfun(th, L2), ths)
    P_('designed rectangle 1x2, tau=%s: split at %s, jump h\'(+)-h\'(-) = %s' % (
        tq, mp.nstr(ths, 10), mp.nstr(dr - dl, 6)))
open(__file__.replace('.py', '_output.txt'), 'w').write('\n'.join(out) + '\n')
