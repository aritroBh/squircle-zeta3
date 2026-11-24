# Referee check of (f): designed square. Closed form vs full-circle quadrature, tail series,
# h+h'' by sympy differentiation (not the author's formula) at dps 200, jumps, and the tau-expansion.
import mpmath as mp
import sympy as sp

out = []
P_ = lambda *a: (out.append(' '.join(str(t) for t in a)), print(*a))


def h_theta(th, tau):
    # u = signed angle from the nearest edge normal j*pi/2
    j = mp.nint(th / (mp.pi / 2))
    u = th - j * mp.pi / 2
    au = abs(u)
    m = mp.cos(u) + mp.sin(au)
    return m * (1 + (2 + mp.sin(au) - mp.cos(u)) ** (-2 / tau)) ** tau


def closed(tau, nterms=60):
    g = 0 if tau == 1 else mp.gamma(1 - tau / 2) ** 2 / mp.gamma(1 - tau)
    s = mp.fsum(mp.binomial(tau, n) * mp.mpf(4) ** (-n / tau) / (2 * n - tau) for n in range(1, nterms))
    return 16 - 8 * g - 16 * tau * s


mp.mp.dps = 60
for tau in [mp.mpf(1), mp.mpf('0.7'), mp.mpf('0.37'), mp.mpf('0.1'), mp.mpf('0.04')]:
    pts = sorted(set([k * mp.pi / 4 for k in range(9)] +
                     [j * mp.pi / 2 + sg * tau * c for j in range(5) for sg in (-1, 1) for c in (0.5, 2, 8)
                      if 0 < j * mp.pi / 2 + sg * tau * c < 2 * mp.pi]))
    I = mp.quad(lambda t: h_theta(t, tau), pts)
    P_('tau=%s  int_0^{2pi} h = %s  closed form - quadrature = %s' % (
        mp.nstr(tau, 4), mp.nstr(I, 40), mp.nstr(closed(tau) - I, 3)))

# tail T(V) series vs quadrature
tau = mp.mpf('0.6')
V = mp.log(2)
Tq = mp.quad(lambda v: mp.exp(v) * ((1 + mp.exp(-2 * v / tau)) ** tau - 1), [V, 2, 5, 20, mp.inf])
Ts = tau * mp.exp(V) * mp.fsum(mp.binomial(tau, n) * mp.exp(-2 * n * V / tau) / (2 * n - tau) for n in range(1, 200))
P_('T(log2) at tau=0.6: series - quadrature =', mp.nstr(Ts - Tq, 3))

# expansion of 16 - 8 Gamma(1-t/2)^2/Gamma(1-t)
tc = mp.taylor(lambda t: 16 - 8 * mp.gamma(1 - t / 2) ** 2 / mp.gamma(1 - t), 0, 7)
z = mp.zeta
claim = [8, 0, 2 * z(2), 2 * z(3), mp.mpf(9) / 8 * z(4), mp.mpf(3) / 2 * z(5) - z(2) * z(3) / 2,
         79 * mp.pi ** 6 / 120960 - z(3) ** 2 / 4, mp.mpf(9) / 8 * z(7) - mp.mpf(3) / 8 * z(2) * z(5) - mp.mpf(9) / 32 * z(3) * z(4)]
for n in range(8):
    P_('[tau^%d]: taylor - claimed = %s' % (n, mp.nstr(tc[n] - claim[n], 3)))
P_('79 pi^6/120960 - 237/384 zeta(6) =', mp.nstr(79 * mp.pi ** 6 / 120960 - mp.mpf(237) / 384 * z(6), 3))
P_('claimed tau^7 coefficient in sympy variables: -pi^2 z5/16 - pi^4 z3/320 + 9 z7/8 minus claim =',
   mp.nstr(-mp.pi ** 2 * z(5) / 16 - mp.pi ** 4 * z(3) / 320 + 9 * z(7) / 8 - claim[7], 3))

# ---------- convexity ----------
u, t = sp.symbols('u t', positive=True)
Lam = sp.log(2 + sp.sin(u) - sp.cos(u))
hs = (sp.cos(u) + sp.sin(u)) * (1 + sp.exp(-2 * Lam / t)) ** t          # on 0 < u < pi/4
hh = hs + sp.diff(hs, u, 2)
hp = sp.diff(hs, u)
fhh = sp.lambdify((u, t), hh, 'mpmath')
fhp = sp.lambdify((u, t), hp, 'mpmath')


def author_formula(uu, tau):
    M = mp.cos(uu) + mp.sin(uu)
    Mp = mp.cos(uu) - mp.sin(uu)
    L = mp.log(2 + mp.sin(uu) - mp.cos(uu))
    y = 2 * L / tau
    sg = 1 / (1 + mp.exp(y))
    f = (1 + mp.exp(-y)) ** tau
    K = (1 - sg) / tau + sg
    return 2 * sg * f * M * mp.exp(-L) * (-3 * Mp + M ** 2 * mp.exp(-L) * (1 + 2 * K))


mp.mp.dps = 210
for tau in ['1', '0.5', '0.2', '0.1', '0.05', '0.02']:
    tau = mp.mpf(tau)
    Ngrid = 800
    grid = [mp.pi / 4 * k / Ngrid for k in range(Ngrid + 1)]
    grid[0] = mp.mpf(10) ** -60   # right limit at 0
    vals = [fhh(g, tau) for g in grid]
    kmin = min(range(len(vals)), key=lambda k: vals[k])
    # also minimum of the ratio (h+h'')/(sigma*...) is not needed; report min and where
    maxdev = max(abs(fhh(g, tau) - author_formula(g, tau)) / abs(fhh(g, tau)) for g in grid[1::40])
    # bracket lower bound check: bracket - 3e^{-L}*2(1-M') >= 0 on grid (K>=1 step)
    hp0 = fhp(mp.mpf(10) ** -100, tau)
    hpq = fhp(mp.pi / 4, tau)
    jump_num = -2 * hpq
    jump_formula = 4 * (1 + mp.mpf(4) ** (-1 / tau)) ** tau / (1 + mp.mpf(4) ** (1 / tau))
    P_('tau=%s: min_grid(h+h") = %s at u=%s ; all >0: %s ; rel.dev author formula %s' % (
        mp.nstr(tau, 3), mp.nstr(vals[kmin], 6), mp.nstr(grid[kmin], 6), all(vv > 0 for vv in vals), mp.nstr(maxdev, 3)))
    P_('     (h+h")(0+) = %s, predicted 2^tau*(1/tau-1) = %s' % (
        mp.nstr(fhh(mp.mpf(10) ** -100, tau), 8), mp.nstr(2 ** tau * (1 / tau - 1), 8)))
    P_("     h'(0+) = %s ; h'(pi/4-) = %s ; jump -2h'(pi/4-) = %s ; formula 4(1+4^(-1/t))^t/(1+4^(1/t)) = %s" % (
        mp.nstr(hp0, 3), mp.nstr(hpq, 8), mp.nstr(jump_num, 8), mp.nstr(jump_formula, 8)))

# (h+h'')(0+) value: from the factored form 2 sigma f M e^{-L}(1/tau-1) at u=0: sigma=1/2, f=2^tau
mp.mp.dps = 30
for tau in ['1.2', '1.5']:
    tau = mp.mpf(tau)
    P_('outside the claimed range, tau=%s: (h+h")(0+) = %s' % (tau, mp.nstr(fhh(mp.mpf(10) ** -20, tau), 6)))
open(__file__.replace('.py', '_output.txt'), 'w').write('\n'.join(out) + '\n')
