# Referee check of (c): direct theta-quadrature of per(B_p) vs the v-form, and extraction of c6, c7
# from theta-quadrature alone (Neville extrapolation in tau).
import mpmath as mp

mp.mp.dps = 90
out = []
P_ = lambda *a: (out.append(' '.join(str(t) for t in a)), print(*a))
z = mp.zeta
s2 = mp.sqrt(2)
c = {0: 4 * s2, 1: 0, 2: s2 * z(2), 3: s2 * z(3), 4: -mp.mpf(33) / 16 * s2 * z(4),
     5: -s2 * (mp.mpf(3) / 4 * z(5) + mp.mpf(7) / 4 * z(2) * z(3)),
     6: s2 * (mp.mpf(4039) / 256 * z(6) - mp.mpf(7) / 8 * z(3) ** 2),
     7: s2 * (mp.mpf(129) / 16 * z(7) + mp.mpf(141) / 16 * z(2) * z(5) + mp.mpf(1071) / 64 * z(3) * z(4))}


def per_theta(tau):
    # P = 8 * int_0^{pi/4} cos(t) (1 + tan(t)^{1/tau})^tau dt   (h_p = (|cos|^q+|sin|^q)^{1/q}, q=1/tau)
    q = 1 / tau
    f = lambda t: mp.cos(t) * (1 + mp.tan(t) ** q) ** tau
    a = mp.pi / 4
    pts = [0] + sorted(set([a - tau * k for k in (64, 32, 16, 8, 4, 2, 1, mp.mpf(1) / 2) if a - tau * k > 0])) + [a]
    return 8 * mp.quad(f, pts)


def per_v(tau):
    phi = lambda v: mp.exp(v) / (s2 * mp.cosh(2 * v) ** mp.mpf(1.5))
    f = lambda v: phi(v) * (1 + mp.exp(-2 * v / tau)) ** tau
    return 8 * mp.quad(f, [0, tau, 4 * tau, 16 * tau, 1, 4, 20, mp.inf])


mp.mp.dps = 50
for tau in [mp.mpf('0.5'), mp.mpf('0.3'), mp.mpf('0.1'), mp.mpf('0.03')]:
    a, b = per_theta(tau), per_v(tau)
    P_('tau=%s  P_theta=%s  P_v - P_theta=%s' % (tau, mp.nstr(a, 30), mp.nstr(b - a, 3)))
P_('tau=1/2: P_theta - 2pi =', mp.nstr(per_theta(mp.mpf(1) / 2) - 2 * mp.pi, 3))

mp.mp.dps = 90
K = 14
taus = [mp.mpf(1) / (20 * 2 ** i) for i in range(K)]
vals = []
for t in taus:
    Pt = per_theta(t)
    vals.append(Pt)


def neville0(xs, ys):
    ys = list(ys)
    n = len(xs)
    for m in range(1, n):
        for i in range(n - m):
            ys[i] = (xs[i] * ys[i + 1] - xs[i + m] * ys[i]) / (xs[i] - xs[i + m])
    return ys[0]


for target in (6, 7):
    g = [(Pt - mp.fsum(c[n] * t ** n for n in range(target))) / t ** target for t, Pt in zip(taus, vals)]
    for npts in (8, 10, 12, 14):
        est = neville0(taus[K - npts:], g[K - npts:])
        P_('c_%d from theta-quadrature, Neville on %d smallest taus: %s ; minus exact = %s' % (
            target, npts, mp.nstr(est, 40), mp.nstr(est - c[target], 3)))
open(__file__.replace('.py', '_output.txt'), 'w').write('\n'.join(out) + '\n')
