"""Referee PUR2R numerics: perimeter by u-quadrature vs closed form minus tail; h+h'' minima;
polygon reconstruction of the curve from h."""
import sys, time
import mpmath as mp
from r_curve import Body

def out(*a):
    print(*a); sys.stdout.flush()

mp.mp.dps = 60
PI4 = mp.pi / 4
bodies = [('c=1', Body('c1')), ('kappa^2=9/7', Body('imag', mp.mpf(9) / 7))]
for name, B in bodies:
    out("%s: V = %s, Phi(V)-1 = %s, phi(V) = %s" % (name, mp.nstr(B.V, 30), mp.nstr(B.Phi(B.V) - 1, 3), mp.nstr(B.phi(B.V), 20)))
out("c=1 Lambert check: V - (5 - W(e^5))/2 =", mp.nstr(bodies[0][1].V - (5 - mp.lambertw(mp.e**5).real) / 2, 3))
uu = mp.mpf('0.37')
qq = 5 + 4 * mp.sin(uu) - 4 * mp.cos(uu)
out("c=1 Lambert check at u=0.37: Lam - (q - W(e^q))/2 =", mp.nstr(bodies[0][1].Lam(uu) - (qq - mp.lambertw(mp.e**qq).real) / 2, 3))

# Lemma 1 at complex parameters (own check)
al, be = mp.mpc('0.3', '0.7'), mp.mpc('-0.45', '-0.7')
ta = al + be
lhs = mp.quad(lambda s: (mp.exp(2 * al * s) + mp.exp(2 * be * s)) * ((1 + mp.exp(-2 * s)) ** ta - 1), [0, 1, 5, 20, 60, mp.inf])
rhs = (mp.gamma(-al) * mp.gamma(-be) / mp.gamma(-ta) + 1 / al + 1 / be) / 2
out("Lemma 1 at alpha=0.3+0.7i, beta=-0.45-0.7i (tau=-0.15 complex-free test): lhs-rhs =", mp.nstr(abs(lhs - rhs), 3))
al, be = mp.mpc('0.6', '0.9'), mp.mpc('0.2', '-0.4')
ta = al + be
lhs = mp.quad(lambda s: (mp.exp(2 * al * s) + mp.exp(2 * be * s)) * ((1 + mp.exp(-2 * s)) ** ta - 1), [0, 1, 5, 20, 60, 150, mp.inf])
rhs = (mp.gamma(-al) * mp.gamma(-be) / mp.gamma(-ta) + 1 / al + 1 / be) / 2
out("Lemma 1 at alpha=0.6+0.9i, beta=0.2-0.4i: lhs-rhs =", mp.nstr(abs(lhs - rhs), 3))

out("\n== perimeter: 8*int_0^{pi/4} h du and full-circle int_0^{2pi} h du vs closed form - tail ==")
for name, B in bodies:
    for tau in [mp.mpf('0.3'), mp.mpf('0.9')]:
        t0 = time.time()
        pts = [0, mp.mpf('0.02'), mp.mpf('0.08'), mp.mpf('0.25'), mp.mpf('0.5'), PI4]
        I8 = 8 * mp.quad(lambda u: B.h_plus(u, tau), pts)
        cf, E, main = B.closed_minus_tail(tau)
        out("%s tau=%s: 8*int h = %s" % (name, tau, mp.nstr(I8, 50)))
        out("    closed - E      = %s  (imag part %s)" % (mp.nstr(mp.re(cf), 50), mp.nstr(mp.im(cf), 3)))
        out("    difference = %s ; E = %s ; (%.0fs)" % (mp.nstr(mp.re(cf) - I8, 3), mp.nstr(mp.re(E), 10), time.time() - t0))
    # full circle at tau=0.3, dps 40, no symmetry reduction
    mp.mp.dps = 40
    tau = mp.mpf('0.3')
    bps = [k * mp.pi / 4 for k in range(9)]
    pts = []
    for k in range(8):
        a, b = bps[k], bps[k + 1]
        pts += [a, a + (b - a) * mp.mpf('0.03'), a + (b - a) * mp.mpf('0.3'), b - (b - a) * mp.mpf('0.3'), b - (b - a) * mp.mpf('0.03')]
    pts.append(bps[8])
    Ifull = mp.quad(lambda u: B.h(u, tau), pts)
    cf, E, main = B.closed_minus_tail(tau)
    out("%s tau=0.3 full circle (dps 40): int_0^2pi h - (closed - E) = %s" % (name, mp.nstr(Ifull - mp.re(cf), 3)))
    mp.mp.dps = 60
    B.cache.clear()

# c=1 at tau = 1 and near 1
B = bodies[0][1]
I1 = 8 * mp.quad(lambda u: B.h_plus(u, 1), [0, mp.mpf('0.25'), mp.mpf('0.5'), PI4])
out("\nc=1 tau=1: 8*int h - (8+4V+2(1-e^{-2V})) =", mp.nstr(I1 - (8 + 4 * B.V + 2 * (1 - mp.exp(-2 * B.V))), 3), " value", mp.nstr(I1, 30))
for eps in ['1e-10', '1e-20']:
    cf, E, main = B.closed_minus_tail(1 - mp.mpf(eps))
    out("   closed-E at tau = 1-%s minus tau=1 value: %s" % (eps, mp.nstr(cf - I1, 3)))
