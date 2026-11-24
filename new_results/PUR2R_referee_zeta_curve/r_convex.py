"""Referee PUR2R: h+h'' minima (formula and finite differences), h' jump at pi/4,
curve reconstruction from h and polygon length, spot checks of E1/E2/E3."""
import sys, time
import mpmath as mp
from r_curve import Body

def out(*a):
    print(*a); sys.stdout.flush()

PI4 = mp.pi / 4
mp.mp.dps = 30
bodies = [('c=1', Body('c1')), ('kappa^2=9/7', Body('imag', mp.mpf(9) / 7))]

out("== h+h'': my formula vs mp.diff (dps 40) at 12 points ==")
mp.mp.dps = 40
for name, B in bodies:
    B.cache.clear()
    worst = 0
    for tau in [mp.mpf('0.3'), mp.mpf('0.9')]:
        for i in range(1, 13):
            u = PI4 * mp.mpf(i) / 13
            num = B.h_plus(u, tau) + mp.diff(lambda x: B.h_plus(x, tau), u, 2)
            fo = B.hhpp_plus(u, tau)
            worst = max(worst, abs(num - fo) / abs(fo))
            hp_num = mp.diff(lambda x: B.h_plus(x, tau), u, 1)
            worst = max(worst, abs(hp_num - B.hp_plus(u, tau)) / (abs(hp_num) + 1e-30))
    out("  %s: max rel. diff (h+h'' and h') = %s" % (name, mp.nstr(worst, 3)))
    B.cache.clear()

mp.mp.dps = 25
out("\n== min of h+h'' on a grid of 4000 points u = (pi/4) i/4000, i=1..4000 (formula) ==")
for name, B in bodies:
    for tau in ['0.3', '0.9', '1', '0.1']:
        tau = mp.mpf(tau)
        vals = [(B.hhpp_plus(PI4 * i / 4000, tau), i) for i in range(1, 4001)]
        m, i = min(vals)
        out("  %s tau=%s: min h+h'' = %s at u = %s ; value at u->0 limit 2^tau(1/tau-1) = %s" % (
            name, mp.nstr(tau, 3), mp.nstr(m, 8), mp.nstr(PI4 * i / 4000, 6), mp.nstr(2**tau * (1 / tau - 1), 8)))
    # jump of h' at pi/4
    for tau in [mp.mpf('0.3'), mp.mpf('0.9')]:
        mp.mp.dps = 40
        B.cache.clear()
        left = mp.diff(lambda x: B.h_plus(x, tau), PI4, 1, direction=-1)
        V = B.V
        f = (1 + mp.exp(-2 * V / tau)) ** tau
        sig = 1 / (1 + mp.exp(2 * V / tau))
        out("  %s tau=%s: jump -2h'(pi/4-) = %s ; claimed 8 f sigma/phi(V) = %s" % (name, tau, mp.nstr(-2 * left, 20), mp.nstr(8 * f * sig / B.phi(V), 20)))
        mp.mp.dps = 25
        B.cache.clear()
    # tau > 1
    out("  %s tau=1.2: (h+h'')(1e-6) = %s vs 2^1.2(1/1.2-1) = %s" % (name, mp.nstr(B.hhpp_plus(mp.mpf('1e-6'), mp.mpf('1.2')), 8), mp.nstr(mp.mpf(2)**1.2 * (1 / mp.mpf(1.2) - 1), 8)))

out("\n== E1 check c=1: min over grid of (D phi - 6p) ; E2 imaginary: D phi vs 6p - 2 kappa tan ==")
for name, B in bodies:
    worst = mp.inf
    for i in range(1, 2001):
        u = PI4 * i / 2000
        L = B.Lam(u); M = mp.cos(u) + mp.sin(u); Mp = mp.cos(u) - mp.sin(u); ph = B.phi(L)
        D = -3 * Mp + (M * M / ph) * (2 + B.dphi(L) / ph)
        p = 1 + mp.sin(u) - mp.cos(u)
        if B.kind == 'c1':
            worst = min(worst, D * ph - 6 * p)
        else:
            worst = min(worst, D * ph - (6 * p - 2 * B.k * mp.tan(B.k * L)))
    out("  %s: min (D phi - lower bound) = %s (should be >= 0)" % (name, mp.nstr(worst, 6)))

out("\n== curve reconstruction x = h cos u - h' sin u, y = h sin u + h' cos u ==")
def polygon(B, tau, N):
    P = []
    for k in range(8):
        a = k * mp.pi / 4
        for i in range(N + 1):
            u = a + (mp.pi / 4) * i / N
            side = +1 if i == 0 else (-1 if i == N else 0)
            hh = B.h(u, tau); hp = B.hp(u, tau, side)
            P.append((hh * mp.cos(u) - hp * mp.sin(u), hh * mp.sin(u) + hp * mp.cos(u)))
    # drop exact duplicates at octant junctions where h' is continuous
    Q = [P[0]]
    for q in P[1:]:
        if abs(q[0] - Q[-1][0]) + abs(q[1] - Q[-1][1]) > mp.mpf(10) ** (-20):
            Q.append(q)
    if abs(Q[-1][0] - Q[0][0]) + abs(Q[-1][1] - Q[0][1]) < mp.mpf(10) ** (-20):
        Q.pop()
    n = len(Q)
    L = 0; turn = 0; mincross = mp.inf; edges = []
    for i in range(n):
        x0, y0 = Q[i]; x1, y1 = Q[(i + 1) % n]
        edges.append((x1 - x0, y1 - y0))
        L += mp.hypot(x1 - x0, y1 - y0)
    for i in range(n):
        e0, e1 = edges[i], edges[(i + 1) % n]
        cr = e0[0] * e1[1] - e0[1] * e1[0]
        dt = e0[0] * e1[0] + e0[1] * e1[1]
        mincross = min(mincross, cr / (mp.hypot(*e0) * mp.hypot(*e1)))
        turn += mp.atan2(cr, dt)
    return L, turn / (2 * mp.pi), mincross, n

for name, B in bodies:
    for tau in [mp.mpf('0.3'), mp.mpf('0.9')]:
        cf = mp.re(B.closed_minus_tail(tau)[0])
        Ls = []
        for N in [100, 200, 400, 800]:
            t0 = time.time()
            L, wind, mc, n = polygon(B, tau, N)
            Ls.append(L)
            out("  %s tau=%s N/octant=%d: vertices %d, length = %s, L - perimeter = %s, turning/2pi = %s, min normalized cross = %s (%.0fs)" % (
                name, tau, N, n, mp.nstr(L, 18), mp.nstr(L - cf, 4), mp.nstr(wind, 15), mp.nstr(mc, 4), time.time() - t0))
        R1 = [(4 * Ls[i + 1] - Ls[i]) / 3 for i in range(3)]
        R2 = [(16 * R1[i + 1] - R1[i]) / 15 for i in range(2)]
        out("    Richardson: %s ; minus perimeter = %s" % (mp.nstr(R2[-1], 18), mp.nstr(R2[-1] - cf, 3)))
        B.cache.clear()

out("\n== E2(iii) spot checks (numerical, grid 600) ==")
mp.mp.dps = 20
for k2 in ['1.69', '1.96']:
    B = Body('imag', mp.mpf(k2))
    grid = [PI4 * i / 600 for i in range(1, 601)]
    def Dm(u):
        L = B.Lam(u); M = mp.cos(u) + mp.sin(u); Mp = mp.cos(u) - mp.sin(u); ph = B.phi(L)
        return (-3 * Mp + (M * M / ph) * (2 + B.dphi(L) / ph)) * ph / (M * M)
    mD = min(Dm(u) for u in grid)
    ok = lambda tau: min(B.hhpp_plus(u, tau) for u in grid) >= 0
    if ok(mp.mpf(1)):
        ts = 1
    else:
        a, b = mp.mpf('0.01'), mp.mpf(1)
        for _ in range(25):
            m = (a + b) / 2
            if ok(m): a = m
            else: b = m
        ts = a
    t0 = 1 / (1 + B.k * mp.tan(B.k * B.V))
    out("  kappa=%s: phi(V)=%s minDphi/M^2=%s tau0=%s tau*=%s" % (mp.nstr(B.k, 4), mp.nstr(B.phi(B.V), 4), mp.nstr(mD, 4), mp.nstr(t0, 4), mp.nstr(ts, 4)))
