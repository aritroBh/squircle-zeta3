"""Quadrature + Richardson checks.
 (c) exactness of P = 8 int_0^inf phi (1+e^{-2v/tau})^tau dv against theta-quadrature of h_p,
     then Richardson of F(tau) for c6, c7 of the squircle.
 (d) soft-max F with V = 1/sqrt2: Richardson c6, c7.
 (e) test example psi = (1+v) e^{-v^2} (omega = 1, j = 0): Richardson F4.
 (g) p-norm rhombus diag(1,2): exactness + Richardson c2..c5.
"""
import sys, pickle
import mpmath as mp

mp.mp.dps = 120


def Bfac(s, tau):
    return mp.expm1(tau * mp.log1p(mp.exp(-2 * s)))


def F_int(phi, tau, V=mp.inf):
    """F = tau int_0^{V/tau} phi(tau s) B(s,tau) ds."""
    smax = mp.mpf(150) if V == mp.inf else min(mp.mpf(150), V / tau)
    pts = [mp.mpf(0)] + [mp.mpf(2)**k / 8 for k in range(0, 11) if mp.mpf(2)**k / 8 < smax] + [smax]
    return tau * mp.quad(lambda s: phi(tau * s) * Bfac(s, tau), pts)


def neville(taus, Ys):
    """Richardson/Neville extrapolation to tau=0 of Y(tau) = c0 + c1 tau + ...; returns table of last-entry values."""
    n = len(taus)
    P = list(Ys)
    best = [P[-1]]
    for k in range(1, n):
        P = [(taus[i] * P[i + 1] - taus[i + k] * P[i]) / (taus[i] - taus[i + k]) for i in range(n - k)]
        best.append(P[-1])
    return best


def richardson(phi, V, known, nstart, taus):
    """Y = (F - sum_{n<nstart} known[n] tau^n)/tau^nstart -> extrapolate."""
    Fs = [F_int(phi, t, V) for t in taus]
    Ys = [(F - sum(known[n] * t**n for n in range(nstart))) / t**nstart for F, t in zip(Fs, taus)]
    return Fs, Ys, neville(taus, Ys)


if __name__ == '__main__':
    which = sys.argv[1]
    out = {}
    if which == 'squircle':
        # exactness check
        mp.mp.dps = 50
        for tau in ['0.5', '0.2', '0.1', '0.05']:
            tau = mp.mpf(tau)
            q = 1 / tau
            h = lambda th: (abs(mp.cos(th))**q + abs(mp.sin(th))**q)**tau
            # P = 8 int_{pi/4}^{pi/2} h ; cluster at pi/4
            pts = [mp.pi / 4 + mp.pi / 4 * mp.mpf(2)**(-k) for k in range(30, -1, -1)]
            Pth = 8 * mp.quad(h, [mp.pi / 4] + pts)
            phi = lambda w: mp.exp(w) / (mp.sqrt(2) * mp.cosh(2 * w)**mp.mpf(1.5))
            Pv = 4 * mp.sqrt(2) + 8 * F_int(phi, tau)
            print("tau=%s  P(theta-quad)=%s  P(v-form)-P(theta) = %s" % (tau, mp.nstr(Pth, 45), mp.nstr(Pv - Pth, 3)))
        mp.mp.dps = 120
        phi = lambda w: mp.exp(w) / (mp.sqrt(2) * mp.cosh(2 * w)**mp.mpf(1.5))
        r2 = mp.sqrt(2)
        z = mp.zeta
        c = {0: 0, 1: 0, 2: r2 * z(2), 3: r2 * z(3), 4: -mp.mpf(33) / 16 * r2 * z(4),
             5: -r2 * (mp.mpf(3) / 4 * z(5) + mp.mpf(7) / 4 * z(2) * z(3))}
        c6 = r2 * (mp.mpf(4039) / 256 * z(6) - mp.mpf(7) / 8 * z(3)**2)
        c7 = r2 * (mp.mpf(129) / 16 * z(7) + mp.mpf(141) / 16 * z(2) * z(5) + mp.mpf(1071) / 64 * z(3) * z(4))
        known = {n: c[n] / 8 for n in c}  # F_n
        taus = [mp.mpf(1) / (20 * 2**i) for i in range(0, 15)]
        Fs, Ys, best = richardson(phi, mp.inf, known, 6, taus)
        print("\nsquircle: Richardson of (F - sum_{n<=5} F_n tau^n)/tau^6 -> F_6 ; c6 = 8 F6")
        for k, b in enumerate(best):
            print(" order %2d  8*Y = %s   minus c6(exact) = %s" % (k, mp.nstr(8 * b, 45), mp.nstr(8 * b - c6, 3)))
        known[6] = c6 / 8
        Ys7 = [(F - sum(known[n] * t**n for n in range(7))) / t**7 for F, t in zip(Fs, taus)]
        best7 = neville(taus, Ys7)
        print("\nsquircle: Richardson -> F_7 ; c7 = 8 F7")
        for k, b in enumerate(best7):
            print(" order %2d  8*Y = %s   minus c7(exact) = %s" % (k, mp.nstr(8 * b, 45), mp.nstr(8 * b - c7, 3)))
    elif which == 'softmax':
        phi = lambda w: 1 + w / mp.sqrt(1 - w**2)
        V = 1 / mp.sqrt(2)
        z = mp.zeta
        c = {0: 0, 1: 0, 2: mp.pi**2 / 3, 3: 2 * z(3), 4: z(4) / 4, 5: mp.mpf(107) / 32 * z(5) - z(2) * z(3)}
        known = {n: c[n] / 8 for n in c}
        taus = [mp.mpf(1) / (200 * 2**i) for i in range(0, 14)]
        Fs, Ys, best = richardson(phi, V, known, 6, taus)
        pickle.dump([(str(t), str(F)) for t, F in zip(taus, Fs)], open('softmax_F.pkl', 'wb'))
        print("soft-max: Richardson -> c6 = 8 F6")
        for k, b in enumerate(best):
            print(" order %2d  c6 ~ %s" % (k, mp.nstr(8 * b, 45)))
    elif which == 'test_e':
        # psi = (1+v) e^{-v^2}: phi = e^v (1+v) e^{-v^2}
        phi = lambda w: mp.exp(w) * (1 + w) * mp.exp(-w**2)
        known = pickle.load(open('test_e_known.pkl', 'rb'))
        known = {int(k): mp.mpf(vv) for k, vv in known.items()}
        taus = [mp.mpf(1) / (20 * 2**i) for i in range(0, 15)]
        Fs, Ys, best = richardson(phi, mp.inf, known, 4, taus)
        print("test (e), psi=(1+v)e^{-v^2}: Richardson -> F_4")
        for k, b in enumerate(best):
            print(" order %2d  F4 ~ %s" % (k, mp.nstr(b, 45)))
    elif which == 'rhombus':
        a, b = mp.mpf(1), mp.mpf(2)
        qf = lambda w: 2 / (a * b * (mp.exp(2 * w) / a**2 + mp.exp(-2 * w) / b**2)**mp.mpf(1.5))
        phi = lambda w: mp.exp(w) * (qf(w) + qf(-w))
        mp.mp.dps = 50
        for tau in ['0.5', '0.2', '0.1']:
            tau = mp.mpf(tau)
            q = 1 / tau
            h = lambda th: (abs(a * mp.cos(th))**q + abs(b * mp.sin(th))**q)**tau
            th0 = mp.atan(a / b)  # kink where a|cos| = b|sin|
            pts = sorted(set([mp.mpf(0)] + [th0 + s * th0 * mp.mpf(2)**(-k) for k in range(0, 30) for s in (-1, 1)
                                             if 0 <= th0 + s * th0 * mp.mpf(2)**(-k) <= mp.pi / 2] + [th0, mp.pi / 2]))
            Pth = 4 * mp.quad(h, pts)
            r = mp.sqrt(a * a + b * b)
            Pv = 4 * r + 4 * F_int(phi, tau)
            print("rhombus a=1,b=2 tau=%s P(theta)=%s  P(v-form)-P(theta)=%s" % (tau, mp.nstr(Pth, 40), mp.nstr(Pv - Pth, 3)))
        mp.mp.dps = 120
        known = pickle.load(open('rhombus_known.pkl', 'rb'))
        known = {int(k): mp.mpf(vv) for k, vv in known.items()}
        taus = [mp.mpf(1) / (20 * 2**i) for i in range(0, 15)]
        for nst in (4, 5):
            Fs, Ys, best = richardson(phi, mp.inf, known, nst, taus)
            print("rhombus: Richardson -> P coefficient c_%d = 4 F_%d ; predicted %s" % (nst, nst, mp.nstr(4 * known[nst], 40)))
            for k, bb in enumerate(best):
                if k % 2 == 0 or k == len(best) - 1:
                    print(" order %2d  c ~ %s   diff %s" % (k, mp.nstr(4 * bb, 42), mp.nstr(4 * bb - 4 * known[nst], 3)))
