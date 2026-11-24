"""RULQ script 4 (task B): evaluation of the rigorous truncation bound of Theorem B2 (result.md Sec. 5).

 |R_K(p)| <= R0 + R1 + R2,   h = 1/p,
 R0 = 2 (h log 2)^{K+1}/(K+1)!
 R1 = h^2 [ sum_{j=0}^{K-4} h^j rho^{-j} Bnd_{K-2-j}(m1; M1, a1, r1) + 2 M1 h (h/rho)^{K-3}/((1-h/rho)(1-a1 h)) ]
 R2 =       sum_{j=2}^{K-2} h^j rho2^{-j} Bnd_{K-j}(m2; M2, a2, r2) + 2 M2 h (h/rho2)^{K-1}/((1-h/rho2)(1-a2 h))
 Bnd_N(m; M, a, r) = (1/(pi r)) int_R m(x) (x^2+r^2)^{-N/2} h^{N+1} Gamma(N+1, (x-r)^+/h) dx
                   + (M/pi) e^{a r} int_0^oo e^{-y(1/h-a)} y^N (y+r)^{-N} dy          (Lemma R)
 n = 1:  r1 = pi - delta, m1(x) = e^{2x} e^{rho |log(1+e^{x+i r1})|} / |1+e^{x+i r1}|  (exact modulus and argument),
         M1 = e^{rho Lam}/sigma1, a1 = 1+rho, sigma1 = sqrt2 sin(delta/2), Lam = pi + log(2/sigma1).
 n >= 2: r2 = 2pi - delta2, m2(x) = M2 e^{a2 x^+}, M2 = rho2^2 e^{rho2(1+pi)} sigma2^{-1-rho2} (zeta(2-rho2)-1),
         a2 = 1 + rho2/2, sigma2 = sqrt2 sin(delta2/4).
Every choice of rho in (h, 1/(1+... )) with (1+rho)h<1, delta in (0, pi/2], rho2 in (h,1), delta2 in (0,2pi) gives a valid bound;
we take the best of a small grid.  Usage: ru4_bound.py"""
import mpmath as mp
mp.mp.dps = 25
N_ = [mp.mpf(l) for l in open('../PIPDIV_divergence_proof/pd1_N_K100_dps110_deg6_S760.txt')]

def upgam(Nn, z, h):   # h^{N+1} Gamma(N+1, z) for integer N, z >= 0 (finite sum)
    return h ** (Nn + 1) * mp.factorial(Nn) * mp.exp(-z) * mp.fsum(z ** i / mp.factorial(i) for i in range(Nn + 1))
def Bnd(Nn, m, M, a, r, h):
    # horizontal edges on [-XL, XR] by quadrature, plus explicit tail bounds (result.md Sec. 5, Lemma R)
    XL = 50 * r; XR = r + 4 * Nn * h + 40 * h + 10
    f = lambda x: m(x) * (x * x + r * r) ** (-mp.mpf(Nn) / 2) * upgam(Nn, max(x - r, 0) / h, h)
    pts = [-XL, -4 * r, -r, -r / 4, 0, r / 4, r, 2 * r, r + Nn * h, XR]
    H = mp.quad(f, sorted(set(pts)))
    tailL = M * mp.factorial(Nn) * h ** (Nn + 1) * XL ** (1 - Nn) / (Nn - 1)
    th = mp.mpf(9) / 10; b = th / h - a      # Gamma(N+1,z) <= e^{-th z} N! (1-th)^{-N-1}
    assert b > 0
    tailR = M * mp.factorial(Nn) * h ** (Nn + 1) * (1 - th) ** (-Nn - 1) * XR ** (-Nn) * mp.exp(th * r / h) * mp.exp(-b * XR) / b
    H = (H + tailL + tailR) / (mp.pi * r)
    g = lambda y: mp.exp(-y * (1 / h - a)) * y ** Nn * (y + r) ** (-Nn)
    V = M / mp.pi * mp.exp(a * r) * mp.quad(g, [0, r / 4, r, 4 * r, 40 * r + 40 * h * Nn, mp.inf])
    return H + V

def bound(p, K, rho, delta, rho2, delta2):
    p = mp.mpf(p); h = 1 / p
    if not ((1 + rho) * h < mp.mpf(9) / 10 and h < rho and (1 + rho) * h < 1 and h < rho2 < 1 and 0 < delta <= mp.pi / 2): return None
    R0 = 2 * (h * mp.log(2)) ** (K + 1) / mp.factorial(K + 1)
    r1 = mp.pi - delta; s1 = mp.sqrt(2) * mp.sin(delta / 2); Lam = mp.pi + mp.log(2 / s1)
    M1 = mp.exp(rho * Lam) / s1; a1 = 1 + rho
    def m1(x):
        z = 1 + mp.exp(mp.mpc(x, r1))
        return mp.exp(2 * x) * mp.exp(rho * abs(mp.log(z))) / abs(z)
    R1 = mp.fsum(h ** j * rho ** (-j) * Bnd(K - 2 - j, m1, M1, a1, r1, h) for j in range(0, K - 3))
    R1 = h ** 2 * (R1 + 2 * M1 * h * (h / rho) ** (K - 3) / ((1 - h / rho) * (1 - a1 * h)))
    r2 = 2 * mp.pi - delta2; s2 = mp.sqrt(2) * mp.sin(delta2 / 4)
    M2 = rho2 ** 2 * mp.exp(rho2 * (1 + mp.pi)) * s2 ** (-1 - rho2) * (mp.zeta(2 - rho2) - 1); a2 = 1 + rho2 / 2
    m2 = lambda x: M2 * mp.exp(a2 * max(x, 0))
    R2 = mp.fsum(h ** j * rho2 ** (-j) * Bnd(K - j, m2, M2, a2, r2, h) for j in range(2, K - 1))
    R2 += 2 * M2 * h * (h / rho2) ** (K - 1) / ((1 - h / rho2) * (1 - a2 * h))
    return R0 + R1 + R2, R0, R1, R2

if __name__ == '__main__':
    import sys
    PLIST = [int(x) for x in sys.argv[1:]] or [5, 10, 20]
    vals = {2: mp.mpf('0.78539816339744830961566084581987572104929235'),
            3: mp.mpf('0.814941998264748774283876876028429574295518671'),
            5: mp.mpf('0.874614202011362268979564321159152575890780092'),
            10: mp.mpf('0.933490202305707074895068494959379907511514221'),
            20: None}
    print(' p   K   actual R_K      bound          bound/|R_K|   bound*p^(5/2)e^(pi p)   (R0, R1, R2)   params')
    for p in PLIST:
        h = 1 / mp.mpf(p)
        t = [N_[k] * h ** (k + 1) for k in range(len(N_))]
        Ks = min(range(2, len(N_) - 1), key=lambda k: abs(t[k]))
        RK = (vals[p] - 1 - mp.fsum(t[:Ks])) if vals[p] is not None else None
        best = None
        rhos = [h * c for c in (2.5, 3, 3.5, 4, 5, 6, 8) if h * c < 1]
        for rho in rhos:
            for dc in (0.5, 1, 1.5, 2.5):
                delta = min(mp.mpf(dc) / p, mp.pi / 2)
                b = bound(p, Ks, rho, delta, mp.mpf(0.9), mp.mpf(1.5))
                if b is not None and (best is None or b[0] < best[0][0]): best = (b, rho, delta)
        b, rho, delta = best
        for rho2 in (0.6, 0.9):
            for d2 in (1, 2):
                bb = bound(p, Ks, rho, delta, mp.mpf(rho2), mp.mpf(d2))
                if bb is not None and bb[0] < best[0][0]: best = (bb, rho, delta, rho2, d2)
        b = best[0]
        print('%3d %3d  %s  %.4e  %s  %.3e   (%.1e, %.2e, %.2e)  rho=%s delta=%s' % (
            p, Ks, ('% .4e' % float(RK)) if RK is not None else '     n/a   ', float(b[0]),
            ('%.3e' % float(b[0] / abs(RK))) if RK is not None else '   n/a   ',
            float(b[0] * mp.mpf(p) ** 2.5 * mp.exp(mp.pi * p)), float(b[1]), float(b[2]), float(b[3]),
            mp.nstr(best[1], 3), mp.nstr(best[2], 3)), flush=True)
