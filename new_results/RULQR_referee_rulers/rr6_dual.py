"""RULQR script 6: dual ruler. Taylor coefficients D_k of Gamma(1+h)^2/Gamma(1+2h) (= L_{p'}(p)/8) vs the claimed
1 - z(2)h^2 + 2 z(3)h^3 - (pi^4/40)h^4, and D_k(-1)^k -> residue behaviour of the simple pole at h = -1 (radius 1);
plus the support-function identity int h_K(n) ds = 2 Area(K) evaluated numerically for the l^p ball at p = 3, 7."""
import mpmath as mp
mp.mp.dps = 40
D = mp.taylor(lambda h: mp.gamma(1 + h) ** 2 / mp.gamma(1 + 2 * h), 0, 40)
print('D_2 + zeta(2) =', mp.nstr(D[2] + mp.zeta(2), 5), '  D_3 - 2 zeta(3) =', mp.nstr(D[3] - 2 * mp.zeta(3), 5),
      '  D_4 + pi^4/40 =', mp.nstr(D[4] + mp.pi ** 4 / 40, 5))
res = mp.limit(lambda h: (1 + h) * mp.gamma(1 + h) ** 2 * mp.rgamma(1 + 2 * h), -1)
print('residue at h=-1:', mp.nstr(res, 12), ' (-1)^k D_k for k = 10, 20, 40:', [mp.nstr((-1) ** k * D[k], 10) for k in (10, 20, 40)])
for p in [3, 7]:
    p = mp.mpf(p); pp = p / (p - 1)
    # boundary point (cos^{2/p} t, sin^{2/p} t); outward normal direction grad(|x|^p+|y|^p) = (x^{p-1}, y^{p-1})
    def integrand(t):
        x, y = mp.cos(t) ** (2 / p), mp.sin(t) ** (2 / p)
        xd, yd = -2 / p * mp.cos(t) ** (2 / p - 1) * mp.sin(t), 2 / p * mp.sin(t) ** (2 / p - 1) * mp.cos(t)
        ds = mp.hypot(xd, yd); nx, ny = x ** (p - 1), y ** (p - 1); nn = mp.hypot(nx, ny)
        return ((nx / nn) ** pp + (ny / nn) ** pp) ** (1 / pp) * ds     # h_K(n) = ||n||_{p'}
    I = 4 * mp.quad(integrand, [0, mp.pi / 4, mp.pi / 2])
    print('p=%s  int h_K(n) ds = %s   2 Area = 8 G(1+1/p)^2/G(1+2/p) = %s' % (p, mp.nstr(I, 25), mp.nstr(8 * mp.gamma(1 + 1 / p) ** 2 / mp.gamma(1 + 2 / p), 25)))
