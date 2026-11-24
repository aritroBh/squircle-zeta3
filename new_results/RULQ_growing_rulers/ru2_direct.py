"""RULQ script 2: direct (quadrature) checks.
 (a) L_q(p)/8 by two independent quadratures: v-form split  2^{-h} + h int_0^oo e^{-hv}(1+e^{-v})^{-1-h}[(1+e^{-rv})^{1/q}-1] dv
     (r = q(1-h)) and the w-form  int_0^{2^{-1/p}} (1 + (w^p/(1-w^p))^{q(1-1/p)})^{1/q} dw  (t = w^p in the t-integral).
 (b) growing rulers: compare with the optimally truncated series of ru1 (files ru1_N_*.txt): the difference must be
     of the size of the smallest term ~ e^{-pi lambda p}.
 (c) dual ruler q = p/(p-1): L = 8 Gamma(1+1/p)^2/Gamma(1+2/p) = (4/p) B(1/p,1/p), and its convergent series.
Usage: ru2_direct.py"""
import mpmath as mp
mp.mp.dps = 60

def L8_v(p, q):
    p = mp.mpf(p); q = mp.mpf(q); h = 1 / p; r = q * (1 - h)
    f = lambda v: mp.exp(-h * v) * (1 + mp.exp(-v)) ** (-1 - h) * mp.expm1(mp.log1p(mp.exp(-r * v)) / q)
    pts = [0] + [mp.mpf(4) ** i / r for i in range(0, 6)] + [mp.mpf(1), 4, 16, 64, 256, 1024]
    pts = sorted(set(pts))
    return mp.mpf(2) ** (-h) + h * mp.quad(f, pts)

def L8_w(p, q):
    p = mp.mpf(p); q = mp.mpf(q); e = q * (1 - 1 / p); W = mp.mpf(2) ** (-1 / p)
    f = lambda w: (1 + (w ** p / (1 - w ** p)) ** e) ** (1 / q)
    pts = [0, W / 2, W * (1 - mp.mpf(1) / 4), W * (1 - mp.mpf(1) / 16), W * (1 - mp.mpf(1) / 64), W * (1 - mp.mpf(1) / 512), W]
    return mp.quad(f, pts) / 1   # = L/8

print('(a) two quadratures of L_q(p)/8')
for (p, q) in [(3, 6), (3, 1.5), (5, 6), (8, 16), (8, 4), (6, 7), (4, mp.mpf(4) / 3)]:
    a = L8_v(p, q); b = L8_w(p, q)
    print('  p=%s q=%s  L/8 = %s   |v-form - w-form| = %s' % (p, mp.nstr(mp.mpf(q), 6), mp.nstr(a, 40), mp.nstr(abs(a - b), 3)))

print('\n(b) growing rulers: direct value vs series truncated at the smallest term')
fams = [('lam_2', lambda p: 2 * p, 2), ('lam_0.5', lambda p: p / 2, 0.5), ('lam_1', lambda p: p, 1),
        ('lin_1_1', lambda p: p + 1, 1), ('lin_1_0.5', lambda p: p + mp.mpf(1) / 2, 1), ('lam_1.5', lambda p: 1.5 * p, 1.5)]
for tag, qf, lam in fams:
    try:
        N = [mp.mpf(l) for l in open('ru1_N_%s_K90.txt' % tag)]
    except IOError:
        print('  (missing %s)' % tag); continue
    for p in [4, 8, 12, 16]:
        p = mp.mpf(p); h = 1 / p
        if qf(p) < 1: continue
        val = L8_v(p, qf(p))
        terms = [N[k] * h ** (k + 1) for k in range(len(N))]
        kmin = min(range(2, len(N)), key=lambda k: abs(terms[k]))
        S = 1 + mp.fsum(terms[:kmin])
        print('  %-9s p=%2d  K*=%2d  smallest term=%.2e  |L/8 - S_K*|=%.2e  e^{-pi lam p}=%.2e' % (
            tag, int(p), kmin, float(abs(terms[kmin])), float(abs(val - S)), float(mp.exp(-mp.pi * lam * p))))

print('\n(c) dual ruler q = p/(p-1)')
for p in [mp.mpf(3), mp.mpf(1.5), mp.mpf(7), mp.mpf(2)]:
    q = p / (p - 1)
    ex = mp.gamma(1 + 1 / p) ** 2 / mp.gamma(1 + 2 / p)
    print('  p=%s  L/8 quad = %s  Gamma form = %s  diff=%s' % (p, mp.nstr(L8_v(p, q), 30), mp.nstr(ex, 30), mp.nstr(L8_v(p, q) - ex, 3)))
# series of Gamma(1+h)^2/Gamma(1+2h) = exp(sum_{k>=2} (-1)^k (2-2^k) zeta(k) h^k / k)
K = 80
c = [mp.mpf(0)] * (K + 1)
for k in range(2, K + 1): c[k] = (-1) ** k * (2 - mp.mpf(2) ** k) * mp.zeta(k) / k
f = [mp.mpf(1)] + [mp.mpf(0)] * K
for m in range(1, K + 1):
    f[m] = mp.fsum(j * c[j] * f[m - j] for j in range(1, m + 1)) / m
print('  coefficients D_k of L/8 = sum D_k h^k:')
print('   D_2 = %s (-zeta(2) = %s)' % (mp.nstr(f[2], 25), mp.nstr(-mp.zeta(2), 25)))
print('   D_3 = %s (2 zeta(3) = %s)' % (mp.nstr(f[3], 25), mp.nstr(2 * mp.zeta(3), 25)))
print('   D_4 = %s (-pi^4/40 = %s)' % (mp.nstr(f[4], 25), mp.nstr(-mp.pi ** 4 / 40, 25)))
print('   |D_k| for k=20,40,60,80: %s' % [mp.nstr(abs(f[k]), 8) for k in (20, 40, 60, 80)])
print('   residue check: D_k (-1)^k -> 2*Res = %s' % mp.nstr(-2 * mp.mpf(1), 5))
for p in [mp.mpf(10), mp.mpf(2), mp.mpf(1.25)]:
    h = 1 / p
    print('   p=%s: partial sum (80 terms) - exact = %s' % (p, mp.nstr(mp.fsum(f[k] * h ** k for k in range(K + 1)) - mp.gamma(1 + h) ** 2 / mp.gamma(1 + 2 * h), 3)))
