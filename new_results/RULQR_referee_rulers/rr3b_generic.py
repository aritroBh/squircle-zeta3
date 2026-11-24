"""RULQR script 3b: generic-phase family q = p + 1/3 + 1/p (lam=1, c=1/3, b1=2/3; both parities lead) against the
two-term law, and the bounded ruler q = 2 + 1/p (A2: coefficients must grow at most geometrically)."""
import mpmath as mp
mp.mp.dps = 60
rd = lambda f: [mp.mpf(l.split()[0]) for l in open(f) if l.strip()]
y0 = mp.mpc(0, mp.pi)
lam, c, b1 = mp.mpf(1), mp.mpf(1) / 3, mp.mpf(2) / 3
A, B = rd('rr2_N_gen_a.txt'), rd('rr2_N_gen_b.txt')
print('q = p + 1/3 + 1/p.  k, rel|a-b|, nu_k, two-term F_k, lead cos, r_k=(nu-F)k^2/log^2k')
for k in range(6, 51):
    scale = 2 * mp.factorial(k - 2) * mp.pi ** (1 - k)
    two = 2 * mp.re(mp.exp((1 - (c - lam)) * y0) * (mp.factorial(k - 2) * y0 ** (1 - k)
          + mp.factorial(k - 3) * y0 ** (2 - k) * (mp.digamma(k - 2) + mp.euler - mp.log(y0) - b1 * y0 - c)))
    nu, F = A[k] / scale, two / scale
    if k % 4 in (0, 1) or k > 44:
        print('%3d %8.1e %10.6f %10.6f %7.4f %8.3f' % (k, float(abs(A[k] - B[k]) / abs(A[k])), float(nu), float(F),
              float(mp.cos(mp.pi * (1 + lam - c) - mp.pi * (k - 1) / 2)), float((nu - F) * k ** 2 / mp.log(k) ** 2)))
A, B = rd('rr2_N_bdd_a.txt'), rd('rr2_N_bdd_b.txt')
print('bounded ruler q = 2 + 1/p: k, N_k, rel|a-b|, |N_k|^(1/k)')
for k in list(range(0, 8)) + [10, 20, 30, 40, 45, 50]:
    print('%3d %s %8.1e %8.4f' % (k, mp.nstr(A[k], 25).rjust(32), float(abs(A[k] - B[k]) / abs(A[k])), float(abs(A[k]) ** (mp.mpf(1) / max(k, 1)))))
