"""RULQ script 5: tables for result.md from the coefficient files of ru1 (run a, 100-digit J table).
nu_k = N_k / (2 (k-2)! / (pi^{k-1} lambda^k));  F_k = two-term law of Theorem A3;  r_k = (N_k - F_k)/scale * k^2/log^2 k.
Also: |N_k|^{1/k} e pi lambda / k (-> 1), and the q = p^2 growth check."""
import mpmath as mp
mp.mp.dps = 50
fams = [('lam_2', 2, 0), ('lam_0.5', 0.5, 0), ('lin_1_1', 1, 1), ('lin_1_0.5', 1, 0.5), ('lam_1.5', 1.5, 0), ('lam_1', 1, 0)]
y0 = mp.mpc(0, mp.pi)
for tag, L, C in fams:
    L = mp.mpf(L); C = mp.mpf(C)
    try:
        N = [mp.mpf(l) for l in open('ru1_N_%s_K90.txt' % tag)]
    except IOError:
        N = [mp.mpf(l) for l in open('runb/ru1_N_%s_K90.txt' % tag)]
    b0 = C - L; b1 = -C; pre = mp.exp((1 - b0) * y0)
    def F(k):
        br = mp.digamma(k - 2) + mp.euler - mp.log(y0) - b1 * y0 - C / L
        return 2 * mp.re(pre * L ** (-k) * (mp.factorial(k - 2) * y0 ** (1 - k) + L * mp.factorial(k - 3) * y0 ** (2 - k) * br))
    def lead(k):
        return -mp.cos(mp.pi * (L - C) - mp.pi * (k - 1) / 2)
    print('\n== %s  (q = %s p + %s): leading cos(...) per parity: odd k -> %s, even k -> %s' % (
        tag, L, C, mp.nstr(abs(lead(89)), 3), mp.nstr(abs(lead(90)), 3)))
    print('   k    nu_k          F_k/scale     lead    r_k      |N_k|^(1/k) e pi lam/k')
    for k in (30, 31, 50, 51, 70, 71, 88, 89, 90):
        sc = 2 * mp.factorial(k - 2) / (mp.pi ** (k - 1) * L ** k)
        print('  %3d  % .8f  % .8f  % .3f  % .3f   %.5f' % (k, float(N[k] / sc), float(F(k) / sc), float(lead(k)),
              float((N[k] - F(k)) / sc * k ** 2 / mp.log(k) ** 2), float(abs(N[k]) ** (mp.mpf(1) / k) * mp.e * mp.pi * L / k)))
N = [mp.mpf(l) for l in open('runb/ru1_N_sq_K90.txt')]
print('\n== q = p^2:  |N_k|^(2/k) / (k/(2 e pi))  and  log|N_k| / ((k/2) log k)')
for k in (30, 40, 50, 60, 70, 80, 90):
    print('  %3d  %.5f  %.5f' % (k, float(abs(N[k]) ** (mp.mpf(2) / k) / (k / (2 * mp.e * mp.pi))), float(mp.log(abs(N[k])) / (k / 2 * mp.log(k)))))
