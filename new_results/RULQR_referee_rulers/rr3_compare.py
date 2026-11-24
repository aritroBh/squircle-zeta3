"""RULQR script 3: N_k cross-checks (run a vs run b; vs RULQ ru1 files; vs PIPDIVR ball-arithmetic N_k for q=p)
and comparison with the two-term law of RULQ Thm A3(ii).  Usage: rr3_compare.py"""
import mpmath as mp
mp.mp.dps = 60
R = '../RULQ_growing_rulers/'
fams = {  # name: (lam, c, b1, RULQ file)
    'lam2': (2, 0, 0, R + 'ru1_N_lam_2_K90.txt'),
    'lamhalf': (mp.mpf(1) / 2, 0, 0, R + 'ru1_N_lam_0.5_K90.txt'),
    'p1': (1, 1, -1, R + 'ru1_N_lin_1_1_K90.txt'),
    'p': (1, 0, 0, R + 'ru1_N_lam_1_K90.txt')}
rd = lambda f: [mp.mpf(l.split()[0]) for l in open(f) if l.strip()]
y0 = mp.mpc(0, mp.pi)
for name, (lam, c, b1, rf) in fams.items():
    lam = mp.mpf(lam); c = mp.mpf(c)
    A, B, RU = rd('rr2_N_%s_a.txt' % name), rd('rr2_N_%s_b.txt' % name), rd(rf)
    b0 = c - lam
    print('=== family', name, ' lambda =', lam, ' c =', c, ' b1 =', b1)
    print(' k   N_k (run a)                 rel|a-b|   rel|a-RULQ|  nu_k        two-term F_k  lead      r_k=(nu-F)k^2/log^2k  |N_k|^(1/k) e pi lam/k')
    worst_ab = worst_ru = 0
    for k in range(2, len(A)):
        ab = abs(A[k] - B[k]) / abs(A[k]); ru = abs(A[k] - RU[k]) / abs(A[k])
        if k <= 45: worst_ab = max(worst_ab, ab); worst_ru = max(worst_ru, ru)
        if k < 4: continue
        scale = 2 * mp.factorial(k - 2) * mp.pi ** (1 - k) * lam ** (-k)
        two = 2 * mp.re(mp.exp((1 - b0) * y0) * lam ** (-k) * (mp.factorial(k - 2) * y0 ** (1 - k)
              + lam * mp.factorial(k - 3) * y0 ** (2 - k) * (mp.digamma(k - 2) + mp.euler - mp.log(y0) - b1 * y0 - c / lam)))
        nu, F = A[k] / scale, two / scale
        lead = mp.cos(mp.pi * (1 + lam - c) - mp.pi * (k - 1) / 2)
        rk = (nu - F) * k ** 2 / mp.log(k) ** 2
        g = abs(A[k]) ** (mp.mpf(1) / k) * mp.e * mp.pi * lam / k
        if k in (10, 11, 20, 21, 30, 31, 38, 39, 40, 41, 44, 45):
            print('%3d %s %9.1e %9.1e  %10.6f  %10.6f  %6.3f  %8.3f  %6.3f' % (k, mp.nstr(A[k], 22).rjust(28), float(ab), float(ru),
                  float(nu), float(F), float(lead), float(rk), float(g)))
    print(' worst rel |a-b| for k<=45: %.1e ; worst rel |a - RULQ| for k<=45: %.1e' % (worst_ab, worst_ru))
    if name == 'p':
        PR = [mp.mpf(l.split('[')[1].split('+/-')[0]) for l in open('../PIPDIVR_referee_divergence/rf1_N_K151_p700_S400.txt')]
        print(' worst rel |a - PIPDIVR ball N_k| (k<=45): %.1e' % max(abs(A[k] - PR[k]) / abs(PR[k]) for k in range(46)))
