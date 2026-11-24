"""RULQR script 3c: q = p^2 (A4). Gevrey-1/2 test |N_k|^{2/k}/(k/(2 e pi)) and comparison with RULQ ru1_N_sq_K90.txt."""
import mpmath as mp
mp.mp.dps = 60
rd = lambda f: [mp.mpf(l.split()[0]) for l in open(f) if l.strip()]
A, B, RU = rd('rr2_N_sq_a.txt'), rd('rr2_N_sq_b.txt'), rd('../RULQ_growing_rulers/ru1_N_sq_K90.txt')
print('k  N_k  rel|a-b|  rel|a-RULQ|  |N_k|^(2/k)/(k/(2 e pi))')
for k in [4, 5, 10, 11, 20, 21, 30, 31, 40, 41, 50, 51, 58, 59, 60]:
    print('%3d %s %8.1e %8.1e %7.3f' % (k, mp.nstr(A[k], 20).rjust(28), float(abs(A[k] - B[k]) / abs(A[k])), float(abs(A[k] - RU[k]) / abs(A[k])),
          float(abs(A[k]) ** (mp.mpf(2) / k) / (k / (2 * mp.e * mp.pi)))))
