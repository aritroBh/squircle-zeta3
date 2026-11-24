# 20_ndim_check.py -- the n-dimensional superball |x_1|^p+...+|x_n|^p = 1.
# Orthant hypersurface area (same Gaussian argument as for n=3):
#   A_n(eps) = sum_{m=1}^n (-1)^(m-1) C(n,m) Gamma(1+eps)^(n-m) eps^(m-1) H_m(eps) / Gamma(1+(n-1)eps),
#   S_n(p) = 2^n A_n(1/p).
# Checks: p=2 (unit sphere S^{n-1}: area 2 pi^{n/2}/Gamma(n/2)) and p=1 (cross-polytope: 2^n sqrt(n)/(n-1)!).
import mpmath as mp
from superlib import H
mp.mp.dps = 30
def A_n(n, eps):
    eps = mp.mpf(eps); G1 = mp.gamma(1+eps)
    tot = 0
    for m in range(1, n+1):
        Hm = 1 if m == 1 else H(m, eps)
        tot += (-1)**(m-1)*mp.binomial(n, m)*G1**(n-m)*eps**(m-1)*Hm
    return tot/mp.gamma(1+(n-1)*eps)
if __name__ == "__main__":
  for n in (2, 3, 4, 5):
      S2 = 2**n*A_n(n, mp.mpf(1)/2); S2ex = 2*mp.pi**(mp.mpf(n)/2)/mp.gamma(mp.mpf(n)/2)
      S1 = 2**n*A_n(n, 1); S1ex = 2**n*mp.sqrt(n)/mp.factorial(n-1)
      print("n=%d  p=2: %s vs %s (diff %s)   p=1: %s vs %s (diff %s)" % (n, mp.nstr(S2, 20), mp.nstr(S2ex, 20), mp.nstr(S2-S2ex, 3),
            mp.nstr(S1, 20), mp.nstr(S1ex, 20), mp.nstr(S1-S1ex, 3)))
