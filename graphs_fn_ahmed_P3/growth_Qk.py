"""growth_Qk.py -- size of the expansion coefficients Q_k (at Xi = 0) for k <= 30:
factorial growth => the expansion of Theorem 2 is divergent (asymptotic only) for generic f;
for f = x it is convergent (Ahmed's exact formula: radius |mu| < 1)."""
import mpmath as mp
from fnlib import TESTS, jets, solve_T, moments_xi
mp.mp.dps = 80
K = 30
for name in ['x', 'expm1', 'mobius']:
    A, B = jets(TESTS[name], K + 2)
    a = A[1]
    # Q_k(Xi) at Xi = 0 via compact form: sum_i theta_{i,k-i} m_i (central moments) + T_k ; use engine P_k at L = -C
    from fnlib import expansion_polys
    P = expansion_polys(A, B, K)
    C = mp.log(a) + 1 - mp.log(2)
    vals = [mp.polyval(P[k][::-1], -C) for k in range(1, K + 1)]   # Q_k(0)
    print(name)
    for k in [5, 10, 15, 20, 25, 30]:
        q = vals[k - 1]
        print('   k=%2d  Q_k(0)=%-12s  |Q_k(0)|^(1/k)=%-10s  |Q_k(0)|/k! = %s' % (k, mp.nstr(q, 5), mp.nstr(abs(q)**(mp.mpf(1)/k), 5), mp.nstr(abs(q)/mp.factorial(k), 5)))
