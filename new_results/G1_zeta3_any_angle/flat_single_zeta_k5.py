"""flat_single_zeta_k5.py -- single-zeta coefficient of the flat-corner Taylor coefficients mu_{k,n} for k = 4, 5
(numeric n = 1..NMAX, exact rationals).  The coefficient of the single zeta(k) comes only from the mixed cumulants
kappa_{r,s} = -psi^{(k-1)}(p+q), r+s = k, r,s >= 1, each = -(-1)^k (k-1)! zeta(k) + rational, entering D_{r,s}
linearly (products of lower cumulants carry products of lower zetas).  Prints (4n^2-1) * coefficient."""
import sys
from fractions import Fraction as Fr
from math import factorial
import sympy as sp

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
u, y, a, b = sp.symbols('u y a b')
for k in (2, 3, 4, 5):
    zc = Fr(-(-1)**k*factorial(k - 1))
    row = []
    for n in range(1, NMAX + 1):
        A = -(n - 1)*a + n*b
        l = b - a
        # [u^k]E = sum_{i+j=k} A^i/i! * sum_{m<=j} binom(1-2n,m) y^m m! S2(j,m) l^j/j!   (Stirling 2nd kind)
        Ek = 0
        for i in range(k + 1):
            j = k - i
            inner = sum(sp.binomial(1 - 2*n, m)*y**m*sp.factorial(m)*sp.functions.combinatorial.numbers.stirling(j, m)
                        for m in range(0, j + 1)) if j > 0 else 1
            Ek += A**i/sp.factorial(i)*inner*l**j/sp.factorial(j)
        E = Ek
        P = sp.Poly(sp.expand(E), y, a, b)
        tot = Fr(0)
        for (j, i, ll), co in P.terms():
            if i + ll == k and i >= 1 and ll >= 1:
                r = Fr(1)
                for t in range(j):
                    r *= Fr(n + t, 2*n + t)
                cr = sp.Rational(co)
                tot += Fr(int(cr.p), int(cr.q))*r*zc
        coef = -tot/(n*(2*n - 1))
        row.append(coef*(4*n*n - 1))
    print(f'k={k}: (4n^2-1) * [single zeta({k}) coefficient of mu_(k,n)], n=1..{NMAX}:', [str(x) for x in row])
