"""flat_symbolic_proof.py -- symbolic-n proof of the zeta content of the flat-corner Taylor coefficients.

mu_{k,n} = -(1/(n(2n-1))) * sum_{(j,i,l)} e^{(k)}_{j,i,l}(n) * (n)_j/(2n)_j * D_{i,l}(n+j, n),
where [u^k]E_n = sum e_{j,i,l}(n) y^j a^i b^l (a = log y, b = log(1-y)) and D_{i,l} = d_p^i d_q^l B / B.
D_{i,l} is the moment polynomial of the cumulants kappa_{r,s} = d_p^r d_q^s log B(p,q):
  kappa_{r,0} = psi^{(r-1)}(p) - psi^{(r-1)}(p+q), kappa_{0,s} = psi^{(s-1)}(q) - psi^{(s-1)}(p+q),
  kappa_{r,s} = -psi^{(r+s-1)}(p+q)  (r,s >= 1).
At integers psi^{(m)}(N) = (-1)^(m+1) m! (zeta(m+1) - H^{(m+1)}_{N-1}) (m >= 1), so the pure kappa_{r,0},
kappa_{0,s} are RATIONAL (zeta cancels), the mixed ones carry zeta(r+s), and psi(N)+gamma = H_{N-1}.
Consequences proved here for all n (symbolic n):
 (P1) mu_{2,n} in Q + Q zeta(2): monomials have i+l <= 2, only kappa_{1,1} = -psi'(p+q) carries zeta(2).
      zeta(2)-coefficient = (1/(n(2n-1))) * sum_j e^{(2)}_{j,1,1}(n) (n)_j/(2n)_j  ==  1/(4n^2-1).
      In particular NO zeta(3) and no log 2 in any mu_{2,n}.
 (P2) mu_{3,n} in Q + Q zeta(2) + Q zeta(3); zeta(3) only from kappa_{2,1}, kappa_{1,2} (each = 2 zeta(3) + rational):
      zeta(3)-coefficient = -(2/(n(2n-1))) * sum_j [e^{(3)}_{j,2,1} + e^{(3)}_{j,1,2}](n) (n)_j/(2n)_j == -1/(4n^2-1).
 (P3) same check for k = 4, 5: coefficient of the single zeta(k) (from kappa_{r,s}, r+s = k, r,s>=1; each
      kappa carries (-1)^k (k-1)! zeta(k) ... with sign from -psi^{(k-1)}) -- printed as a rational function of n.
"""
import sympy as sp
from math import factorial

u, y, a, b, n = sp.symbols('u y a b n')
import sys
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3
A = -(n - 1)*a + n*b
l = b - a
logE = u*A + (1 - 2*n)*sp.log(1 + y*(sp.exp(u*l) - 1))
E_ser = sp.series(sp.exp(sp.series(logE, u, 0, KMAX + 1).removeO()), u, 0, KMAX + 1).removeO()


def poch(x, j):
    return sp.Mul(*[x + t for t in range(j)]) if j > 0 else sp.Integer(1)


for k in range(2, KMAX + 1):
    P = sp.Poly(sp.expand(E_ser.coeff(u, k)), y, a, b)
    # coefficient of the single zeta(k): from kappa_{r,s} with r+s=k, r,s>=1 appearing linearly in D_{r,s}
    # (D_{i,l} with i+l = k contains kappa_{i,l} with coefficient 1; products of lower cumulants carry only
    #  lower zetas).  kappa_{r,s} = -psi^{(k-1)}(p+q) = -(-1)^k (k-1)! (zeta(k) - H) -> zeta(k) coeff
    zc = -(-1)**k*factorial(k - 1)
    tot = 0
    for (j, i, ll), co in P.terms():
        if i + ll == k and i >= 1 and ll >= 1:
            tot += co*poch(n, j)/poch(2*n, j)*zc
    coef = sp.factor(sp.simplify(-tot/(n*(2*n - 1))))
    print(f'k={k}: coefficient of the single zeta({k}) in mu_(k,n) = {coef}')
    print(f'      times (4n^2-1): {sp.simplify(coef*(4*n**2 - 1))}')
    if k == 2:
        # also: no monomial with i+l >= 3 (so no zeta(3) possible) and the maximal total log degree
        degs = sorted({i + ll for (j, i, ll) in P.monoms()})
        print('      total log-degrees of monomials in [u^2]E:', degs, '(<= 2 => only psi, psi\' => Q + Q zeta(2))')
