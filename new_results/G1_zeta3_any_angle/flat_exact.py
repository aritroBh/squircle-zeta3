"""flat_exact.py -- EXACT Taylor coefficients mu_{k,n} of M_k(alpha) in powers of c^2 = cos^2(alpha/2).

Derivation (result.md Sec. 4).  With y = g/(1+g) and the g -> 1/g symmetry of each term,
  M_k(alpha) = sum_{n>=1} mu_{k,n} c^(2n),
  mu_{k,n} = -(1/(n(2n-1))) * [u^k] ( J~_n(u) / B(n,n) ),
  J~_n(u) = int_0^1 y^(n-1)(1-y)^(n-1) E_n(y;u) dy,
  E_n = exp(u A) (1 + y(e^(u l) - 1))^(1-2n),  A = -(n-1) log y + n log(1-y),  l = log(1-y) - log y.
[u^k]E_n is a polynomial in y, a = log y, b = log(1-y); each monomial y^j a^i b^l integrates to
  d_p^i d_q^l B(p,q) at p = n+j, q = n,
i.e. B(n,n) (n)_j/(2n)_j times a Bell polynomial in the cumulants of log B (polygammas at integers).
Polygammas at integers: psi(k)+gamma = H_{k-1}, psi'(k) = z2 - H^(2)_{k-1}, psi''(k) = -2 z3 + 2 H^(3)_{k-1},
psi'''(k) = 6 z4 - 6 H^(4)_{k-1}.  Hence mu_{k,n} is EXACTLY a polynomial in z2, z3, z4 with rational
coefficients (no log 2, no other constants).  This script computes it exactly (fractions) and prints
  mu_{k,n} = r + q2 z2 + q3 z3 (+ q4 z4 + q22 z2^2)
for k <= KMAX, n <= NMAX, compares with the PSLQ results of flat_expansion.py, and saves a JSON table.
Usage: python flat_exact.py NMAX KMAX
"""
import sys, json
from fractions import Fraction as Fr
import sympy as sp

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3

u, y, a, b, n = sp.symbols('u y a b n')
A = -(n - 1)*a + n*b
l = b - a
# series of E in u up to u^KMAX (exact, symbolic n)
logE = u*A + (1 - 2*n)*sp.log(1 + y*(sp.exp(u*l) - 1))
logE_ser = sp.series(logE, u, 0, KMAX + 1).removeO()
E_ser = sp.series(sp.exp(logE_ser), u, 0, KMAX + 1).removeO()
Ek = [sp.Poly(sp.expand(E_ser.coeff(u, k)), y, a, b) for k in range(KMAX + 1)]

# ---- exact arithmetic in Q[z2,z3,z4] : dict {(e2,e3,e4): Fraction}
def padd(P, Q, s=1):
    R = dict(P)
    for k, v in Q.items():
        R[k] = R.get(k, 0) + s*v
        if R[k] == 0:
            del R[k]
    return R

def pmul(P, Q):
    R = {}
    for k1, v1 in P.items():
        for k2, v2 in Q.items():
            k = tuple(x + y_ for x, y_ in zip(k1, k2))
            R[k] = R.get(k, 0) + v1*v2
            if R[k] == 0:
                del R[k]
    return R

def const(x):
    return {(0, 0, 0): Fr(x)} if x != 0 else {}

Z2 = {(1, 0, 0): Fr(1)}; Z3 = {(0, 1, 0): Fr(1)}; Z4 = {(0, 0, 1): Fr(1)}

_H = {}
def H(r, m):  # generalized harmonic number H^{(r)}_m
    key = (r, m)
    if key not in _H:
        _H[key] = sum((Fr(1, i**r) for i in range(1, m + 1)), Fr(0))
    return _H[key]

def psi0(k):  return const(H(1, k - 1))                       # psi(k) + gamma
def psi1(k):  return padd(Z2, const(H(2, k - 1)), -1)          # psi'(k)
def psi2(k):  return padd(pmul(const(-2), Z3), const(2*H(3, k - 1)))
def psi3(k):  return padd(pmul(const(6), Z4), const(-6*H(4, k - 1)))
PSI = [psi0, psi1, psi2, psi3]

def cumulant(i, l, p, q):
    """d_p^i d_q^l log B(p,q) for i+l >= 1 (exact, in Q[z2,z3,z4])."""
    m = i + l
    f = PSI[m - 1]
    out = pmul(const(-1), f(p + q))
    if l == 0:
        out = padd(out, f(p))
    if i == 0:
        out = padd(out, f(q))
    return out

def bell(i, l, p, q):
    """d_p^i d_q^l B / B : sum over set partitions of the multiset {p^i q^l} of products of cumulants."""
    from math import comb
    memo = {}

    def DD(i, l):
        if (i, l) in memo:
            return memo[(i, l)]
        if i == 0 and l == 0:
            res = const(1)
        elif i > 0:
            # differentiate wrt p once more: D_{i,l} = sum_{r=0}^{i-1} sum_{s=0}^{l} C(i-1,r) C(l,s) k_{r+1,s} D_{i-1-r, l-s}
            res = {}
            for r in range(i):
                for s in range(l + 1):
                    term = pmul(const(comb(i - 1, r)*comb(l, s)), pmul(cumulant(r + 1, s, p, q), DD(i - 1 - r, l - s)))
                    res = padd(res, term)
        else:
            res = {}
            for s in range(l):
                term = pmul(const(comb(l - 1, s)), pmul(cumulant(0, s + 1, p, q), DD(0, l - 1 - s)))
                res = padd(res, term)
        memo[(i, l)] = res
        return res
    return DD(i, l)

def poch_ratio(nv, j):
    r = Fr(1)
    for t in range(j):
        r *= Fr(nv + t, 2*nv + t)
    return r

def mu(k, nv):
    tot = {}
    for (j, i, l), co in Ek[k].terms():
        cr = sp.Rational(sp.expand(co.subs(n, nv))); cval = Fr(int(cr.p), int(cr.q))
        term = pmul(const(cval*poch_ratio(nv, j)), bell(i, l, nv + j, nv))
        tot = padd(tot, term)
    return pmul(const(Fr(-1, nv*(2*nv - 1))), tot)

def fmt(P):
    names = {(0, 0, 0): '1', (1, 0, 0): 'z2', (0, 1, 0): 'z3', (0, 0, 1): 'z4', (2, 0, 0): 'z2^2'}
    return ' + '.join(f'({v}){names.get(k, str(k))}' for k, v in sorted(P.items()))

if __name__ == '__main__':
    table = {}
    for k in range(KMAX + 1):
        table[k] = {}
        for nv in range(1, NMAX + 1):
            P = mu(k, nv)
            table[k][nv] = {str(kk): str(vv) for kk, vv in P.items()}
            if nv <= 6:
                print(f'mu_{{{k},{nv}}} = {fmt(P)}', flush=True)
        # the canonical single-zeta coefficient pattern
        key = {2: (1, 0, 0), 3: (0, 1, 0)}.get(k)
        if key:
            ok = all(Fr(table[k][nv].get(str(key), '0')) * (4*nv*nv - 1) == (1 if k == 2 else -1)
                     for nv in range(1, NMAX + 1))
            print(f'k={k}: coefficient of {"z2" if k==2 else "z3"} equals {"+" if k==2 else "-"}1/(4n^2-1) for all n<={NMAX}: {ok}')
        if k == 2:
            print('k=2: no z3 term for any n <=', NMAX, ':', all(str((0, 1, 0)) not in table[2][nv] for nv in table[2]))
    with open(f'flat_exact_N{NMAX}_K{KMAX}.json', 'w') as f:
        json.dump(table, f)
    print('saved', f'flat_exact_N{NMAX}_K{KMAX}.json')
