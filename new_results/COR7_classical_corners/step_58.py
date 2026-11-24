import mpmath as mp
import sympy as sp
from math import gcd

mp.mp.dps = 80

def sqfree_list(n):
    return [d for d in range(2, n+1) if all(e == 1 for e in sp.factorint(d).values())]

def num_rank_quad(D, a, b, e):
    rD = mp.sqrt(D)
    c1 = (a + b*rD)/e
    if not (0 < c1 < 1):
        return None
    c2 = (a - b*rD)/e
    letters = [(mp.mpf(2), mp.mpf(2)), (c1, c2), (1-c1, 1-c2), (1+c1, 1+c2)]
    def norm(p, q, ee):
        return sp.Rational(p*p - D*q*q, ee*ee)
    Ns = [sp.Integer(4), norm(a, b, e), norm(e-a, -b, e), norm(e+a, b, e)]
    primes = set()
    for N in Ns:
        for p in list(sp.factorint(abs(N.numerator)).keys()) + list(sp.factorint(N.denominator).keys()):
            if p > 2:
                primes.add(p)
    primes = sorted(primes)
    LV = []
    for (x1, x2), N in zip(letters, Ns):
        v = [mp.log(abs(x1)), mp.log(abs(x2))]
        for p in primes:
            eN = sp.factorint(abs(N.numerator)).get(p, 0) - sp.factorint(N.denominator).get(p, 0)
            v.append(eN * mp.log(p))
        LV.append(v)
    A = [row[:] for row in LV]
    R = 0
    tol = mp.mpf(10)**(-50)
    for cc in range(len(A[0])):
        if R >= len(A): break
        p = max(range(R, len(A)), key=lambda r: abs(A[r][cc]))
        if abs(A[p][cc]) < tol: continue
        A[R], A[p] = A[p], A[R]
        for r in range(R+1, len(A)):
            if A[r][cc]:
                f = A[r][cc]/A[R][cc]
                A[r] = [x - f*y for x, y in zip(A[r], A[R])]
        R += 1
    return R

cands = []
BND = 14
for D in sqfree_list(30):
    rD = sp.sqrt(D)
    for e in range(1, BND+1):
        for a in range(-BND, BND+1):
            for b in range(1, BND+1):
                if gcd(gcd(abs(a), b), e) != 1: continue
                val = (a + b*rD)/e
                if not (0 < val < 1): continue
                r = num_rank_quad(D, a, b, e)
                if r is not None and r <= 2:
                    cands.append((D, a, b, e, r))
print("candidates with numerical rank <= 2:", len(cands))
for D, a, b, e, r in cands:
    print(f"  Q(sqrt{D}): c0 = ({a}+{b}*sqrt{D})/{e}  ~ {mp.nstr((a+b*mp.sqrt(D))/e, 12)}  rank {r}")
