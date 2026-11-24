from fractions import Fraction as Fr
from math import gcd
import sympy as sp

def letter_rank_rational(c):
    """rank of <2, c, 1-c, 1+c> in Q* via exponent vectors over the union support."""
    letters = [Fr(2), c, 1-c, 1+c]
    primes = set()
    facs = []
    for L in letters:
        f = {}
        for q in (abs(L.numerator), L.denominator):
            for p, e in sp.factorint(q).items():
                f[p] = f.get(p, 0) + (e if q == abs(L.numerator) else -e)
                primes.add(p)
        facs.append(f)
    primes = sorted(primes)
    M = [[f.get(p, 0) for p in primes] for f in facs]
    # exact rank
    M = [row[:] for row in M]; R = 0
    for cc in range(len(M[0])):
        p = next((r for r in range(R, len(M)) if M[r][cc]), None)
        if p is None: continue
        M[R], M[p] = M[p], M[R]
        q = M[R][cc]; M[R] = [x / q for x in M[R]]
        for r in range(len(M)):
            if r != R and M[r][cc]:
                M[r] = [a - b*M[R][cc] for a, b in zip(M[r], M[R])]
        R += 1
    return R

sols = []
B = 5000
for b in range(2, B+1):
    for a in range(1, b):
        if gcd(a, b) != 1: continue
        r = letter_rank_rational(Fr(a, b))
        if r <= 2:
            sols.append((a, b, r))
print("rational c0 in (0,1) with denom <= 5000 and letter-rank <= 2:")
for a, b, r in sols:
    print(f"  {a}/{b}  rank {r}")
