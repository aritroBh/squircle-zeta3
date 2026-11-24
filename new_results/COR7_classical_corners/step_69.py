import sys, time
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp
from cor7_pt import QuadField, QuadPoint, factorint

mp.mp.dps = 80
H1 = [2,3,5,6,7,11,13,17,19,21,23,29,33,37,41]  # real quadratic fields with h=1, D <= 50

def sqfree(D):
    return all(e == 1 for e in factorint(D).values())

def lv_rank(D, a, b, e):
    """cheap filter: rank of LV lattice (embeddings + norm valuations); <= true rank."""
    rD = mp.sqrt(D)
    c1 = (mp.mpf(a) + mp.mpf(b)*rD)/e
    if not (0 < c1 < 1): return None
    c2 = (mp.mpf(a) - mp.mpf(b)*rD)/e
    letters = [(mp.mpf(2), mp.mpf(2)), (c1, c2), (1-c1, 1-c2), (1+c1, 1+c2)]
    Ns = [(4,1), (a*a - D*b*b, e*e), ((e-a)*(e-a) - D*b*b, e*e), ((e+a)*(e+a) - D*b*b, e*e)]
    primes = set()
    facN = []
    for num, den in Ns:
        f = {}
        if abs(num) > 1:
            for p, ee in factorint(abs(num)).items(): f[p] = f.get(p, 0) + ee
        for p, ee in factorint(den).items(): f[p] = f.get(p, 0) - ee
        facN.append(f)
        primes |= set(f.keys())
    primes = sorted(primes)
    LV = []
    for (x1, x2), f in zip(letters, facN):
        v = [mp.log(abs(x1)), mp.log(abs(x2))]
        for p in primes:
            v.append(f.get(p, 0) * mp.log(p))
        LV.append(v)
    A = [row[:] for row in LV]; R = 0
    tol = mp.mpf(10)**(-45)
    for cc in range(len(A[0])):
        if R >= len(A): break
        p = max(range(R, len(A)), key=lambda r: abs(A[r][cc]))
        if abs(A[p][cc]) < tol: continue
        A[R], A[p] = A[p], A[R]
        for r in range(R+1, len(A)):
            if A[r][cc]:
                ff = A[r][cc]/A[R][cc]
                A[r] = [x - ff*y for x, y in zip(A[r], A[R])]
        R += 1
    return R, primes

t0 = time.time()
BND = 16
survivors = []
skipped_gen = 0
QFs = {}
for D in H1:
    QFs[D] = QuadField(D)
    for e in range(1, BND+1):
        for a in range(-BND, BND+1):
            for b in range(-BND, BND+1):
                if b == 0: continue
                if gcd(gcd(abs(a), abs(b)), e) != 1: continue
                r = lv_rank(D, a, b, e)
                if r is None or r[0] > 2: continue
                primes = [p for p in r[1]]
                P = QuadPoint(QFs[D], f"D{D} ({a},{b})/{e}", (Fr(a, e), Fr(b, e)), primes=primes)
                rr, labs, vecs = P.letter_span_rank()
                if rr is not None and rr <= 2:
                    survivors.append((D, a, b, e, rr, vecs))
print(f"{time.time()-t0:.0f}s; exact rank<=2 survivors: {len(survivors)}")
for D, a, b, e, rr, vecs in survivors:
    print(f"Q(sqrt{D}): c0 = ({a}{'+' if b>=0 else ''}{b}s{D})/{e}  rank {rr}   vecs {vecs}")
import json
json.dump([[D, a, b, e, rr] for D, a, b, e, rr, _ in survivors], open('census_rank2.json', 'w'))
