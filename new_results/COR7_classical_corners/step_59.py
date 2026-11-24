import mpmath as mp
from math import gcd
import time

mp.mp.dps = 80

def factorint_small(n):
    f = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def sqfree_list(n):
    out = []
    for d in range(2, n+1):
        f = factorint_small(d)
        if all(e == 1 for e in f.values()):
            out.append(d)
    return out

def num_rank_quad(D, a, b, e):
    rD = mp.sqrt(D)
    c1 = mp.mpf(a + 0)/e + mp.mpf(b)*rD/e
    if not (0 < c1 < 1):
        return None
    c2 = mp.mpf(a)/e - mp.mpf(b)*rD/e
    letters = [(mp.mpf(2), mp.mpf(2)), (c1, c2), (1-c1, 1-c2), (1+c1, 1+c2)]
    Ns = [(4,1), (a*a - D*b*b, e*e), ((e-a)*(e-a) - D*b*b, e*e), ((e+a)*(e+a) - D*b*b, e*e)]
    primes = set()
    facN = []
    for num, den in Ns:
        fn = factorint_small(abs(num)) if num != 0 else {}
        fd = factorint_small(den)
        f = dict(fn)
        for p, ee in fd.items():
            f[p] = f.get(p, 0) - ee
        facN.append(f)
        primes |= {p for p in f if p > 2}
    primes = sorted(primes)
    LV = []
    for (x1, x2), f in zip(letters, facN):
        v = [mp.log(abs(x1)), mp.log(abs(x2))]
        for p in primes:
            v.append(f.get(p, 0) * mp.log(p))
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
                ff = A[r][cc]/A[R][cc]
                A[r] = [x - ff*y for x, y in zip(A[r], A[R])]
        R += 1
    return R

t0 = time.time()
cands = []
BND = 16
for D in sqfree_list(50):
    for e in range(1, BND+1):
        for a in range(-BND, BND+1):
            for b in range(-BND, BND+1):
                if b == 0: continue
                if gcd(gcd(abs(a), abs(b)), e) != 1: continue
                r = num_rank_quad(D, a, b, e)
                if r is not None and r <= 2:
                    cands.append((D, a, b, e, r))
print(f"{time.time()-t0:.0f}s, candidates: {len(cands)}")
for D, a, b, e, r in cands:
    print(f"  Q(sqrt{D}): c0 = ({a}{'+' if b>=0 else ''}{b}*sqrt{D})/{e}  ~ {mp.nstr((a+b*mp.sqrt(D))/e, 12)}  rank {r}")
