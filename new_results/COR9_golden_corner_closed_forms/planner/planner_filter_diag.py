# Diagnose the survivor-count difference: COR9 step_57 row 0 ("2") subtracts 2 v_p(e) (the "f is not ({2:2})" test is always True),
# so it filters {2/e, c, 1-c, 1+c} instead of {2, c, 1-c, 1+c}. Recompute both filters on the same candidates.
import math
from fractions import Fraction as Fr
from collections import Counter
LIM = 400000
spf = list(range(LIM + 1))
for i in range(2, int(LIM ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, LIM + 1, i):
            if spf[j] == j: spf[j] = i
def fac(n):
    n = abs(n); out = {}
    while n > 1:
        p = spf[n]; k = 0
        while n % p == 0: n //= p; k += 1
        out[p] = k
    return out
def qrank(rows):
    M = [[Fr(x) for x in r] for r in rows]; rk = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
good, buggy = set(), set()
for D in [D for D in range(2, 201) if all(v == 1 for v in fac(D).values())]:
    sD = math.sqrt(D)
    for e in range(1, 41):
        fe = fac(e)
        for a in range(-40, 41):
            for b in range(-40, 41):
                if b == 0 or math.gcd(math.gcd(abs(a), abs(b)), e) != 1 or not (0 < (a + b * sD) / e < 1): continue
                f = [fac(a*a - b*b*D), fac((e-a)**2 - b*b*D), fac((e+a)**2 - b*b*D)]
                P = sorted(set().union(*f) | set(fe) | {2})
                rows = [[g.get(p, 0) - 2 * fe.get(p, 0) for p in P] for g in f]
                r2 = [2 if p == 2 else 0 for p in P]
                r2b = [(2 if p == 2 else 0) - 2 * fe.get(p, 0) for p in P]
                if qrank([r2] + rows) <= 2: good.add((D, a, b, e))
                if qrank([r2b] + rows) <= 2: buggy.add((D, a, b, e))
print("correct filter:", len(good), " COR9-style filter:", len(buggy), " in correct but not COR9-style:", len(good - buggy), " in COR9-style only:", len(buggy - good))
print("dropped by COR9-style filter (first 10):", sorted(good - buggy)[:10])
