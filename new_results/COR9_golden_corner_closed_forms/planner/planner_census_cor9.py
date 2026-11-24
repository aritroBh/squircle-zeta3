# Planner census check of COR9 Task 3, independent of COR9's code.
# Candidates c = (a + b sqrt D)/e, D squarefree in [2,200], |a|,|b| <= 40, 1 <= e <= 40, b != 0, gcd(|a|,|b|,e) = 1, 0 < c < 1.
# Filter: rank of the rational-prime norm-valuation matrix of {2, c, 1-c, 1+c} <= 2 (necessary for group rank <= 2).
# Exact group rank of survivors: rank(prime-ideal valuation matrix V) + delta, delta = 1 iff the units in Gamma are not all +-1
# (decided by log|sigma_1| of a Z-basis of the left kernel of V; a non-torsion unit has |log| >= log(golden ratio)).
import math, sys
from fractions import Fraction as Fr
from collections import Counter
from sympy import sqrt_mod, Matrix, legendre_symbol
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
    M = [[Fr(x) for x in r] for r in rows]; rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
def sqfree(n): return all(v == 1 for v in fac(n).values())
Ds = [D for D in range(2, 201) if sqfree(D)]
cands = []; surv = []
for D in Ds:
    sD = math.sqrt(D)
    for e in range(1, 41):
        for a in range(-40, 41):
            for b in range(-40, 41):
                if b == 0 or math.gcd(math.gcd(abs(a), abs(b)), e) != 1: continue
                c = (a + b * sD) / e
                if not (0 < c < 1): continue
                cands.append((D, a, b, e))
                # elements (A, B, E): value (A + B sqrt D)/E
                els = [(2, 0, 1), (a, b, e), (e - a, -b, e), (e + a, b, e)]
                vecs = []
                for A, B, E in els:
                    f = Counter(fac(A * A - B * B * D)); 
                    for p, k in fac(E).items(): f[p] -= 2 * k
                    vecs.append(f)
                primes = sorted(set().union(*vecs))
                if qrank([[v.get(p, 0) for p in primes] for v in vecs]) <= 2: surv.append((D, a, b, e))
print("squarefree D:", len(Ds), " candidates:", len(cands), " filter survivors:", len(surv))
K = 80
def vp_int(n, p):
    if n == 0: return 10 ** 9
    k = 0
    while n % p == 0: n //= p; k += 1
    return k
def ideal_vals(A, B, E, D, primes):
    # returns {(p, tag): v_P(x)} for x = (A + B sqrt D)/E, prime ideals above the given rational primes
    out = {}
    N = A * A - B * B * D
    for p in primes:
        vN = vp_int(N, p) - 2 * vp_int(E, p)
        if p == 2: typ = 'ram' if D % 4 in (2, 3) else ('split' if D % 8 == 1 else 'inert')
        else: typ = 'ram' if D % p == 0 else ('split' if legendre_symbol(D % p, p) == 1 else 'inert')
        if typ == 'ram': out[(p, 0)] = vN
        elif typ == 'inert': assert vN % 2 == 0; out[(p, 0)] = vN // 2
        else:
            mod = p ** K
            r = sqrt_mod(D, mod)
            prec = K - 1 if p == 2 else K
            for tag, rr in ((0, r), (1, -r)):
                t = (A + B * rr) % (p ** prec)
                v = vp_int(t, p) if t != 0 else prec
                assert v < prec - 5, "precision"
                out[(p, tag)] = v - vp_int(E, p)
            assert out[(p, 0)] + out[(p, 1)] == vN
    return out
def nullspace_int(rows):
    M = Matrix(rows).T  # columns = elements; left kernel of rows = kernel of M
    ns = M.nullspace()
    vecs = []
    for v in ns:
        den = math.lcm(*[int(x.q) for x in v]); vecs.append([int(x * den) for x in v])
    return vecs
dist = Counter(); rank2 = []
for D, a, b, e in surv:
    els = [(2, 0, 1), (a, b, e), (e - a, -b, e), (e + a, b, e)]
    primes = sorted(set().union(*[set(fac(A * A - B * B * D)) | set(fac(E)) for A, B, E in els]) | {2})
    vals = [ideal_vals(A, B, E, D, primes) for A, B, E in els]
    keys = sorted(set().union(*vals))
    V = [[v.get(k, 0) for k in keys] for v in vals]
    rV = Matrix(V).rank()
    logs = [math.log(abs((A + B * math.sqrt(D)) / E)) for A, B, E in els]
    delta = 0
    for m in nullspace_int(V):
        s = sum(mi * li for mi, li in zip(m, logs))
        if abs(s) > 0.4: delta = 1
        else: assert abs(s) < 1e-9, ("ambiguous", D, a, b, e, s)
    rk = rV + delta; dist[rk] += 1
    if rk <= 2: rank2.append((D, a, b, e, (a + b * math.sqrt(D)) / e, rV, delta))
print("exact rank distribution over survivors:", dict(sorted(dist.items())))
for r in rank2: print("rank<=2:", r)
