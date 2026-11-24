import time, pickle
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
from math import gcd, isqrt

BND, H = 200, 40
# smallest-prime-factor sieve up to H*H*200+... norms |a^2 - D b^2| <= 200*1600 = 320000
LIM = 320000
spf = list(range(LIM+1))
for i in range(2, isqrt(LIM)+1):
    if spf[i] == i:
        for j in range(i*i, LIM+1, i):
            if spf[j] == j: spf[j] = i
def vpairs(n):
    """factorization of |n| as dict."""
    f = {}
    n = abs(n)
    while n > 1:
        p = spf[n]; c = 0
        while n % p == 0: n //= p; c += 1
        f[p] = f.get(p, 0) + c
    return f
FACE = {e: vpairs(e) for e in range(1, H+1)}

def rank_int(M):
    """rank of integer matrix (exact, via Fraction-free Gaussian elimination)."""
    M = [row[:] for row in M if any(row)]
    if not M: return 0
    R = 0
    for c in range(len(M[0])):
        p = next((r for r in range(R, len(M)) if M[r][c]), None)
        if p is None: continue
        M[R], M[p] = M[p], M[R]
        piv = M[R][c]
        for r in range(len(M)):
            if r != R and M[r][c]:
                f = M[r][c]
                M[r] = [a*piv - b*f for a, b in zip(M[r], M[R])]
        R += 1
        if R == len(M): break
    return R

sqfree = [D for D in range(2, BND+1) if all(e == 1 for e in vpairs(D).values())]
t0 = time.time()
survivors = {}
tot = 0
for D in sqfree:
    sq = D**0.5
    surv = []
    for e in range(1, H+1):
        fe = FACE[e]
        for b in range(-H, H+1):
            if b == 0: continue
            Db2 = D*b*b
            for a in range(-H, H+1):
                g = gcd(gcd(abs(a), abs(b)), e)
                if g != 1: continue
                if not (0 < a + b*sq < e): continue
                tot += 1
                n0 = a*a - Db2
                n1 = (e-a)*(e-a) - Db2
                n2 = (e+a)*(e+a) - Db2
                f0 = vpairs(n0); f1 = vpairs(n1); f2 = vpairs(n2)
                primes = sorted(set(f0) | set(f1) | set(f2) | set(fe) | {2})
                rows = []
                for f in ({2: 2}, f0, f1, f2):
                    rows.append([f.get(p, 0) - 2*fe.get(p, 0) if f is not ({2:2}) else f.get(p, 0) for p in primes])
                # note: for letters c, 1-c, 1+c subtract 2 v_p(e)
                rows[1] = [f0.get(p, 0) - 2*fe.get(p, 0) for p in primes]
                rows[2] = [f1.get(p, 0) - 2*fe.get(p, 0) for p in primes]
                rows[3] = [f2.get(p, 0) - 2*fe.get(p, 0) for p in primes]
                if rank_int(rows) <= 2:
                    surv.append((a, b, e))
    if surv:
        survivors[D] = surv
    print(f"D={D}: {len(surv)} norm-rank<=2 survivors   ({time.time()-t0:.0f}s)", flush=True)
print("total candidates:", tot)
pickle.dump(survivors, open(WD+"census_filter_survivors.pkl", "wb"))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_56: integer norm-rank filter, D<=200, heights<=40.\n")
