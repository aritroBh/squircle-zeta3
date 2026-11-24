import json, sympy as sp
from cyclo import Cyclo, mobius_images
from itertools import permutations

def prime_for(N, kmin=50):
    """smallest prime p = kN+1 with k >= kmin (so the mu_N subgroup is <= 1/50 of F_p^*)."""
    k = kmin
    while True:
        p = k*N + 1
        if sp.isprime(p): return int(p)
        k += 1

def primroot_power(p, N):
    g = sp.primitive_root(p)
    return pow(int(g), (p-1)//N, p)

results = {}
for N in (6, 12, 18, 24, 30, 36):
    p = prime_for(N)
    z = primroot_power(p, N)
    roots = [pow(z, k, p) for k in range(N)]
    rset = set(roots)
    strict_c = []; relaxed_extra = []
    for i, j, k in permutations(range(N), 3):
        a, b, c = roots[i], roots[j], roots[k]
        bma = (b-a) % p; cmb = (c-b) % p
        t1 = c*bma % p; t2 = a*cmb % p
        cls = []
        for num, den in (((t1 + 3*t2) % p, (bma + 3*cmb) % p), ((t2 - t1) % p, (cmb - bma) % p)):
            if den == 0: cls.append('inf')
            elif num == 0: cls.append('zero')
            else:
                u = num*pow(den, -1, p) % p
                cls.append('root' if u in rset else 'other')
        if 'other' not in cls:
            if cls == ['root','root']:
                strict_c.append((i,j,k))
            else:
                relaxed_extra.append((i,j,k, cls))
    # exact verification (mod-p filter is a complete superset: no false negatives)
    C = Cyclo(N)
    strict_ok = []
    for i,j,k in strict_c:
        im = list(mobius_images(C, i, j, k))
        if all(x[0]=='root' for x in im):
            exps = (i, im[0][1], j, im[1][1], k)
            if len(set(exps)) == 5:
                strict_ok.append(exps)
    relaxed_ok = []
    for i,j,k,cls in relaxed_extra:
        im = list(mobius_images(C, i, j, k))
        if all(x[0] in ('root','zero','inf') for x in im):
            relaxed_ok.append((i,j,k, im))
    results[N] = dict(p=p, k=(p-1)//N, modp_candidates=len(strict_c)+len(relaxed_extra),
                      strict=len(strict_ok), relaxed_extra=len(relaxed_ok),
                      exps=[list(e) for e in strict_ok])
    print(f"N={N}: p={p} (k={results[N]['k']}) modp_cand={results[N]['modp_candidates']} "
          f"EXACT strict={len(strict_ok)} extra_relaxed={len(relaxed_ok)}")
    # structure: all maps = shifts of two base maps?
    S = set(strict_ok)
    bases = set()
    for e in strict_ok:
        bases.add(tuple((x - e[0]) % N for x in e))
    print(f"   normalized (subtract image of 0): {sorted(bases)}")
json.dump(results, open('map_counts.json','w'), indent=1)
with open('progress.log','a') as f:
    f.write('step_22 (redo, sharp primes): exact strict map counts ' +
            str({n: r['strict'] for n,r in results.items()}) +
            '; extra relaxed (0/inf) ' + str({n: r['relaxed_extra'] for n,r in results.items()}) + '\n')
