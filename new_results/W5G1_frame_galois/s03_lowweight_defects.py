"""s03: Galois defects of the proved lower-weight bases (weights 1..4), mod p.
For each base row b of weight n and u in {7,3,5}: is u.b in P_n = span(base_n)?  If not, is it in
P_n + e sh Q^{words(n-1)}, e = (2)-(6)?  Also dim P_n vs dim of its largest Galois-stable subspace."""
from common import *
from flint import nmod_mat
p = 1000000007
B = bases()
for n in (1, 2, 3, 4):
    cols = words(n); ids = {w: i for i, w in enumerate(cols)}
    def vec(d):
        v = [0] * len(cols)
        for w, c in d.items(): v[ids[w]] = (c.numerator * pow(c.denominator, -1, p)) % p
        return v
    Pm = nmod_mat([vec(b) for b in B[n]], p)
    rP = Pm.rank()
    erows = [mul(E1, {w: 1}) for w in words(n - 1)] if n > 1 else [E1]
    PEm = nmod_mat([vec(b) for b in B[n]] + [vec({w: F(c) for w, c in r.items()}) for r in erows], p)
    rPE = PEm.rank()
    bad = {u: 0 for u in (7, 3, 5)}; badE = {u: 0 for u in (7, 3, 5)}
    for b in B[n]:
        for u in (7, 3, 5):
            gb = vec(actv(b, u))
            if nmod_mat([vec(x) for x in B[n]] + [gb], p).rank() > rP if n < 4 else None:
                pass
    # faster membership: stack all images at once per u and compare ranks incrementally via rref of P
    from flint import nmod_mat as M
    rr, rk = Pm.rref()
    piv = []; j = 0
    for i in range(rk):
        while rr[i, j] == 0: j += 1
        piv.append(j); j += 1
    def inspan(v, rr=rr, piv=piv):
        v = list(v)
        for i, pc in enumerate(piv):
            a = v[pc]
            if a:
                for jj in range(len(cols)):
                    x = rr[i, jj]
                    if x: v[jj] = (v[jj] - a * int(x)) % p
        return not any(v)
    rr2, rk2 = PEm.rref()
    piv2 = []; j = 0
    for i in range(rk2):
        while rr2[i, j] == 0: j += 1
        piv2.append(j); j += 1
    # vectorised membership via numpy
    import numpy as np
    R1 = np.array([[int(rr[i, jj]) for jj in range(len(cols))] for i in range(rk)], dtype=object) if n < 4 else None
    def member(v, rr, piv, rk):
        v = np.array(v, dtype=np.int64)
        for i, pc in enumerate(piv):
            a = int(v[pc])
            if a:
                row = ROWS[id(rr)][i]
                v = (v - a * row) % p
        return not v.any()
    ROWS = {id(rr): np.array([[int(rr[i, jj]) for jj in range(len(cols))] for i in range(rk)], dtype=np.int64),
            id(rr2): np.array([[int(rr2[i, jj]) for jj in range(len(cols))] for i in range(rk2)], dtype=np.int64)}
    for b in B[n]:
        for u in (7, 3, 5):
            gb = vec(actv(b, u))
            if not member(gb, rr, piv, rk):
                bad[u] += 1
                if not member(gb, rr2, piv2, rk2): badE[u] += 1
    log(f'weight {n}: cols={len(cols)} |base|={len(B[n])} rank P_n={rP} rank(P_n + e sh words(n-1))={rPE}; '
        f'u.b not in P_n: {bad}; not in P_n + e sh *: {badE}', 's03.log')
