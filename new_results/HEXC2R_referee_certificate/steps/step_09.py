"""step_09.py -- referee: quotient dimensions with MY OWN vectorized GF(p) elimination.
Quotient = #admissible words - rank(relation rows), weights 1-4, two primes each.
Claim: 3, 8, 23, 69 vs Deligne's (conjectured) 3, 8, 21, 55 -> missing 0,0,2,14 relations."""
import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import numpy as np
import cmzv6
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

def rank_modp_fast(int_rows, ncols, p):
    """my own GF(p) Gaussian elimination, row-vectorized (numpy int64)."""
    M = np.array(int_rows, dtype=np.int64) % p
    nr = M.shape[0]
    rank = 0
    for j in range(ncols):
        if rank >= nr: break
        col = M[rank:, j]
        nz = np.nonzero(col)[0]
        if len(nz) == 0: continue
        piv = rank + int(nz[0])
        if piv != rank:
            M[[rank, piv]] = M[[piv, rank]]
        inv = pow(int(M[rank, j]), p - 2, p)
        M[rank] = (M[rank] * inv) % p
        f = M[rank + 1:, j].copy()
        M[rank + 1:] = (M[rank + 1:] - f[:, None] * M[rank]) % p
        rank += 1
    return rank

DELIGNE = {1: 3, 2: 8, 3: 21, 4: 55}
lower = {}
res = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    cols = list(cmzv6.words(n))
    ids = {w: j for j, w in enumerate(cols)}
    int_rows = [[0] * len(cols) for _ in range(len(rows))]
    for i, row in enumerate(rows):
        for w, c in row:
            int_rows[i][ids[w]] = c
    ranks = []
    for p in (1000000007, 1000000009):
        t1 = time.monotonic()
        r = rank_modp_fast(int_rows, len(cols), p)
        ranks.append(r)
        log(f"weight {n}: {len(rows)} rows, {len(cols)} words, rank mod {p} = {r} "
            f"({time.monotonic()-t1:.1f}s)")
    assert ranks[0] == ranks[1]
    q = len(cols) - ranks[0]
    res[n] = q
    log(f"  -> quotient dim weight {n}: {q}  (Deligne {DELIGNE[n]}; missing {q - DELIGNE[n]})")
assert res == {1: 3, 2: 8, 3: 23, 4: 69}
log("quotients 3, 8, 23, 69 vs Deligne 3, 8, 21, 55  =>  0, 0, 2, 14 relations beyond the")
log("standard ones (frame double-shuffle type) — NOT used by the certificate.")
log("STEP 09 PASS")
with open("progress.log", "a") as f:
    f.write("step_09: quotient dims 3,8,23,69 (my GF(p) elimination, 2 primes) vs Deligne 3,8,21,55\n")
