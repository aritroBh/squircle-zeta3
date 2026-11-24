import time, sys
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/HEXC1_hexagon_certificate')
import numpy as np
import cmzv6
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

def rank_modp_fast(int_rows, ncols, p):
    M = np.array(int_rows, dtype=np.int64) % p
    nr, rank = M.shape[0], 0
    for j in range(ncols):
        if rank >= nr: break
        col = M[rank:, j]
        nz = np.nonzero(col)[0]
        if len(nz) == 0: continue
        piv = rank + int(nz[0])
        if piv != rank: M[[rank, piv]] = M[[piv, rank]]
        inv = pow(int(M[rank, j]), p-2, p)
        M[rank] = (M[rank]*inv) % p
        f = M[rank+1:, j].copy()
        M[rank+1:] = (M[rank+1:] - f[:, None]*M[rank]) % p
        rank += 1
    return rank

DELIGNE = {1: 3, 2: 8, 3: 21, 4: 55}
lower = {}
res = {}
ranks_seen = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    cols = list(cmzv6.words(n))
    ids = {w: j for j, w in enumerate(cols)}
    int_rows = [[0]*len(cols) for _ in range(len(rows))]
    for i, row in enumerate(rows):
        for w, c in row:
            int_rows[i][ids[w]] = c
    ranks = [rank_modp_fast(int_rows, len(cols), p) for p in (1000000007, 1000000009)]
    assert ranks[0] == ranks[1]
    ranks_seen[n] = ranks[0]
    res[n] = len(cols) - ranks[0]
    log(f"weight {n}: {len(rows)} rows, {len(cols)} words, rank {ranks[0]} -> quotient {res[n]} (Deligne {DELIGNE[n]}, missing {res[n]-DELIGNE[n]})")
assert res == {1: 3, 2: 8, 3: 23, 4: 69}
assert ranks_seen == {1: 3, 2: 28, 3: 229, 4: 1695}
log("quotient dims 3,8,23,69 = certificate row counts 28/229/1695 (all independent proved rows used); missing vs Deligne: 0,0,2,14")
with open('progress.log','a') as f:
    f.write("referee step(dims): quotient dims 3,8,23,69 (my GF(p), 2 primes); ranks 3/28/229/1695\n")
log("STEP PASS")
