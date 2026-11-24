
import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import numpy as np
import cmzv6

def rank_modp_fast(int_rows, ncols, p):
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
        f = M[rank+1:, j].copy()
        M[rank+1:] = (M[rank+1:] - f[:, None] * M[rank]) % p
        rank += 1
    return rank

lower = {}
t0 = time.monotonic()
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    cols = list(cmzv6.words(n))
    ids = {w: j for j, w in enumerate(cols)}
    int_rows = [[0]*len(cols) for _ in range(len(rows))]
    for i, row in enumerate(rows):
        for w, c in row:
            int_rows[i][ids[w]] = c
    for p in (1000000007, 1000000009):
        t1 = time.monotonic()
        r = rank_modp_fast(int_rows, len(cols), p)
        print(f"weight {n}: {len(cols)} words, rank mod {p} = {r}  -> quotient {len(cols)-r}  ({time.monotonic()-t1:.1f}s)")
print(f"total {time.monotonic()-t0:.1f}s")
