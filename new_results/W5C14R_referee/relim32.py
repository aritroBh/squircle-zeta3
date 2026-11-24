"""relim.elim (W5G2R's own blocked elimination mod p) with the generator matrix STORED as float32 to halve memory.
Residues are < p < 2^22 < 2^24, so float32 holds them exactly; every arithmetic step is done in float64 exactly as in
relim.elim (panel width <= 128, so k*p^2 < 2^53). Only the storage type differs."""
import numpy as np, time
from relim import inv_small

def elim32(A, Tm, p, b=128, chunk=2048, logf=None):
    assert p < 2 ** 22 and b <= 128 and A.dtype == np.float32 and Tm.dtype == np.float64
    m, n = A.shape
    r = 0; t0 = time.time()
    for j0 in range(0, n, b):
        j1 = min(n, j0 + b)
        P0 = A[r:, j0:j1].astype(np.float64); P = P0.copy()
        T0 = Tm[:, j0:j1].copy(); PT = T0.copy()
        used = np.zeros(m - r, bool); piv = []
        for c in range(j1 - j0):
            col = P[:, c]
            cand = np.flatnonzero((col != 0) & ~used)
            if not len(cand): continue
            i = cand[0]
            P[i, c:] = np.mod(P[i, c:] * pow(int(col[i]), p - 2, p), p)
            used[i] = True; piv.append((i, c))
            f = P[:, c] * (~used)
            nz = np.flatnonzero(f)
            if len(nz) > 0.3 * len(f):
                P[:, c:] = np.mod(P[:, c:] - np.outer(f, P[i, c:]), p)
            elif len(nz):
                P[nz, c:] = np.mod(P[nz, c:] - f[nz, None] * P[i, c:][None, :], p)
            ft = PT[:, c]; nzt = np.flatnonzero(ft)
            if len(nzt):
                PT[nzt, c:] = np.mod(PT[nzt, c:] - ft[nzt, None] * P[i, c:][None, :], p)
        Tm[:, j0:j1] = PT
        k = len(piv)
        A[r:, j0:j1] = 0
        if k == 0: continue
        Rl = np.array([i for i, _ in piv]); Jl = np.array([c for _, c in piv])
        Minv = inv_small(P0[np.ix_(Rl, Jl)], p)
        if j1 < n:
            X = np.mod(Minv @ A[r + Rl, j1:].astype(np.float64), p)
            Tm[:, j1:] = np.mod(Tm[:, j1:] - np.mod(T0[:, Jl] @ X, p), p)
            L = np.ascontiguousarray(P0[:, Jl])
            for s in range(r, m, chunk):
                e = min(m, s + chunk)
                blk = A[s:e, j1:].astype(np.float64)
                blk -= np.mod(L[s - r:e - r] @ X, p)
                np.mod(blk, p, out=blk)
                A[s:e, j1:] = blk
        tg = np.arange(r, r + k); src = r + Rl
        S = np.union1d(tg, src)
        rest_src = np.setdiff1d(S, src)
        rest_pos = np.setdiff1d(S, tg)
        newA = A[np.concatenate([src, rest_src])]
        A[np.concatenate([tg, rest_pos])] = newA
        r += k
        if logf and (j0 // b) % 20 == 0:
            logf(f'    elim: cols {j1}/{n}, rank {r}, {time.time()-t0:.0f}s')
    return r, Tm
