"""Own dense elimination mod p (p < 2^22) with float64 BLAS, blocked (panel width b <= 128 so k*p^2 < 2^53).
elim(A, Tm, p): A (m x n) generator rows, Tm (k x n) test rows (never used as pivots).
Returns rank of A and the residues of the test rows modulo rowspace(A) (zero on pivot columns).
A test row lies in rowspace(A) mod p iff its residue is zero. A and Tm are overwritten.
Per panel: pivots found by sequential elimination inside the panel; the trailing columns are updated with the
Schur complement  A_o <- A_o - A_{o,J} A_{R,J}^{-1} A_R  (pre-panel values), which equals sequential elimination."""
import numpy as np, time

def inv_small(M, p):
    k = M.shape[0]
    A = np.concatenate([np.mod(M, p), np.eye(k)], axis=1)
    for c in range(k):
        nz = np.flatnonzero(A[c:, c] != 0)
        assert len(nz), 'singular pivot block'
        i = c + nz[0]
        if i != c: A[[c, i]] = A[[i, c]]
        A[c] = np.mod(A[c] * pow(int(A[c, c]), p - 2, p), p)
        f = A[:, c].copy(); f[c] = 0
        A -= np.outer(f, A[c]); np.mod(A, p, out=A)
    return A[:, k:]

def elim(A, Tm, p, b=128, chunk=4096, logf=None):
    assert p < 2 ** 22 and b <= 128
    m, n = A.shape
    r = 0; t0 = time.time()
    for j0 in range(0, n, b):
        j1 = min(n, j0 + b)
        P0 = A[r:, j0:j1].copy(); P = P0.copy()
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
            X = np.mod(Minv @ A[r + Rl, j1:], p)                     # k x n_trail
            Tm[:, j1:] = np.mod(Tm[:, j1:] - np.mod(T0[:, Jl] @ X, p), p)
            L = np.ascontiguousarray(P0[:, Jl])
            for s in range(r, m, chunk):
                e = min(m, s + chunk)
                blk = A[s:e, j1:]
                blk -= np.mod(L[s - r:e - r] @ X, p)
                np.mod(blk, p, out=blk)
            # pivot rows are now zero on trailing columns as well (A_{R,J} Minv = I)
        # bring the k pivot rows (now all-zero) to positions r..r+k-1 so the active block stays contiguous
        tg = np.arange(r, r + k); src = r + Rl
        S = np.union1d(tg, src)
        rest_src = np.setdiff1d(S, src)          # non-pivot rows inside S, in order
        rest_pos = np.setdiff1d(S, tg)
        newA = A[np.concatenate([src, rest_src])]
        A[np.concatenate([tg, rest_pos])] = newA
        r += k
        if logf and (j0 // b) % 20 == 0:
            logf(f'    elim: cols {j1}/{n}, rank {r}, {time.time()-t0:.0f}s')
    return r, Tm
