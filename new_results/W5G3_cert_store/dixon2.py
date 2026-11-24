"""Dixon p-adic lifting for M Y = B over Q, M r x r sparse integer (small entries), B r x m integer.
Ainv = M^{-1} mod p (dense) is computed once with flint and cached on disk as uint32.
Lifting: X_k = Ainv (B_k mod p) mod p with float64 BLAS on 12-bit halves (exact: r*p*2^12 < 2^53),
B_{k+1} = (B_k - M X_k)/p.  Digits are combined two at a time (base p^2 fits int64) into Python-int accumulators.
Termination: per-column common denominator from probe entries, numerators by symmetric residue, then an
EXACT integer check M n_j == d_j B_j (Python ints).  Returns (n: object array r x m, d: list of m ints)."""
import time, math, os, sys
sys.set_int_max_str_digits(0)
import numpy as np, scipy.sparse as sps
import flint

def _log(s, logf):
    s = time.strftime('%Y-%m-%d %H:%M:%S ') + s
    print(s, flush=True)
    if logf:
        with open(logf, 'a') as fh: fh.write(s + '\n')

def modinv_dense(Mcsr, p, cache, logf=None):
    if cache and os.path.exists(cache):
        return np.load(cache).astype(np.float64)
    t0 = time.time()
    r = Mcsr.shape[0]
    D = (Mcsr.toarray() % p)
    flint.ctx.threads = 4
    Am = flint.nmod_mat(r, r, [int(x) for x in D.ravel()], p)
    del D
    _log(f'dixon: nmod_mat built ({time.time()-t0:.0f}s), inverting r={r}', logf)
    Ai = Am.inv(); del Am
    _log(f'dixon: inverse done ({time.time()-t0:.0f}s)', logf)
    ent = Ai.entries(); del Ai
    out = np.fromiter((int(x) for x in ent), dtype=np.uint32, count=r * r).reshape(r, r); del ent
    if cache: np.save(cache, out)
    _log(f'dixon: inverse extracted/saved ({time.time()-t0:.0f}s)', logf)
    return out.astype(np.float64)

def mulmod(Af, Bm, p):
    lo = (Bm & 4095).astype(np.float64); hi = (Bm >> 12).astype(np.float64)
    a = np.fmod(Af @ lo, p).astype(np.int64)
    b = np.fmod(Af @ hi, p).astype(np.int64)
    return (a + b * 4096) % p

def ratrec(x, m, N):
    if x == 0: return (0, 1)
    r0, r1, t0, t1 = m, x % m, 0, 1
    while r1 > N:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; t0, t1 = t1, t0 - q * t1
    if t1 == 0 or abs(t1) > N or math.gcd(r1, abs(t1)) != 1: return None
    if t1 < 0: r1, t1 = -r1, -t1
    return (r1, t1)

def exact_check(Mcsr, n, d, B):
    """M n_j == d_j B_j exactly (Python ints)."""
    indptr, indices, data = Mcsr.indptr, Mcsr.indices, Mcsr.data
    r, m = n.shape
    for j in range(m):
        col = n[:, j]; dj = d[j]
        for i in range(r):
            s = 0
            for k in range(indptr[i], indptr[i + 1]):
                s += int(data[k]) * col[indices[k]]
            if s != dj * int(B[i, j]): return False
    return True

def dixon_solve(Mcsr, B, p, Af, logf=None, check_every=16, nprobe=64, maxsteps=100000):
    r, m = B.shape
    Mcsr = Mcsr.tocsr().astype(np.int64)
    maxA = int(abs(Mcsr).max())
    assert r * maxA * p < 2 ** 62
    t0 = time.time()
    big = max(abs(int(x)) for x in np.asarray(B, dtype=object).ravel()) >= 2 ** 62
    Bcur = np.array(B, dtype=object) if big else np.array(B, dtype=np.int64)
    acc = np.zeros((r, m), dtype=object); acc[:] = 0
    M = 1; step = 0; pend = None
    rng = np.random.default_rng(7)
    probe = rng.integers(0, r, nprobe)
    prev = None
    while step < maxsteps:
        Xk = mulmod(Af, (Bcur % p).astype(np.int64), p)
        AX = Mcsr @ Xk
        diff = Bcur - AX
        assert not (diff % p).any()
        Bcur = diff // p
        if Bcur.dtype == object and max(abs(int(x)) for x in Bcur.ravel()) < 2 ** 61:
            Bcur = Bcur.astype(np.int64)
        if pend is None:
            pend = Xk
        else:
            blk = pend + Xk * p                # < p^2 < 2^48
            acc += blk.astype(object) * M
            M *= p * p; pend = None
        step += 1
        if pend is None and step % check_every == 0:
            Nb = math.isqrt(M // 2)
            dens = []
            ok = True
            for j in range(m):
                dj = 1
                for i in probe:
                    x = int(acc[i, j]) * dj % M
                    rc = ratrec(x, M, Nb)
                    if rc is None: ok = False; break
                    dj *= rc[1]
                if not ok: break
                dens.append(dj)
            if ok:
                key = tuple(dens)
                if prev == key:
                    # full candidate
                    n = np.empty((r, m), dtype=object)
                    good = True
                    half = M // 2
                    for j in range(m):
                        dj = dens[j]
                        colj = acc[:, j]
                        i = 0
                        while i < r:
                            x = int(colj[i]) * dj % M
                            if x > half: x -= M
                            if abs(x) > Nb:
                                rc = ratrec(x % M, M, Nb)
                                if rc is None or rc[1] == 1: good = False; break
                                dj *= rc[1]; i = 0; continue      # enlarge the common denominator, restart column
                            n[i, j] = x; i += 1
                        if not good: break
                        dens[j] = dj
                    if good:
                        # quick check mod a 61-bit prime before the exact check
                        q = 2147483647
                        nm = np.array([[int(n[i, j]) % q for j in range(m)] for i in range(r)], dtype=np.int64)
                        lhs = (Mcsr @ nm) % q
                        rhs = np.array([[dens[j] % q * (int(B[i, j]) % q) % q for j in range(m)] for i in range(r)], dtype=np.int64)
                        if not ((lhs - rhs) % q).any():
                            _log(f'dixon: candidate at step {step} (M {len(str(M))} digits) passes mod 2^31-1; exact check ...', logf)
                            if exact_check(Mcsr, n, dens, B):
                                _log(f'dixon: EXACT solution, step {step}, {time.time()-t0:.0f}s, max den digits {max(len(str(x)) for x in dens)}', logf)
                                return n, dens, step
                            _log('dixon: exact check failed; continuing', logf)
                prev = key
            if step % 64 == 0:
                _log(f'dixon: step {step}, M digits {len(str(M))}, probe ok {ok}, {time.time()-t0:.0f}s', logf)
    raise RuntimeError('dixon: no convergence')
