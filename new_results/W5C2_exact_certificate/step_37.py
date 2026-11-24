import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
dixon = r'''"""W5C2 custom Dixon p-adic solver: solve A X = B over Q, A r×r integer (small entries),
B r×m integer. Returns (Xnum, Xden) with X = Xnum/Xden as Fractions per entry.

Method: Dixon lifting mod p with early termination via rational reconstruction,
stability check across extra steps, and modular cross-check at fresh primes.
Uses numpy int64 for exact A@X products (bounds checked) and flint nmod for mod-p matmul.
"""
import sys, os, time, math
import numpy as np
import flint

def log(s):
    print(time.strftime('%Y-%m-%d %H:%M:%S ')+s, flush=True)

def rat_recon(x, m, N):
    """a/b = x mod m with |a|<=N, 0<b<=N, or None. Standard extended-gcd reconstruction."""
    if x == 0: return (0, 1)
    r0, r1 = m, x % m
    t0, t1 = 0, 1
    while r1 > N:
        q = r0 // r1
        r0, r1 = r1, r0 - q*r1
        t0, t1 = t1, t0 - q*t1
    if r1 == 0: return None
    if t1 < 0: r1, t1 = -r1, -t1
    if t1 == 0 or t1 > N: return None
    if math.gcd(r1, t1) != 1: return None
    return (r1, t1)

def dixon_solve(A, B, p, logcb=log, recon_every=1, max_steps=100000):
    """A: (r,r) numpy int64 (exact). B: (r,m) numpy int64 (exact). p prime with det(A)!=0 mod p.
    Returns (Y, steps, M) with Y[r][m] = (num,den) Fraction pairs, or None if not converged."""
    r, m = B.shape
    assert A.shape == (r, r)
    maxA = int(np.abs(A).max()); maxB = int(np.abs(B).max())
    logcb(f'dixon: r={r} m={m} maxA={maxA} maxB_digits={len(str(maxB))} p={p}')
    # overflow bounds: A@Xk entries <= r*maxA*(p-1) ; Ainv@Bm <= r*(p-1)^2
    b1 = r*maxA*(p-1); b2 = r*(p-1)*(p-1)
    logcb(f'dixon: int64 bounds A@Xk<={b1:.3e} Ainv@Bm<={b2:.3e} (int64 max 9.22e18)')
    assert b1 < 4e18 and b2 < 4e18, 'int64 overflow risk'
    t0 = time.monotonic()
    Amod = flint.nmod_mat((A % p).astype(object).tolist(), p)
    logcb(f'dixon: nmod A built {time.monotonic()-t0:.1f}s; inverting')
    t0 = time.monotonic()
    Ainv = Amod.inv()
    logcb(f'dixon: nmod inv done {time.monotonic()-t0:.1f}s')
    Ainv_np = np.array(Ainv.tolist(), np.int64)
    A_nmod = flint.nmod_mat((A % p).astype(object).tolist(), p)  # reuse? not needed
    del A_nmod
    Xacc = np.zeros((r, m), object)
    Xacc[:] = 0
    Bcur = B.astype(np.int64).copy()
    M = 1
    pk = 1
    prevY = None
    stable = 0
    step = 0
    while step < max_steps:
        Bm = (Bcur % p).astype(np.int64)
        Xk = (Ainv_np @ Bm) % p                      # exact: entries < p^2*r
        AXk = (A @ Xk)                                # exact int64
        diff = Bcur - AXk
        assert (diff % p == 0).all()
        Bcur = diff // p
        Xacc = Xacc + pk * Xk.astype(object)
        M *= p; pk *= p; step += 1
        if step % recon_every == 0 or step >= 3:
            N = math.isqrt(M // 2)
            Y = np.empty((r, m), object)
            ok = True
            for i in range(r):
                for j in range(m):
                    rc = rat_recon(int(Xacc[i, j]) % M, M, N)
                    if rc is None: ok = False; break
                    Y[i, j] = rc
                if not ok: break
            if ok:
                if prevY is not None and all(Y[i, j] == prevY[i, j] for i in range(r) for j in range(m)):
                    stable += 1
                    logcb(f'dixon: step={step} reconstructed ALL and stable #{stable} M_digits={len(str(M))} el={time.monotonic()-t0:.0f}s')
                    if stable >= 1:
                        return Y, step, M
                else:
                    logcb(f'dixon: step={step} reconstructed ALL entries (first) M_digits={len(str(M))}')
                prevY = Y
            else:
                prevY = None
        if step % 10 == 0:
            logcb(f'dixon: step={step} M_digits={len(str(M))} el={time.monotonic()-t0:.0f}s')
    return None, step, M

def check_mod(Y, A, B, q):
    """Check A @ (num/den) == B mod q for all entries."""
    r, m = B.shape
    Am = A % q
    for j in range(m):
        lhs = np.zeros(r, np.int64)
        for i in range(r):
            n, d = Y[i, j]
            lhs[i] = (n % q) * pow(d % q, -1, q) % q
        prod = (Am @ lhs.reshape(-1,1)).reshape(-1) % q
        if not ((prod - B[:, j] % q) % q == 0).all():
            return False
    return True
'''
open('w5c2_dixon.py','w').write(dixon)
print('dixon module written')
import subprocess
print(subprocess.run(['tail','-3','step_01.out'], capture_output=True, text=True).stdout)