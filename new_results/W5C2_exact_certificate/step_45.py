import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
dixon = r'''"""W5C2 custom Dixon p-adic solver: solve A X = B over Q, A r×r integer (small entries),
B r×m integer. Dixon lifting mod p with early termination via rational reconstruction,
stability requirement, and checkpoint/resume. Returns Y[i][j]=(num,den) with A Y = B."""
import time, math, pickle, os
import numpy as np
import flint

def log(s):
    print(time.strftime('%Y-%m-%d %H:%M:%S ')+s, flush=True)

def rat_recon(x, m, N):
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

def dixon_solve(A, B, p, ckpt_path=None, recon_probe=256):
    r, m = B.shape
    maxA = int(np.abs(A).max())
    b1 = r*maxA*(p-1); b2 = r*(p-1)*(p-1)
    log(f'dixon: r={r} m={m} maxA={maxA} int64 bounds {b1:.2e} {b2:.2e} (<9.2e18)')
    assert b1 < 4e18 and b2 < 4e18
    t0 = time.monotonic()
    Amod = flint.nmod_mat((A % p).astype(object).tolist(), p)
    log(f'dixon: nmod A built {time.monotonic()-t0:.1f}s; inverting')
    t0 = time.monotonic()
    Ainv = np.array(Amod.inv().tolist(), np.int64)
    log(f'dixon: nmod inv done {time.monotonic()-t0:.1f}s')
    # state
    step = 0; M = 1; pk = 1
    Xacc = np.zeros((r, m), object); Xacc[:] = 0
    Bcur = B.astype(np.int64).copy()
    prevY = None
    if ckpt_path and os.path.exists(ckpt_path):
        ck = pickle.load(open(ckpt_path,'rb'))
        step, M, pk = ck['step'], ck['M'], ck['pk']
        Xacc, Bcur, prevY = ck['Xacc'], ck['Bcur'], ck['prevY']
        log(f'dixon: resumed at step={step} M_digits={len(str(M))}')
    probe_idx = np.random.default_rng(3).integers(0, r, min(recon_probe, r))
    while True:
        Bm = (Bcur % p).astype(np.int64)
        Xk = (Ainv @ Bm) % p
        AXk = A @ Xk
        diff = Bcur - AXk
        assert (diff % p == 0).all()
        Bcur = diff // p
        Xacc = Xacc + pk * Xk.astype(object)
        M *= p; pk *= p; step += 1
        # probe-based reconstruction filter
        N = math.isqrt(M // 2)
        ok = True
        for i in probe_idx:
            for j in range(m):
                if rat_recon(int(Xacc[i, j]) % M, M, N) is None: ok = False; break
            if not ok: break
        if ok:
            Y = np.empty((r, m), object)
            for i in range(r):
                for j in range(m):
                    rc = rat_recon(int(Xacc[i, j]) % M, M, N)
                    if rc is None: ok = False; break
                    Y[i, j] = rc
                if not ok: break
            if ok:
                if prevY is not None and all(Y[i, j] == prevY[i, j] for i in range(r) for j in range(m)):
                    log(f'dixon: CONVERGED step={step} stable-repeat, M_digits={len(str(M))} el={time.monotonic()-t0:.0f}s')
                    if ckpt_path and os.path.exists(ckpt_path): os.remove(ckpt_path)
                    return Y, step, M
                log(f'dixon: step={step} full reconstruction succeeded (need repeat) M_digits={len(str(M))}')
                prevY = Y
        if step % 25 == 0:
            log(f'dixon: step={step} M_digits={len(str(M))} el={time.monotonic()-t0:.0f}s')
            if ckpt_path:
                pickle.dump({'step':step,'M':M,'pk':pk,'Xacc':Xacc,'Bcur':Bcur,'prevY':prevY},
                            open(ckpt_path,'wb'), protocol=4)

def check_mod(Y, A, B, q):
    """A @ Y == B mod q, Y[i][j]=(num,den)."""
    r, m = B.shape
    Am = A % q
    for j in range(m):
        yj = np.array([(Y[i][j][0] % q) * pow(Y[i][j][1] % q, -1, q) % q for i in range(r)], np.int64)
        prod = (Am @ yj) % q
        if not ((prod - B[:, j] % q) % q == 0).all(): return False
    return True
'''
open('w5c2_dixon.py','w').write(dixon)
# fix step_03 cross-check call
s = open('step_03.py').read()
s = s.replace("ok = dx.check_mod([[ (F(Y[i][j][0],Y[i][j][1])/Ls[('rat','sq')[j]]) for j in range(2)] for i in range(r)], A, Bs, P2)",
              "ok = dx.check_mod(Y, A, Bs, P2)")
s = s.replace("# note: check_mod checks A@y == B/L mod p2 i.e. A@(Y/L)... construct Bscaled = B (since y=Y/L, A Y = B)\n", "")
open('step_03.py','w').write(s)
print('dixon rewritten, step_03 fixed')
import subprocess
print(subprocess.run(['tail','-2','step_01.out'], capture_output=True, text=True).stdout)