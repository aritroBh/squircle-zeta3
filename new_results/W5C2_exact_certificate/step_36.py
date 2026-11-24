import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
step03test = r'''"""W5C2 step_03_test: timing/memory probes for the exact solve at weight-5 sizes.
1) nmod_mat.inv at 12100 (for custom Dixon); 2) fmpz_mat.solve at increasing sizes
(entries +-300, 2 RHS) to detect Dixon vs FFLU scaling. Logs RSS via /usr/bin/time not
available mid-run; we print timestamps; external watchdog may kill."""
import sys, os, time
sys.dont_write_bytecode=True
import numpy as np, flint
def log(s):
    print(time.strftime('%H:%M:%S ')+s, flush=True)
rng = np.random.default_rng(1)
mode = sys.argv[1] if len(sys.argv)>1 else 'all'
if mode in ('inv','all'):
    n = 12100; p = 16777213
    A = rng.integers(-300, 300, (n,n)) % p
    M = flint.nmod_mat(n, n, p)
    # fill via tolist for speed
    M = flint.nmod_mat(A.astype(object).tolist(), p)
    log(f'nmod 12100 built; starting inv')
    t0=time.monotonic()
    Mi = M.inv()
    log(f'nmod_mat.inv(12100) took {time.monotonic()-t0:.1f}s')
    t0=time.monotonic()
    P = M*Mi
    ok = all(P[i,i]==1 for i in range(n))
    log(f'inv check diag-ones={ok} ({time.monotonic()-t0:.1f}s)')
if mode in ('solve','all'):
    for n in (2000, 4000, 6000):
        A = rng.integers(-300, 300, (n,n))
        B = rng.integers(-10**6, 10**6, (n,2))
        Af = flint.fmpz_mat(A.astype(object).tolist())
        Bf = flint.fmpz_mat(B.astype(object).tolist())
        log(f'fmpz_mat.solve at n={n}: starting')
        t0=time.monotonic()
        X = Af.solve(Bf)
        dt = time.monotonic()-t0
        h = max(max(abs(X[i,j].numerator.bit_length()), X[i,j].denominator.bit_length()) for i in range(n) for j in range(2))
        log(f'fmpz_mat.solve n={n} took {dt:.1f}s maxheight_bits={h}')
log('TEST DONE')
'''
open('step_03_test.py','w').write(step03test)
print('step_03_test.py written')