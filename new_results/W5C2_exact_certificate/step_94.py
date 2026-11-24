import os, subprocess
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
print(subprocess.run(['bash','-c','grep "step_02 pp" step_02.out'], capture_output=True, text=True).stdout)
st = r'''"""step_03_test: (1) nmod_mat.inv timing at 12156 (custom Dixon cost); (2) fmpz_mat.solve scaling at 4000/6000."""
import sys, time
import numpy as np, flint
def log(s): print(time.strftime('%H:%M:%S ')+s, flush=True)
rng = np.random.default_rng(1)
n = 12156; p = 16777213
A = rng.integers(-300, 300, (n,n)) % p
M = flint.nmod_mat(A.astype(object).tolist(), p)
log('nmod 12156 built; inv start')
t0=time.monotonic(); Mi = M.inv(); log(f'nmod inv(12156) {time.monotonic()-t0:.1f}s')
del M, Mi
for n2 in (4000, 6000):
    A2 = rng.integers(-300,300,(n2,n2)); B2 = rng.integers(-10**6,10**6,(n2,2))
    Af = flint.fmpz_mat(A2.astype(object).tolist()); Bf = flint.fmpz_mat(B2.astype(object).tolist())
    t0=time.monotonic(); X = Af.solve(Bf); dt=time.monotonic()-t0
    h = max(max(abs(X[i,j].numerator.bit_length()), X[i,j].denominator.bit_length()) for i in range(0,n2,97) for j in range(2))
    log(f'fmpz solve n={n2} {dt:.1f}s sample_height_bits={h}')
log('TEST DONE')
'''
open('step_03_test.py','w').write(st)
print('test written')