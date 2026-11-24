import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
test = r'''"""Unit test: Echelon rank & membership vs flint nmod_mat on random sparse systems."""
import sys, time
sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import numpy as np, flint
from w5c2_lib import Echelon

p = 16777213
rng = np.random.default_rng(7)
for trial,(nr,nc,dens) in enumerate([(300,120,0.06),(500,200,0.03),(260,260,0.05)]):
    # random sparse matrix with exact rank via construction: A = L*U-ish -> just random
    A = (rng.random((nr,nc)) < dens) * rng.integers(1,p,(nr,nc))
    A = A % p
    M = flint.nmod_mat(nr, nc, p)
    for i in range(nr):
        for j in range(nc):
            if A[i,j]: M[i,j] = int(A[i,j])
    rk_flint = M.rank()
    rows = []
    for i in range(nr):
        cols = np.nonzero(A[i])[0].astype(np.int32)
        vals = A[i, cols].astype(np.int64)
        rows.append((cols, vals))
    E = Echelon(nc, p)
    S, C = E.feed(rows, t0=time.monotonic())
    assert E.rank() == rk_flint, (trial, E.rank(), rk_flint)
    # S rows independent & same span: check each original row reduces to 0 against echelon
    for i in [0, nr//2, nr-1]:
        rem = E.reduce(rows[i][0].tolist(), rows[i][1].tolist())
        assert not rem, (trial, i, rem[:3])
    # a random vector: membership agrees with flint solve
    okc = 0
    for t in range(20):
        b = rng.integers(0,p,nc)
        rem = E.reduce(np.nonzero(b)[0].astype(np.int32).tolist(), b[np.nonzero(b)].astype(np.int64).tolist())
        inspan_ech = not rem
        B = flint.nmod_mat(nr,1,p)
        # solve M^T x = b via rref of augmented [M^T | b]
        MA = flint.nmod_mat(nc, nr+1, p)
        for i in range(nr):
            for j in range(nc):
                if A[i,j]: MA[j,i] = int(A[i,j])
        for j in range(nc): MA[j,nr] = int(b[j])
        R, rk2 = MA.rref()
        # inconsistent iff some row has zeros in first nr cols and nonzero last
        incons = False
        for j in range(nc):
            rowzero = all(R[j,k] == 0 for k in range(nr))
            if rowzero and R[j,nr] != 0: incons = True; break
        inspan_flint = not incons
        assert inspan_ech == inspan_flint, (trial, t, inspan_ech, inspan_flint)
        okc += 1
    print(f'trial {trial}: nr={nr} nc={nc} rank={rk_flint} OK (S={len(S)}, membership {okc}/20)')
print('ECHELON UNIT TEST PASSED')
'''
open('step_01_unittest.py','w').write(test)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_01_unittest.py'], capture_output=True, text=True)
print(r.stdout, r.stderr[-2000:] if r.stderr else '')