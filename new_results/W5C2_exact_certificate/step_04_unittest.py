
import sys, time
sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import numpy as np, flint
from w5c2_lib import RrefEchelon
p = 16777213
rng = np.random.default_rng(13)
for trial,(nr,nc,dens) in enumerate([(300,120,0.06),(500,200,0.03),(400,150,0.3),(2000,500,0.02)]):
    A = (rng.random((nr,nc)) < dens) * rng.integers(1,p,(nr,nc)) % p
    M = flint.nmod_mat(A.astype(object).tolist(), p)
    rk = M.rank()
    E = RrefEchelon(nc, p)
    rows = [(np.nonzero(A[i])[0].astype(np.int32), A[i,np.nonzero(A[i])[0]].astype(np.int64)) for i in range(nr)]
    E.feed(rows)
    assert E.r == rk, (trial, E.r, rk)
    for i in (0, nr//2, nr-1):
        buf = E.reduce(*rows[i])
        assert not buf.any(), (trial, i)
    for t in range(10):
        b = rng.integers(0,p,nc)
        buf = E.reduce(np.nonzero(b)[0].astype(np.int32), b[np.nonzero(b)].astype(np.int64))
        ine = not buf.any()
        MA = flint.nmod_mat(nc, nr+1, p)
        At = A.T
        for j in range(nc):
            for i in range(nr):
                if At[j,i]: MA[j,i] = int(At[j,i])
            MA[j,nr] = int(b[j])
        R, _ = MA.rref()
        incons = any(all(R[j,k]==0 for k in range(nr)) and R[j,nr]!=0 for j in range(nc))
        assert ine == (not incons), (trial, t)
    print(f'trial {trial}: rank={rk} OK ops={E.ops}')
print('FAST ECHELON UNIT TEST PASSED')
