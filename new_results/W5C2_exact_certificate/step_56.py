import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
lib = open('w5c2_lib.py').read()
old = """            while True:
                nz = pc[buf[pc] != 0]
                if nz.size == 0: break
                c = int(nz[0])           # smallest creation-time pivot column present
                # find its pivot index: pivot cols are unique, creation order = index order
                t = int(np.searchsorted(pc, c))
                a = int(buf[c])
                buf = (buf - a*self.piv[t]) % p
                elims += 1
            self.elims += elims"""
new = """            lo = 0
            while lo < r:
                nz = pc[lo:][buf[pc[lo:]] != 0]
                if nz.size == 0: break
                c = int(nz[0])           # smallest creation-time pivot column present
                t = lo + int(np.searchsorted(pc[lo:], c))
                a = int(buf[c])
                buf = (buf - a*self.piv[t]) % p
                elims += 1
                lo = t + 1               # pivot rows > t have zeros in pivot cols <= t
            self.elims += elims"""
assert old in lib
open('w5c2_lib.py','w').write(lib.replace(old, new))
print('optimized')
# unit test DenseEchelon vs flint
test = r'''
import sys, time
sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import numpy as np, flint
from w5c2_lib import DenseEchelon
p = 16777213
rng = np.random.default_rng(11)
for trial,(nr,nc,dens) in enumerate([(300,120,0.06),(500,200,0.03),(260,260,0.05),(400,150,0.3)]):
    A = (rng.random((nr,nc)) < dens) * rng.integers(1,p,(nr,nc)) % p
    M = flint.nmod_mat(A.astype(object).tolist(), p)
    rk = M.rank()
    E = DenseEchelon(nc, p)
    rows = [(np.nonzero(A[i])[0].astype(np.int32), A[i,np.nonzero(A[i])[0]].astype(np.int64)) for i in range(nr)]
    E.feed(rows)
    assert E.r == rk, (trial, E.r, rk)
    for i in (0, nr//2, nr-1):
        buf = E.reduce(*rows[i])
        assert not buf.any(), (trial, i)
    # random vector membership vs flint
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
        R, rk2 = MA.rref()
        incons = any(all(R[j,k]==0 for k in range(nr)) and R[j,nr]!=0 for j in range(nc))
        assert ine == (not incons), (trial, t)
    print(f'trial {trial}: rank={rk} OK elims={E.elims}')
print('DENSE ECHELON UNIT TEST PASSED')
'''
open('step_02_unittest.py','w').write(test)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_02_unittest.py'], capture_output=True, text=True)
print(r.stdout, r.stderr[-1500:] if r.stderr else '')