import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
lib = open('w5c2_lib.py').read()
fast = '''

class FastEchelon:
    """Sparse row echelon over GF(p), batched-round elimination via np.bincount.
    Pivot rows stored sparse (idx, vals), normalized to lead 1.
    Invariant: pivot row created at time t has zeros in all pivot cols of time < t.
    Rounds eliminate ALL currently-hit pivot cols at once; min creation time of hits
    strictly increases each round, so termination <= r rounds (typically few)."""
    def __init__(self, ncols, p):
        self.NC = ncols; self.p = p
        self.piv_idx = []; self.piv_val = []
        self.pivcol = np.full(ncols, -1, np.int64)
        self.pivot_of = np.full(ncols, -1, np.int64)
        self.colcnt = np.zeros(ncols, np.int64)
        self.pivrow = []
        self.r = 0
        self.rounds_tot = 0
        self.ops = 0
    def reduce(self, cols, vals):
        p = self.p; NC = self.NC
        buf = np.zeros(NC, np.int64)
        buf[cols] = np.asarray(vals, np.int64) % p
        r = self.r
        if r:
            pc = self.pivcol[:r]
            for _rnd in range(1000000):
                hit = pc[buf[pc] != 0]
                if hit.size == 0: break
                coeffs = buf[hit]           # < p already (per-round mod)
                tidx = self.pivot_of[hit]
                parts_i = [self.piv_idx[t] for t in tidx]
                parts_v = [self.piv_val[t] for t in tidx]
                cat_i = np.concatenate(parts_i)
                cat_v = np.concatenate(parts_v)
                crep = np.repeat(coeffs, [len(x) for x in parts_i])
                upd = np.bincount(cat_i, weights=crep*cat_v, minlength=NC)
                buf = (buf - upd) % p
                self.rounds_tot += 1
                self.ops += len(cat_i)
        return buf
    def add_pivot(self, buf, src=-1):
        nz = np.nonzero(buf)[0]
        if nz.size == 0: return -1
        # Markowitz: min colcnt column
        c = int(min(nz.tolist(), key=lambda cc: (self.colcnt[cc], cc)))
        lead = int(buf[c]) % self.p
        inv = pow(lead, -1, self.p)
        self.piv_idx.append(nz.astype(np.int32))
        self.piv_val.append((buf[nz]*inv) % self.p)
        self.pivcol[self.r] = c
        self.pivot_of[c] = self.r
        self.colcnt[nz] += 1
        self.pivrow.append(src)
        self.r += 1
        return c
    def feed(self, rows, label='', t0=None, logevery=4000, stopcheck=None):
        import time as _time
        for ri in range(len(rows)):
            cols, vals = rows[ri]
            buf = self.reduce(cols, vals)
            if buf.any(): self.add_pivot(buf, src=ri)
            if (ri+1) % logevery == 0:
                el = _time.monotonic()-t0 if t0 else 0
                nnz = [len(x) for x in self.piv_idx[-2000:]] or [0]
                log(f'{label} row {ri+1}/{len(rows)} pivots={self.r} rounds/row={self.rounds_tot/(ri+1):.1f} nnz_avg={sum(nnz)/len(nnz):.0f} max={max(nnz)} ops={self.ops/1e9:.2f}e9 el={el:.0f}s')
            if stopcheck and (ri+1) % stopcheck[0] == 0:
                if stopcheck[1](self):
                    log(f'{label} early stop at row {ri+1} pivots={self.r}')
                    return self.pivrow, list(self.pivcol[:self.r]), ri+1
        return self.pivrow, list(self.pivcol[:self.r]), None
'''
open('w5c2_lib.py','w').write(lib + fast)
print('FastEchelon added')
# quick unit test vs flint
test = r'''
import sys, time
sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import numpy as np, flint
from w5c2_lib import FastEchelon
p = 16777213
rng = np.random.default_rng(13)
for trial,(nr,nc,dens) in enumerate([(300,120,0.06),(500,200,0.03),(400,150,0.3),(2000,500,0.02)]):
    A = (rng.random((nr,nc)) < dens) * rng.integers(1,p,(nr,nc)) % p
    M = flint.nmod_mat(A.astype(object).tolist(), p)
    rk = M.rank()
    E = FastEchelon(nc, p)
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
    print(f'trial {trial}: rank={rk} OK rounds/row={E.rounds_tot/nr:.1f}')
print('FAST ECHELON UNIT TEST PASSED')
'''
open('step_03_unittest.py','w').write(test)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_03_unittest.py'], capture_output=True, text=True)
print(r.stdout, r.stderr[-1500:] if r.stderr else '')