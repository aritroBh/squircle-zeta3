import os, subprocess
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
# sanity: venv python has flint
r = subprocess.run(['<project>/.venv/bin/python','-c','import flint, numpy, scipy; print(flint.__version__)'], capture_output=True, text=True)
print(r.stdout, r.stderr)

lib = r'''"""W5C2 shared library: sector data loading, mod-p sparse echelon with pivot-row tracking.
Adapted from W5C1 w5_elim3_lib.py (heap-ordered pivot elimination = correct GE)."""
import sys, pickle, time, os, json, heapq
import numpy as np
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
W5C1 = os.path.join(HERE, '..', 'W5C1_certificate_scoping')
SECTORS = [(1,1),(1,-1),(-1,1),(-1,-1)]
SECNAME = {(1,1):'pp',(1,-1):'pm',(-1,1):'mp',(-1,-1):'mm'}
P1 = 16777213
P2 = 16777199

def log(s):
    s = time.strftime('%Y-%m-%d %H:%M:%S ') + s
    print(s, flush=True)
    with open(os.path.join(HERE,'progress.log'),'a') as f: f.write(s+'\n')

def load_sector_rows(sec):
    """Combined row list for sector sec=(e0,e1): [(fam, famidx, ((col,val),...) sorted)], stable-sorted by length."""
    a,b = sec
    eng = pickle.load(open(os.path.join(W5C1,f'sector/rows_engine_{a}_{b}.pkl'),'rb'))
    f234 = pickle.load(open(os.path.join(W5C1,f'sector/rows_frame234_{a}_{b}.pkl'),'rb'))
    fw5 = pickle.load(open(os.path.join(W5C1,f'sector/rows_framew5_{a}_{b}.pkl'),'rb'))
    rows = [('engine', i, r) for i,r in enumerate(eng)]
    rows += [('frame234', i, r) for i,r in enumerate(f234)]
    rows += [('framew5', i, tuple(sorted(r.items()))) for i,r in enumerate(fw5)]
    rows.sort(key=lambda t: (len(t[2]), t[0]))  # length first; family name tiebreak for full determinism
    return rows

def load_D(sec, comp):
    a,b = sec
    return pickle.load(open(os.path.join(W5C1,f'sector/D_{comp}_{a}_{b}.pkl'),'rb'))

def ncols(sec):
    meta = pickle.load(open(os.path.join(W5C1,'sector_meta.pkl'),'rb'))
    return len(meta['sec_cols'][sec])

class Echelon:
    """Sparse row echelon over GF(p); tracks pivot-creating source rows."""
    def __init__(self, ncols, p):
        self.NC = ncols; self.p = p
        self.pivot_of = np.full(ncols, -1, np.int32)
        self.colcnt = np.zeros(ncols, np.int32)
        self.pivots = []     # (idx int32, vals int64) normalized pivot rows
        self.pivcol = []     # pivot column per pivot
        self.pivrow = []     # source row index per pivot (-1 for ad-hoc)
        self.buf = np.zeros(ncols, np.int64)
        self.in_touched = np.zeros(ncols, bool)
        self.ops = 0
    def reduce(self, cols, vals):
        p = self.p
        buf, in_touched, pivot_of = self.buf, self.in_touched, self.pivot_of
        touched = []
        for c,a in zip(cols, vals):
            if not in_touched[c]: in_touched[c]=True; touched.append(c)
            buf[c] = (buf[c]+a) % p
        h = []
        for c in touched:
            if pivot_of[c] >= 0 and buf[c] != 0:
                h.append((int(pivot_of[c]), c))
        heapq.heapify(h)
        ops = 0
        while h:
            pi, c = heapq.heappop(h)
            a = int(buf[c] % p)
            if a == 0: continue
            idx, vals2 = self.pivots[pi]
            buf[idx] -= a*vals2
            buf[idx] %= p
            ops += len(idx)
            new = idx[(~in_touched[idx]) & (buf[idx] != 0)]
            if len(new):
                in_touched[new] = True
                touched.extend(new.tolist())
            pcs = idx[(pivot_of[idx] >= 0) & (buf[idx] != 0)]
            for c2 in pcs.tolist():
                heapq.heappush(h, (int(pivot_of[c2]), c2))
            buf[c] = 0
        self.ops += ops
        rem = [(c, int(buf[c])) for c in touched if buf[c] != 0]
        for c in touched: buf[c] = 0; in_touched[c] = False
        return rem
    def add_pivot(self, rem, src=-1):
        p = self.p
        c = min(rem, key=lambda ca: (self.colcnt[ca[0]], ca[0]))[0]
        idx = np.array([cc for cc,_ in rem], np.int32)
        vals = np.array([vv for _,vv in rem], np.int64) % p
        inv = pow(int(dict(rem)[c]), -1, p)
        vals = (vals*inv) % p
        self.pivots.append((idx, vals)); self.pivcol.append(c); self.pivrow.append(src)
        self.pivot_of[c] = len(self.pivots)-1
        self.colcnt[idx] += 1
        return c
    def rank(self): return len(self.pivots)
    def feed(self, rows, label='', t0=None, logevery=8000):
        """rows: list of (cols int array, vals int array) mod p. Returns S, C."""
        for ri in range(len(rows)):
            cols, vals = rows[ri]
            rem = self.reduce(cols.tolist(), vals.tolist())
            if rem: self.add_pivot(rem, src=ri)
            if (ri+1) % logevery == 0:
                el = time.monotonic()-t0 if t0 else 0
                nnz = [len(x[0]) for x in self.pivots[-2000:]] or [0]
                log(f'{label} row {ri+1}/{len(rows)} pivots={self.rank()} nnz_avg={sum(nnz)/len(nnz):.0f} max={max(nnz)} ops={self.ops/1e9:.2f}e9 el={el:.0f}s')
        return self.pivrow, self.pivcol

def to_modp(a, p):
    if isinstance(a, (int, np.integer)): return int(a) % p
    a = F(a)
    d = a.denominator % p
    assert d != 0, 'denominator divisible by p'
    return a.numerator % p * pow(d, -1, p) % p

def vec_modp(dct, p):
    cols = np.array(sorted(dct.keys()), np.int32)
    vals = np.array([to_modp(dct[int(c)], p) for c in cols], np.int64)
    keep = vals != 0
    return cols[keep], vals[keep]
'''
open('w5c2_lib.py','w').write(lib)
print("lib written", len(lib))