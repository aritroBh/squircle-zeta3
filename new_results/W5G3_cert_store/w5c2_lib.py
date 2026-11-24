"""W5C2 shared library: sector data loading, mod-p sparse echelon with pivot-row tracking.
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
    def feed(self, rows, label='', t0=None, logevery=8000, stopcheck=None):
        """rows: list of (cols int array, vals int array) mod p. Returns S, C."""
        for ri in range(len(rows)):
            cols, vals = rows[ri]
            rem = self.reduce(cols.tolist(), vals.tolist())
            if rem: self.add_pivot(rem, src=ri)
            if (ri+1) % logevery == 0:
                el = time.monotonic()-t0 if t0 else 0
                nnz = [len(x[0]) for x in self.pivots[-2000:]] or [0]
                log(f'{label} row {ri+1}/{len(rows)} pivots={self.rank()} nnz_avg={sum(nnz)/len(nnz):.0f} max={max(nnz)} ops={self.ops/1e9:.2f}e9 el={el:.0f}s')
            if stopcheck and (ri+1) % stopcheck[0] == 0:
                if stopcheck[1](self):
                    log(f'{label} early stop at row {ri+1} pivots={self.rank()}')
                    return self.pivrow, self.pivcol, ri+1
        return self.pivrow, self.pivcol, None

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


class DenseEchelon:
    """Row echelon over GF(p) with DENSE pivot rows (fast at high density).
    Invariant: pivot row created at time t has zeros in all pivot columns of time < t
    (rem contains no pivot columns at creation). Reduction sweeps pivot columns in
    increasing creation time; eliminating pivot t only introduces nonzeros in pivot
    columns of time > t, so a rescan after each elimination is correct."""
    def __init__(self, ncols, p, maxpiv=None):
        self.NC = ncols; self.p = p
        self.maxpiv = maxpiv or ncols
        self.piv = np.zeros((self.maxpiv, ncols), np.int64)
        self.pivcol = np.full(self.maxpiv, -1, np.int64)
        self.pivot_of = np.full(ncols, -1, np.int64)
        self.pivrow = []
        self.r = 0
        self.elims = 0
    def reduce(self, cols, vals):
        p = self.p; NC = self.NC
        buf = np.zeros(NC, np.int64)
        buf[cols] = np.asarray(vals, np.int64) % p
        r = self.r
        if r:
            pc = self.pivcol[:r]
            elims = 0
            lo = 0
            while lo < r:
                nz = pc[lo:][buf[pc[lo:]] != 0]
                if nz.size == 0: break
                c = int(nz[0])           # smallest creation-time pivot column present
                t = int(self.pivot_of[c])
                a = int(buf[c])
                buf = (buf - a*self.piv[t]) % p
                elims += 1
                lo = t + 1               # pivot rows > t have zeros in pivot cols <= t
            self.elims += elims
        return buf
    def add_pivot(self, buf, src=-1):
        nz = np.nonzero(buf)[0]
        assert nz.size
        c = int(nz[0])
        lead = int(buf[c]) % self.p
        inv = pow(lead, -1, self.p)
        self.piv[self.r] = (buf*inv) % self.p
        self.pivcol[self.r] = c
        self.pivot_of[c] = self.r
        self.pivrow.append(src)
        self.r += 1
        return c
    def feed(self, rows, label='', t0=None, logevery=4000, stopcheck=None):
        for ri in range(len(rows)):
            cols, vals = rows[ri]
            buf = self.reduce(cols, vals)
            if buf.any(): self.add_pivot(buf, src=ri)
            if (ri+1) % logevery == 0:
                el = time.monotonic()-t0 if t0 else 0
                log(f'{label} row {ri+1}/{len(rows)} pivots={self.r} elims={self.elims/1e6:.1f}e6 el={el:.0f}s')
            if stopcheck and (ri+1) % stopcheck[0] == 0:
                if stopcheck[1](self):
                    log(f'{label} early stop at row {ri+1} pivots={self.r}')
                    return ri+1
        return None


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


class RrefEchelon:
    """Reduced row-echelon over GF(p) with DENSE basis rows (int64).
    Invariant: basis in RREF - pivot row t has 1 at pivcol[t] and 0 at all other pivot cols.
    reduce() is cascade-free: subtract hit pivot rows in ONE gathered matmul.
    add_pivot() back-eliminates the new pivot column from all existing rows (batched 2D op)."""
    def __init__(self, ncols, p):
        self.NC = ncols; self.p = p
        self.piv = np.zeros((ncols, ncols), np.int64)
        self.pivcol = np.full(ncols, -1, np.int64)
        self.pivot_of = np.full(ncols, -1, np.int64)
        self.pivrow = []
        self.r = 0
        self.ops = 0
    def reduce(self, cols, vals):
        p = self.p; NC = self.NC
        buf = np.zeros(NC, np.int64)
        buf[cols] = np.asarray(vals, np.int64) % p
        r = self.r
        if r:
            pc = self.pivcol[:r]
            hit = pc[buf[pc] != 0]
            if hit.size:
                coeffs = buf[hit] % p
                rows = self.piv[self.pivot_of[hit]]
                buf = (buf - coeffs @ rows) % p
                self.ops += int(hit.size)*NC
        return buf
    def add_pivot(self, buf, src=-1):
        p = self.p; NC = self.NC; r = self.r
        nz = np.nonzero(buf)[0]
        if nz.size == 0: return -1
        c = int(nz[0])
        lead = int(buf[c]) % p
        P = (buf*pow(lead, -1, p)) % p
        # back-eliminate c from existing rows
        colv = self.piv[:r, c]
        rw = np.nonzero(colv)[0]
        if rw.size:
            coeffs = colv[rw] % p
            self.piv[:r][rw] = (self.piv[:r][rw] - coeffs[:, None]*P[None, :]) % p
            self.ops += int(rw.size)*NC
        self.piv[r] = P
        self.pivcol[r] = c
        self.pivot_of[c] = r
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
                log(f'{label} row {ri+1}/{len(rows)} pivots={self.r} ops={self.ops/1e9:.2f}e9 el={el:.0f}s')
            if stopcheck and (ri+1) % stopcheck[0] == 0:
                if stopcheck[1](self):
                    log(f'{label} early stop at row {ri+1} pivots={self.r}')
                    return self.pivrow, list(self.pivcol[:self.r]), ri+1
        return self.pivrow, list(self.pivcol[:self.r]), None


class LeanEchelon:
    """Sparse row echelon over GF(p), monotone time-sweep reduction (no heap).
    Invariant: pivot row created at time t has zeros in all pivot cols of time < t,
    so during reduce() only pivot times > current can become newly nonzero; a single
    upward sweep with a dirty-mark array eliminates each pivot col at most once per row.
    buf entries kept unreduced between elims (bounded: < p + Nelims*p^2 << int64);
    multiplier always reduced mod p first; pivot col zeroed explicitly after its elim."""
    def __init__(self, ncols, p):
        self.NC = ncols; self.p = p
        self.piv_idx = []; self.piv_val = []
        self.pivcol = np.full(ncols, -1, np.int64)
        self.pivot_of = np.full(ncols, -1, np.int64)
        self.colcnt = np.zeros(ncols, np.int64)
        self.pivrow = []
        self.r = 0
        self.elims = 0
    def reduce(self, cols, vals):
        p = self.p; NC = self.NC
        buf = np.zeros(NC, np.int64)
        buf[cols] = np.asarray(vals, np.int64) % p
        r = self.r
        if r == 0:
            return buf
        dirty = np.zeros(r, bool)
        ts0 = self.pivot_of[cols]
        init = ts0[(ts0 >= 0) & (buf[cols] != 0)]
        if init.size:
            dirty[init] = True
        cur = 0
        elims = 0
        while cur < r:
            seg = dirty[cur:]
            if not seg.any(): break
            t = cur + int(np.argmax(seg))
            c = int(self.pivcol[t])
            a = int(buf[c] % p)
            if a:
                idx = self.piv_idx[t]; vals2 = self.piv_val[t]
                buf[idx] -= a*vals2
                buf[c] = 0
                ts = self.pivot_of[idx]
                newt = ts[(ts > t) & (buf[idx] != 0)]
                if newt.size:
                    dirty[newt] = True
                elims += 1
            dirty[t] = False
            cur = t + 1
        self.elims += elims
        buf %= p
        return buf
    def add_pivot(self, buf, src=-1):
        p = self.p
        nz = np.nonzero(buf)[0]
        if nz.size == 0: return -1
        c = int(min(nz.tolist(), key=lambda cc: (self.colcnt[cc], cc)))
        lead = int(buf[c] % p)
        inv = pow(lead, -1, p)
        self.piv_idx.append(nz.astype(np.int32))
        self.piv_val.append((buf[nz] % p)*inv % p)
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
                log(f'{label} row {ri+1}/{len(rows)} pivots={self.r} elims={self.elims/1e6:.1f}e6 nnz_avg={sum(nnz)/len(nnz):.0f} max={max(nnz)} el={el:.0f}s')
            if stopcheck and (ri+1) % stopcheck[0] == 0:
                if stopcheck[1](self):
                    log(f'{label} early stop at row {ri+1} pivots={self.r}')
                    return self.pivrow, list(self.pivcol[:self.r]), ri+1
        return self.pivrow, list(self.pivcol[:self.r]), None
