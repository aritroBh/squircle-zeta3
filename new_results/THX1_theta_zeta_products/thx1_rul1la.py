"""rul1la.py -- exact rational linear algebra: sparse Fraction Gaussian elimination with history."""
from fractions import Fraction as Fr

class Row:
    __slots__ = ('d', 'h', 'pivot')
    def __init__(self, d, h=None, pivot=None):
        self.d = {k: v for k, v in d.items() if v}
        self.h = h or {}
        self.pivot = pivot

def eliminate(rows, columns):
    """RREF of Row list over columns (in given order). History tracked.
    Returns (rref, pivot_cols): rref = list of Row with .pivot set; each pivot col has coeff 1."""
    rows = [Row(dict(r.d), dict(r.h)) for r in rows]
    used = set()
    out = []
    for col in columns:
        cand = None
        for ri, r in enumerate(rows):
            if ri in used: continue
            if col in r.d:
                cand = ri; break
        if cand is None: continue
        r = rows[cand]
        pc = r.d[col]
        if pc != 1:
            r.d = {k: v/pc for k, v in r.d.items()}
            r.h = {k: v/pc for k, v in r.h.items()}
        for rj, r2 in enumerate(rows):
            if rj == cand: continue
            if col not in r2.d: continue
            f = r2.d[col]
            r2.d = {k: v - f*r.d.get(k, Fr(0)) for k, v in r2.d.items()}
            for k, v in r.d.items():
                if k not in r2.d:
                    r2.d[k] = -f*v
            r2.d = {k: v for k, v in r2.d.items() if v}
            r2.h = {k: v - f*r.h.get(k, Fr(0)) for k, v in r2.h.items()}
            for k, v in r.h.items():
                if k not in r2.h:
                    r2.h[k] = -f*v
            r2.h = {k: v for k, v in r2.h.items() if v}
        used.add(cand)
        r.pivot = col
        out.append(r)
    return out, [r.pivot for r in out]

def reduce_vector(vec, rref, hist_scale=1):
    """Reduce vector (dict) by rref rows (in order). Returns (residue, history).
    history maps original-row-id -> coeff: vec - residue = sum hist[id]*orig_row[id]."""
    v = {k: x for k, x in vec.items() if x}
    hist = {}
    for r in rref:
        pc = r.pivot
        if pc not in v: continue
        f = v[pc]
        v = {k: x - f*r.d.get(k, Fr(0)) for k, x in v.items()}
        for k, x in r.d.items():
            if k not in v: v[k] = -f*x
        v = {k: x for k, x in v.items() if x}
        for k, x in r.h.items():
            hist[k] = hist.get(k, Fr(0)) + f*x
    return {k: x for k, x in v.items() if x}, {k: x for k, x in hist.items() if x}
