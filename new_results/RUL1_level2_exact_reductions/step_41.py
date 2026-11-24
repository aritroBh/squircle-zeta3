import pathlib
code = '''"""rul1la.py -- exact rational linear algebra: sparse Fraction Gaussian elimination with history."""
from fractions import Fraction as Fr

class Row:
    __slots__ = ('d', 'h')
    def __init__(self, d, h=None):
        self.d = {k: v for k, v in d.items() if v}   # vector
        self.h = h or {}                              # history: orig row id -> coeff

def eliminate(rows, columns):
    """Gauss elimination of Row list over the given column order.
    Returns (rref_rows, pivot_cols): rref_rows in row-echelon (reduced) form with history."""
    colidx = {c: i for i, c in enumerate(columns)}
    rows = [Row(dict(r.d), dict(r.h)) for r in rows]
    pivots = []   # (col, row)
    used = set()
    for ci, col in enumerate(columns):
        # find row with pivot in col among rows not used, with no pivot columns set
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
        # eliminate this col from all other rows
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
        pivots.append((col, cand))
    rref = [rows[cand] for (col, cand) in pivots]
    return rref, [col for (col, cand) in pivots]

def reduce_vector(vec, rref):
    """Reduce a vector (dict) by the rref rows; returns (residue dict, history dict)."""
    v = {k: x for k, x in vec.items() if x}
    hist = {}
    for r in rref:
        # pivot col = the one appearing... we stored rows reduced; pivot is the col where r has coeff 1 and others 0
        # find pivot: the column of r with coefficient 1 that no previous row has (we recompute: use any col present)
        pc = None
        for k in r.d:
            if r.d[k] == 1 and k in v:
                pc = k; break
        if pc is None:
            # check all cols
            for k in r.d:
                if k in v:
                    pc = k; break
        if pc is None: continue
        f = v[pc]/r.d[pc]
        v = {k: x - f*r.d.get(k, Fr(0)) for k, x in v.items()}
        for k, x in r.d.items():
            if k not in v: v[k] = -f*x
        v = {k: x for k, x in v.items() if x}
        for k, x in r.h.items():
            hist[k] = hist.get(k, Fr(0)) - f*x
    return v, {k: x for k, x in hist.items() if x}
'''
pathlib.Path('rul1la.py').write_text(code)
print("written rul1la.py")
