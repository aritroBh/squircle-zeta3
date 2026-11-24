"""thx1_elim.py -- sparse elimination: exact Fraction with history, and modular."""
from fractions import Fraction as Fr

class ERow:
    __slots__ = ('d','h')
    def __init__(self, d, h=None):
        self.d = {k:v for k,v in d.items() if v}
        self.h = h

def eliminate_exact(rows, columns, track_history=True):
    rows = [ERow(dict(r.d), dict(r.h) if track_history else None) for r in rows]
    colidx = {}
    for i, r in enumerate(rows):
        for k in r.d:
            colidx.setdefault(k, []).append(i)
    alive = [True]*len(rows)
    pivots = []
    pivot_rows = []
    for col in columns:
        ids = [i for i in colidx.get(col, []) if alive[i] and col in rows[i].d]
        if not ids: continue
        pi = min(ids, key=lambda i: len(rows[i].d))
        pr = rows[pi]
        pc = pr.d[col]
        if pc != 1:
            inv = 1/pc
            pr.d = {k: v*inv for k, v in pr.d.items()}
            if track_history: pr.h = {k: v*inv for k, v in pr.h.items()}
        piv_cols = list(pr.d.keys())
        for i in ids:
            if i == pi or not alive[i] or col not in rows[i].d: continue
            r = rows[i]
            f = r.d[col]
            oldkeys = set(r.d.keys())
            nd = {}
            for k, v in r.d.items():
                nv = v - f*pr.d.get(k, Fr(0))
                if nv: nd[k] = nv
            for k, v in pr.d.items():
                if k not in oldkeys:
                    nd[k] = -f*v
            r.d = nd
            if track_history:
                nh = dict(r.h)
                for k, v in pr.h.items():
                    nv = nh.get(k, Fr(0)) - f*v
                    if nv: nh[k] = nv
                    elif k in nh: del nh[k]
                r.h = nh
            for k in nd.keys():
                if k not in oldkeys:
                    colidx.setdefault(k, []).append(i)
        alive[pi] = False
        for k in piv_cols:
            if k in colidx:
                colidx[k] = [i for i in colidx[k] if alive[i]]
        pivots.append(col)
        pivot_rows.append(pr)
    return pivot_rows, pivots

def eliminate_mod(rowsD, columns, p, want_rows=False):
    """rowsD: list of dict col->int (mod p). Returns (pivots, rank) and optionally RREF rows."""
    rows = [dict(r) for r in rowsD]
    colidx = {}
    for i, r in enumerate(rows):
        for k in r:
            colidx.setdefault(k, []).append(i)
    alive = [True]*len(rows)
    pivots = []
    pivot_rows = []
    for col in columns:
        ids = [i for i in colidx.get(col, []) if alive[i] and col in rows[i]]
        if not ids: continue
        pi = min(ids, key=lambda i: len(rows[i]))
        pr = rows[pi]
        pc = pr[col]
        if pc != 1:
            inv = pow(pc, -1, p)
            for k in list(pr.keys()):
                pr[k] = (pr[k]*inv) % p
        piv_cols = list(pr.keys())
        for i in ids:
            if i == pi or not alive[i] or col not in rows[i]: continue
            r = rows[i]
            f = r[col] % p
            nd = {}
            for k, v in r.items():
                nv = (v - f*pr.get(k, 0)) % p
                if nv: nd[k] = nv
            for k, v in pr.items():
                if k not in r:
                    nv = (-f*v) % p
                    if nv: nd[k] = nv
            rows[i] = nd
            for k in nd.keys():
                if k not in r:
                    colidx.setdefault(k, []).append(i)
        alive[pi] = False
        for k in piv_cols:
            if k in colidx:
                colidx[k] = [i for i in colidx[k] if alive[i]]
        pivots.append(col)
        if want_rows: pivot_rows.append(pr)
    return pivots, pivot_rows
