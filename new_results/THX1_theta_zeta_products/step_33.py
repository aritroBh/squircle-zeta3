
import time, pickle, sys
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

class ERow:
    __slots__ = ('d','h')
    def __init__(self, d, h=None):
        self.d = {k:v for k,v in d.items() if v}
        self.h = h

def eliminate_indexed(rows, columns, track_history=True):
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
            if i == pi: continue
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
        # prune dead from index for pivot row's columns
        for k in piv_cols:
            if k in colidx:
                colidx[k] = [i for i in colidx[k] if alive[i]]
        pivots.append(col)
        pivot_rows.append(pr)
    return pivot_rows, pivots

W = 6
rels = rul1.relations_at_weight(W, families=('DS','DUAL'))
cols = sorted(rul1.conv_words_of_weight(W))
t=time.time()
rows = [ERow(r, {i: Fr(1)}) for i,(fam,tag,r) in enumerate(rels)]
pr, piv = eliminate_indexed(rows, cols)
print("w6 exact indexed: rank", len(piv), "codim", len(cols)-len(piv), "time %.1fs"%(time.time()-t))
log("indexed exact eliminator w6: codim %d (matches DG 13)" % (len(cols)-len(piv)))
