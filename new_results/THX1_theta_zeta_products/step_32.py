
import time, pickle, sys
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# ---------- fast sparse elimination with column index, exact Fraction and modular ----------
class ERow:
    __slots__ = ('d','h')
    def __init__(self, d, h=None):
        self.d = {k:v for k,v in d.items() if v}
        self.h = h or {}

def eliminate_indexed(rows, columns, track_history=True):
    """Gaussian elimination over Fractions. rows: list of ERow (d=vec, h=history).
    Returns (pivot_rows, pivots): pivot_rows have pivot coeff normalized to 1.
    Column index: col -> list of row ids currently containing it."""
    rows = [ERow(dict(r.d), dict(r.h) if track_history else None) for r in rows]
    colidx = {}
    for i, r in enumerate(rows):
        for k in r.d:
            colidx.setdefault(k, []).append(i)
    alive = [True]*len(rows)
    pivots = []
    pivot_rows = []
    for col in columns:
        ids = [i for i in colidx.get(col, []) if alive[i]]
        if not ids: continue
        # choose sparsest
        pi = min(ids, key=lambda i: len(rows[i].d))
        pr = rows[pi]
        pc = pr.d[col]
        if pc != 1:
            inv = 1/pc
            pr.d = {k: v*inv for k, v in pr.d.items()}
            if track_history: pr.h = {k: v*inv for k, v in pr.h.items()}
        # eliminate col from all other rows
        for i in ids:
            if i == pi: continue
            r = rows[i]
            f = r.d[col]
            nd = {}
            for k, v in r.d.items():
                nv = v - f*pr.d.get(k, Fr(0))
                if nv: nd[k] = nv
            for k, v in pr.d.items():
                if k not in r.d:
                    if -f*v: nd[k] = -f*v
            r.d = nd
            if track_history:
                nh = dict(r.h)
                for k, v in pr.h.items():
                    nv = nh.get(k, Fr(0)) - f*v
                    if nv: nh[k] = nv
                    elif k in nh: del nh[k]
                r.h = nh
            # update colidx
            for k in list(set(list(r.d.keys()) + [col])):
                pass
            # simpler: rebuild index lazily below
        alive[pi] = False
        pivots.append(col)
        pivot_rows.append(pr)
        # rebuild colidx for affected columns: columns of pr.d
        for k in list(pr.d.keys()):
            colidx[k] = [i for i in colidx.get(k, []) if alive[i]]
    return pivot_rows, pivots

# test at weight 6 vs rul1la
import rul1la
W = 6
rels = rul1.relations_at_weight(W, families=('DS','DUAL'))
cols = sorted(rul1.conv_words_of_weight(W))
t=time.time()
rows = [ERow(r, {i: Fr(1)}) for i,(fam,tag,r) in enumerate(rels)]
pr, piv = eliminate_indexed(rows, cols)
print("w6 exact indexed: rank", len(piv), "codim", len(cols)-len(piv), "time %.1fs"%(time.time()-t))
pickle.dump({'rank': len(piv), 'codim': len(cols)-len(piv)}, open('w6_test.pkl','wb'))
log("indexed exact eliminator tested at w6: codim %d" % (len(cols)-len(piv)))
