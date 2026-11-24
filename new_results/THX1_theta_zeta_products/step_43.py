
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_elim

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

data = pickle.load(open('rels7_dsdual.pkl','rb'))
rels = [(fam, tag, {k: Fr(n, d) for k, (n, d) in r.items()}) for fam, tag, r in data]
cols7 = sorted(rul1.conv_words_of_weight(7))
pickle.dump(cols7, open('cols7.pkl','wb'))

# exact elimination with history, time-capped, checkpointing pivot rows
def eliminate_exact_cp(rows_in, columns, cap=None, cp_every=100, cpfile='w7_exact_cp.pkl'):
    rows = [thx1_elim.ERow(dict(r.d), dict(r.h)) for r in rows_in]
    colidx = {}
    for i, r in enumerate(rows):
        for k in r.d:
            colidx.setdefault(k, []).append(i)
    alive = [True]*len(rows)
    pivots = []
    pivot_rows = []
    t0 = time.time()
    for col in columns:
        ids = [i for i in colidx.get(col, []) if alive[i] and col in rows[i].d]
        if not ids: continue
        pi = min(ids, key=lambda i: len(rows[i].d))
        pr = rows[pi]
        pc = pr.d[col]
        if pc != 1:
            inv = 1/pc
            pr.d = {k: v*inv for k, v in pr.d.items()}
            pr.h = {k: v*inv for k, v in pr.h.items()}
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
        if len(pivots) % cp_every == 0:
            pickle.dump({'pivots': pivots, 'pivot_rows': [(p.d, p.h) for p in pivot_rows]}, open(cpfile,'wb'))
            print(f"piv={len(pivots)} t={time.time()-t0:.0f}s", flush=True)
        if cap and time.time()-t0 > cap:
            pickle.dump({'pivots': pivots, 'pivot_rows': [(p.d, p.h) for p in pivot_rows]}, open(cpfile,'wb'))
            print("CAP at", len(pivots))
            return None, None
    pickle.dump({'pivots': pivots, 'pivot_rows': [(p.d, p.h) for p in pivot_rows]}, open(cpfile,'wb'))
    return pivot_rows, pivots

rows = [thx1_elim.ERow(r, {i: Fr(1)}) for i, (fam, tag, r) in enumerate(rels)]
pr, piv = eliminate_exact_cp(rows, cols7, cap=1150, cp_every=150)
if pr is not None:
    print("w7 exact RREF with history: rank", len(piv), "codim", len(cols7)-len(piv))
    log("w7 exact RREF complete: codim %d" % (len(cols7)-len(piv)))
else:
    log("w7 exact RREF capped; will switch to modular certificates")
