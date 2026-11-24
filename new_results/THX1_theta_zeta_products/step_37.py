
import sys, time, pickle
sys.dont_write_bytecode = True
import rul1, thx1_elim

data = pickle.load(open('rels7_dsdual.pkl','rb'))
P = 2**61-1
rowsD = []
for fam, tag, r in data:
    rowsD.append({k: (n * pow(d, -1, P)) % P for k, (n, d) in r.items()})
cols = sorted(rul1.conv_words_of_weight(7))
print("rows", len(rowsD), "cols", len(cols), flush=True)

# instrumented elimination with progress
def eliminate_mod_progress(rowsD, columns, p, print_every=100):
    rows = [dict(r) for r in rowsD]
    colidx = {}
    for i, r in enumerate(rows):
        for k in r:
            colidx.setdefault(k, []).append(i)
    alive = [True]*len(rows)
    pivots = []
    t0 = time.time()
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
        if len(pivots) % print_every == 0:
            avg = sum(len(r) for r in rows if r)/max(1,sum(1 for r in rows if r))
            print(f"pivots={len(pivots)} t={time.time()-t0:.0f}s avg_rowlen={avg:.1f} alive={sum(alive)}", flush=True)
    return pivots

piv = eliminate_mod_progress(rowsD, cols, P)
print("w7 DS+DUAL: rank", len(piv), "codim", len(cols)-len(piv), flush=True)
