
import sys, time, pickle
sys.dont_write_bytecode = True
import rul1

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

t=time.time()
rels8 = rul1.relations_at_weight(8, families=('DS','DUAL'))
print("build w8: %.1fs rels=%d"%(time.time()-t, len(rels8)), flush=True)
data = [(fam, tag, {k: (v.numerator, v.denominator) for k, v in r.items()}) for fam, tag, r in rels8]
pickle.dump(data, open('rels8_dsdual.pkl','wb'))
cols8 = sorted(rul1.conv_words_of_weight(8))
print("cols", len(cols8), flush=True)

P = 2**61-1
rowsD = [{k: (n * pow(d, -1, P)) % P for k, (n, d) in r.items()} for fam, tag, r in data]

def eliminate_mod_simple(rowsD, columns, p, print_every=250, cap=None):
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
        if print_every and len(pivots)%print_every==0:
            nnz = sum(len(r) for r in rows if r)
            print(f"piv={len(pivots)} t={time.time()-t0:.0f}s nnz={nnz} alive={sum(alive)}", flush=True)
        if cap and time.time()-t0 > cap:
            print("CAP"); return None
    return pivots

piv = eliminate_mod_simple(rowsD, cols8, P, cap=1200)
if piv:
    print("w8 DS+DUAL: rank", len(piv), "codim", len(cols8)-len(piv), flush=True)
    log("w8 DS+DUAL codim %d (DG 34)" % (len(cols8)-len(piv)))
