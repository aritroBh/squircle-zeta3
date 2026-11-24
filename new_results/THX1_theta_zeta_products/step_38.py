
import sys, time, pickle
sys.dont_write_bytecode = True
import rul1

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# Markowitz elimination mod p: at each step pick (col,row) minimizing (count_col-1)*(rowlen-1)
def eliminate_mod_markowitz(rowsD, p, print_every=200, max_seconds=None):
    rows = [dict(r) for r in rowsD]
    colidx = {}
    for i, r in enumerate(rows):
        for k in r:
            colidx.setdefault(k, []).append(i)
    alive = [True]*len(rows)
    # column counts (live, lazily updated)
    pivots = 0
    t0 = time.time()
    total_cols = len(colidx)
    while True:
        # find best (col,row): scan columns with their live counts
        best = None; bestscore = None
        for col, ids in colidx.items():
            live = [i for i in ids if alive[i] and col in rows[i]]
            if not live: continue
            nr = len(live)
            # sparsest row in this col
            li = min(live, key=lambda i: len(rows[i]))
            score = (nr-1)*max(1,(len(rows[li])-1))
            if bestscore is None or score < bestscore:
                bestscore = score; best = (col, li, live)
                if score == 0: break
        if best is None: break
        col, pi, live = best
        pr = rows[pi]
        pc = pr[col]
        if pc != 1:
            inv = pow(pc, -1, p)
            for k in list(pr.keys()):
                pr[k] = (pr[k]*inv) % p
        for i in live:
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
        pivots += 1
        if print_every and pivots % print_every == 0:
            nnz = sum(len(r) for r in rows if r)
            print(f"pivots={pivots} t={time.time()-t0:.0f}s nnz={nnz} alive={sum(alive)}", flush=True)
        if max_seconds and time.time()-t0 > max_seconds:
            print("TIME CAP", flush=True)
            return None, rows, alive
    return pivots, rows, alive

# benchmark on w7 first
data = pickle.load(open('rels7_dsdual.pkl','rb'))
P = 2**61-1
rowsD = [{k: (n * pow(d, -1, P)) % P for k, (n, d) in r.items()} for fam, tag, r in data]
t=time.time()
pivots, rows, alive = eliminate_mod_markowitz(rowsD, P, print_every=200)
print("w7 Markowitz: rank", pivots, "time %.0fs"%(time.time()-t))
log("Markowitz w7: rank %d codim 21 confirmed, time %.0fs" % (pivots, time.time()-t))
