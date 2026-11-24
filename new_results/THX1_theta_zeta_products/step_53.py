
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1
import sympy as sp

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

def elim_mod_hist(rowsD, columns, p, cpfile):
    rows = [(dict(d), {i: 1}) for i, d in enumerate(rowsD)]
    colidx = {}
    for i, (d, h) in enumerate(rows):
        for k in d:
            colidx.setdefault(k, []).append(i)
    alive = [True]*len(rows)
    pivots = []; pivot_rows = []
    t0 = time.time()
    for col in columns:
        ids = [i for i in colidx.get(col, []) if alive[i] and col in rows[i][0]]
        if not ids: continue
        pi = min(ids, key=lambda i: len(rows[i][0]))
        pd, ph = rows[pi]
        pc = pd[col]
        if pc != 1:
            inv = pow(pc, -1, p)
            for k in list(pd.keys()): pd[k] = pd[k]*inv % p
            for k in list(ph.keys()): ph[k] = ph[k]*inv % p
        piv_cols = list(pd.keys())
        for i in ids:
            if i == pi or not alive[i] or col not in rows[i][0]: continue
            rd, rh = rows[i]
            f = rd[col] % p
            oldkeys = set(rd.keys())
            nd = {}
            for k, v in rd.items():
                nv = (v - f*pd.get(k, 0)) % p
                if nv: nd[k] = nv
            for k, v in pd.items():
                if k not in oldkeys:
                    nv = (-f*v) % p
                    if nv: nd[k] = nv
            nh = dict(rh)
            for k, v in ph.items():
                nv = (nh.get(k, 0) - f*v) % p
                if nv: nh[k] = nv
                elif k in nh: del nh[k]
            rows[i] = (nd, nh)
            for k in nd.keys():
                if k not in oldkeys:
                    colidx.setdefault(k, []).append(i)
        alive[pi] = False
        for k in piv_cols:
            if k in colidx:
                colidx[k] = [i for i in colidx[k] if alive[i]]
        pivots.append(col); pivot_rows.append((col, dict(pd), dict(ph)))
    pickle.dump({'pivots': pivots, 'pivot_rows': pivot_rows}, open(cpfile,'wb'))
    return len(pivots)

data = pickle.load(open('rels7_dsdual.pkl','rb'))
primes = [int(sp.prevprime(2**61)), int(sp.prevprime(2**60))]
cols7 = sorted(rul1.conv_words_of_weight(7))
for idx, P in enumerate(primes):
    t0=time.time()
    rowsD = [{k: (n * pow(d, -1, P)) % P for k, (n, d) in r.items()} for fam, tag, r in data]
    n = elim_mod_hist(rowsD, cols7, P, f'w7_hist_p{idx+2}.pkl')
    print(f"prime {P}: {n} pivots, {time.time()-t0:.0f}s", flush=True)
    log("w7 RREF+history prime %d done" % (idx+2))
