
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, rul1la, rul1_named

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

t0 = time.time()
log("rul1_step02: build relations (DS+DUAL), RREF with history, named-basis inverse; save rul1_data.pkl")

data = {}
for w in range(1, 6):
    tw = time.time()
    rels = rul1.relations_at_weight(w, families=('DS','DUAL'))
    cols = sorted(rul1.conv_words_of_weight(w))
    rows = [rul1la.Row(r, {i: Fr(1)}) for i, (fam, tag, r) in enumerate(rels)]
    rref, piv = rul1la.eliminate(rows, cols)
    free = [c for c in cols if c not in piv]
    # named combos -> free coordinates
    M = []          # rows: named constants
    hist_named = []
    for nm, combo in rul1_named.NAMED[w]:
        res, h = rul1la.reduce_vector(combo, rref)
        M.append(res)
        hist_named.append(h)
    # matrix named x free; check square & invertible
    fnames = [nm for nm, _ in rul1_named.NAMED[w]]
    assert len(free) == len(fnames), (w, len(free), len(fnames))
    import sympy as sp
    Ms = sp.Matrix([[row.get(f, Fr(0)) for f in free] for row in M])
    det = Ms.det()
    Minv = Ms.inv()
    data[w] = {
        'rels': [(fam, tag, r) for fam, tag, r in rels],
        'cols': cols, 'rref': [(r.d, r.h, r.pivot) for r in rref],
        'free': free, 'M': M, 'Minv': [[Minv[i,j] for j in range(len(free))] for i in range(len(free))],
        'hist_named': hist_named, 'fnames': fnames,
        'det': det,
    }
    log("w=%d: rels %d rank %d codim %d det(named)=%s  (%.1fs)" % (w, len(rels), len(piv), len(free), det, time.time()-tw))

with open('rul1_data.pkl','wb') as f:
    pickle.dump(data, f)
log("saved rul1_data.pkl, total %.1fs" % (time.time()-t0))
