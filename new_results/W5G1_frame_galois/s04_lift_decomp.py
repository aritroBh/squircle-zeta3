"""s04: EXACT decomposition, for every lift row used by W5C2 and u in {3,5}:
   u.b_j = sum_i t_i b_i + e sh y   (b_i: proved base of weight n, e=(2)-(6), y in Q^{words(n-1)}),
verified in exact rational arithmetic. Writes lift_decomp.pkl {(n, j, u): (t as {i: Fraction}, y as dict)}."""
from common import *
import collections
used = pickle.load(open(HERE + '/used_rows.pkl', 'rb'))
R = {f: pickle.load(open(W5C1 + f'/rels_{f}.pkl', 'rb')) for f in FAM_N}
B = bases()
need = collections.defaultdict(set)
for tag, L in used.items():
    for f, i, *_ in L:
        if f in FAM_N: need[FAM_N[f]].add(R[f][i][0][1])
log('base rows needed per weight: ' + str({n: len(v) for n, v in need.items()}), 's04.log')
out = {}
for n in (1, 2, 3, 4):
    t0 = time.time()
    cols = words(n)
    erows = [{w2: F(c) for w2, c in mul(E1, {w: 1}).items()} for w in words(n - 1)] if n > 1 else [dict(E1)]
    gens = B[n] + erows
    keys = [(j, u) for j in sorted(need[n]) for u in (3, 5)]
    rhss = [actv(B[n][j], u) for j, u in keys]
    xs, rk = exact_combination(gens, rhss, cols, tag=f'w{n}')
    nb = len(B[n]); bad = 0; nzy = 0
    for (j, u), x in zip(keys, xs):
        assert x is not None, (n, j, u)
        t = {i: x[i] for i in range(nb) if x[i]}
        y = {}
        ww = words(n - 1) if n > 1 else [()]
        for k, wv in enumerate(ww):
            if x[nb + k]: y[wv] = x[nb + k]
        if y: nzy += 1
        out[(n, j, u)] = (t, y)
    log(f'weight {n}: rank(gens)={rk}; decomposed {len(keys)} (j,u) pairs EXACTLY; with nonzero e-part: {nzy}; {time.time()-t0:.0f}s', 's04.log')
pickle.dump(out, open(HERE + '/lift_decomp.pkl', 'wb'))
log('saved lift_decomp.pkl', 's04.log')
