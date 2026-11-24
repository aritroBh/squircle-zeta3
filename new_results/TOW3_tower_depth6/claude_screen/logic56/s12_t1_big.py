# s12: main pointwise T1 run: all keys, d=2..60, k=d+1..d+10 (K=0..9), delta<=min(6,2d-4); candidate = closed forms only.
import sys, time
sys.dont_write_bytecode = True
from candlib import *
t0 = time.time()
cand = Cand(base_exprs())
tot = 0; bad = []; per = {}; types = {}
for dd in range(2, 61):
    tested, b, pr = t1_test(cand, [dd], kspan=10)
    tot += tested; bad += b
    for dl, c in pr.items(): per[dl] = per.get(dl, 0)+c
    for p in range(dd, (dd-1)//2, -1):
        q = dd-p
        typ = 'diag' if p == q else ('adjacent' if q == p-1 else 'nonadjacent')
        for dl in pr:
            types[(typ, dl)] = types.get((typ, dl), 0) + 10 - (1 if typ == 'diag' else 0) if False else types.get((typ, dl), 0)
    if dd in (5, 6, 7, 20, 40, 60): print(f'  through d={dd}: cells {tot}, bad {len(bad)} ({round(time.time()-t0,1)}s)', flush=True)
print(f'TOTAL pointwise T1 cells {tot}, bad {len(bad)}; per delta {dict(sorted(per.items()))}')
# smallest-degree detail
t_, b_, p_ = t1_test(cand, [5], kspan=10)
print('d=5 only (all keys (5,0),(4,1),(3,2); includes depth 5 and 6; adjacent p=3): cells', t_, 'bad', len(b_), 'per delta', dict(sorted(p_.items())))
