# s05: POINTWISE T1 test of the candidate tower (closed forms only) -- all keys, d=2..40, k=d+1..d+6, delta<=min(6,2d-4)
import sys, time
sys.dont_write_bytecode = True
from candlib import *
t0 = time.time()
cand = Cand(base_exprs())
print('compiled families', round(time.time()-t0, 1), flush=True)
# per-degree report
tot = 0; allbad = []; per_all = {}
mind = {}
for dd in range(2, 41):
    tested, bad, per = t1_test(cand, [dd])
    tot += tested; allbad += bad
    for dl, c in per.items():
        per_all[dl] = per_all.get(dl, 0)+c
        if dl not in mind: mind[dl] = dd
    if dd in (4, 5, 6, 7, 8, 10, 15, 20, 30, 40) or bad:
        print(f'  d={dd}: cells {tested}, bad {len(bad)}  ({round(time.time()-t0,1)}s)', flush=True)
        if bad: print('    first bad:', bad[:6])
print(f'TOTAL T1 cells {tot}, bad {len(allbad)}; cells per delta {per_all}; smallest degree d at which depth delta exists: {mind}')
