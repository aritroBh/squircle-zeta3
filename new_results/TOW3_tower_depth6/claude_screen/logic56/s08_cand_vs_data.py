# s08: candidate (closed forms, my evaluator/e00 formula) vs the exact TOWN1 data: all keys d<=52, k=d..min(d+8,52), delta<=min(6,2d-4)
import sys, time
sys.dont_write_bytecode = True
from candlib import *
t0 = time.time()
cand = Cand(base_exprs())
H = load_data()
keys = sorted([k for k in H if k != '_meta'], key=lambda x: (x[0]+x[1], x))
tested = 0; bad = []; per = {}; rho_cells = {}; rho_bad = []
for (p, q) in keys:
    dd = p+q
    if dd < 2 or dd > 52: continue
    for k in range(dd, min(dd+8, 52)+1):
        for delta in range(0, min(6, 2*dd-4)+1):
            if not data_has(p, q, k, delta): continue
            if p == q and k < 2*p: continue
            a = cand.get(p, q, k, delta); b = data_get(p, q, k, delta)
            tested += 1; per[delta] = per.get(delta, 0)+1
            kind = 'rho(k=d or diag 2p)' if (k == dd) else 'k>d'
            if k == dd:
                rho_cells[delta] = rho_cells.get(delta, 0)+1
                if a != b: rho_bad.append((p, q, k, delta))
            if a != b: bad.append((p, q, k, delta))
print(f'candidate vs TOWN1 data: {tested} cells, {len(bad)} bad; per delta {dict(sorted(per.items()))}')
print(f'  of which initial-level cells (k=d; diag k=2p) per delta {dict(sorted(rho_cells.items()))}, bad {len(rho_bad)}')
print('  first bad:', bad[:5], round(time.time()-t0, 1), 's')
