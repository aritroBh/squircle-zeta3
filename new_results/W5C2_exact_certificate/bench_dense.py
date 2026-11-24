"""Benchmark DenseEchelon on first N rows of pp sector."""
import sys, os, time, pickle
sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import numpy as np
from w5c2_lib import *
sec = (1,1); p = P1; N = int(sys.argv[1]) if len(sys.argv)>1 else 20000
t0 = time.monotonic()
rows = load_sector_rows(sec)
rrows = []
for fam, i, r in rows[:N]:
    cols = np.array([c for c,_ in r], np.int32)
    vals = np.array([int(a) % p for _,a in r], np.int64)
    rrows.append((cols, vals))
NC = ncols(sec)
E = DenseEchelon(NC, p)
E.feed(rrows, label='bench', t0=t0, logevery=2000)
print(f'BENCH {N} rows: pivots={E.r} elims={E.elims} time={time.monotonic()-t0:.0f}s')
# membership timing of D
D = load_D(sec, 'rat')
cols, vals = vec_modp(D, p)
t1 = time.monotonic()
buf = E.reduce(cols, vals)
print(f'D_rat reduce: {time.monotonic()-t1:.2f}s residual={int(buf.any())}')
