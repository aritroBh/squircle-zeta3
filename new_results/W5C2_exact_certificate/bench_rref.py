"""Benchmark RrefEchelon on real pp sector rows (length-sorted feed)."""
import sys, os, time, pickle
sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import numpy as np
from w5c2_lib import *
sec = (1,1); p = P1; N = int(sys.argv[1]) if len(sys.argv)>1 else 12000
t0 = time.monotonic()
rows = load_sector_rows(sec)
rrows = []
for fam, i, r in rows[:N]:
    rrows.append((np.array([c for c,_ in r], np.int32), np.array([int(a)%p for _,a in r], np.int64)))
NC = ncols(sec)
E = RrefEchelon(NC, p)
E.feed(rrows, label='rrefbench', t0=t0, logevery=2000)
print(f'BENCH {N}: pivots={E.r} ops={E.ops/1e9:.2f}e9 time={time.monotonic()-t0:.0f}s')
D = load_D(sec,'rat'); cols, vals = vec_modp(D, p)
t1=time.monotonic(); buf=E.reduce(cols,vals)
print(f'D_rat reduce {time.monotonic()-t1:.2f}s residual={int(buf.any())}')
