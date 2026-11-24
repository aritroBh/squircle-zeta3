import sys, time, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, '../HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from cor13c_common import target4
import exact_linalg as XL
t0 = time.monotonic()
lower = {}
for n in (1,2,3):
    rows, _ = cmzv6.make_relations(n, lower)
    lower[n] = rows
rows4, _ = cmzv6.make_relations(4, lower)
cols = cmzv6.words(4)
print(f'rows={len(rows4)} cols={len(cols)} ({time.monotonic()-t0:.1f}s)', flush=True)
prow, r = XL.select_independent_rows(rows4, cols)
print(f'rank={r} ({time.monotonic()-t0:.1f}s)', flush=True)
sel = [rows4[i] for i in prow]
pcols = XL.select_pivot_cols(sel, cols)
print(f'pcols={len(pcols)} ({time.monotonic()-t0:.1f}s)', flush=True)
pickle.dump(dict(prow=prow, pcols=pcols), open('w4_pivots_13.pkl','wb'))
T = target4()
y = XL.exact_solve(sel, pcols, cols, T)
print(f'solve done ({time.monotonic()-t0:.1f}s)', flush=True)
pickle.dump(y, open('w4_y_13.pkl','wb'))
hn = [max(len(str(abs(c.numerator))), len(str(c.denominator))) for c in y]
nz = sum(1 for c in y if c != 0)
print(f'multipliers: {len(y)} total, {nz} nonzero, max height {max(hn)} digits, median {sorted(hn)[len(hn)//2]}')
with open('progress.log','a') as f:
    f.write(f'step_29: weight-4 exact solve done, rank {r}, max height {max(hn)} digits\n')
