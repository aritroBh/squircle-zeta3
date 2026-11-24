import sys, time, json
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from hexc2_common import target4
import exact_linalg as XL
from flint import fmpq_mat, fmpz_mat

t0 = time.monotonic()
lower = {}
for n in (1,2,3):
    rows, _ = cmzv6.make_relations(n, lower)
    lower[n] = rows
rows4, counts4 = cmzv6.make_relations(4, lower)
cols = cmzv6.words(4)
print(f'rows={len(rows4)} cols={len(cols)} ({time.monotonic()-t0:.1f}s)')
prow, r = XL.select_independent_rows(rows4, cols)
print(f'rank={r} ({time.monotonic()-t0:.1f}s)')
sel = [rows4[i] for i in prow]
pcols = XL.select_pivot_cols(sel, cols)
print(f'pcols={len(pcols)} ({time.monotonic()-t0:.1f}s)')
import pickle
pickle.dump(dict(prow=prow, pcols=pcols), open('w4_pivots.pkl','wb'))
T = target4()
y = XL.exact_solve(sel, pcols, cols, T)
print(f'solve done ({time.monotonic()-t0:.1f}s)')
pickle.dump(y, open('w4_y.pkl','wb'))
print('saved w4_pivots.pkl, w4_y.pkl')
