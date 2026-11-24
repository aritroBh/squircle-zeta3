"""r06: claim 1, numbers. Values of u.r (u = 1,7,3,5) for ALL 389 used lift1 rows on the two weight-1 branch rows,
computed word by word with the referee's own double-precision evaluator (rcommon.H), not with eval8.
Also the scaled value |u.r . V| / max|coeff| that W5G1's s02 reports (s02 printed |value| unscaled)."""
from rcommon import *
import pickle
L = 'r06.log'
W5 = words(5)
lift1 = pickle.load(open(T + '/W5C1_certificate_scoping/rels_lift1.pkl', 'rb'))
ex = pickle.load(open(HERE + '/r02_extract.pkl', 'rb'))['exact']
idxs = sorted({i for t, f, i, *_ in ex if f == 'lift1'})
res = {u: [] for u in (1, 7, 3, 5)}
for i in idxs:
    r = {W5[k]: c for k, c in lift1[i][1]}
    for u in (1, 7, 3, 5):
        res[u].append(abs(val(actv(r, u))))
for u in (1, 7, 3, 5):
    v = res[u]
    log(f'u={u}: {len(v)} rows; min |value| {min(v):.3g}, max |value| {max(v):.3g}', L)
pickle.dump(res, open(HERE + '/r06_values.pkl', 'wb'))
