import sys, time, pickle, json
sys.dont_write_bytecode = True
sys.path.insert(0, '../HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from cor13c_common import target4
from fractions import Fraction as Fr
from flint import fmpq_mat
t0 = time.monotonic()
lower = {}
for n in (1,2,3):
    rows, _ = cmzv6.make_relations(n, lower)
    lower[n] = rows
rows4, _ = cmzv6.make_relations(4, lower)
cols = cmzv6.words(4)
prow = pickle.load(open('w4_pivots_13.pkl','rb'))['prow']
sel = [rows4[i] for i in prow]
y = pickle.load(open('w4_y_13.pkl','rb'))
T = target4()
ids = {w:i for i,w in enumerate(cols)}
r = len(sel); m = len(cols)
print(f'r={r} m={m}; building fmpq_mat ({time.monotonic()-t0:.1f}s)', flush=True)
S = fmpq_mat(r, m)
for i, row in enumerate(sel):
    for w, c in row:
        S[i, ids[w]] = c
Y = fmpq_mat(r, 1)
for i, c in enumerate(y):
    Y[i, 0] = fmpq_mat(1,1,[str(c)])[0,0]
print(f'built ({time.monotonic()-t0:.1f}s); multiplying...', flush=True)
res = (Y.transpose() * S)
print(f'multiplied ({time.monotonic()-t0:.1f}s); comparing...', flush=True)
bad = []
for w, j in ids.items():
    got = Fr(str(res[0, j]))
    if got != T.get(w, Fr(0)):
        bad.append((w, got, T.get(w, Fr(0))))
print(f'exact verification: mismatches = {len(bad)} ({time.monotonic()-t0:.1f}s)')
assert not bad, bad[:5]
print('WEIGHT-4 EXACT CERTIFICATE VERIFIED')
out = dict(cols=[list(w) for w in cols], prow=prow, y=[str(c) for c in y])
json.dump(out, open('cert13_w4.json','w'))
import os
print('cert13_w4.json', os.path.getsize('cert13_w4.json'), 'bytes')
with open('progress.log','a') as f:
    f.write(f'step_30: weight-4 exact certificate verified on all 1764 columns ({time.monotonic()-t0:.0f}s)\n')
