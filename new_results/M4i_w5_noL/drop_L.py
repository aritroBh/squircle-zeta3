"""M4f: can L(chi8,5) be dropped entirely? Exact test in the coordinate system:
allowed columns = kv (70) + u_S (14), NO L/s2L. Does tbar lie in the span? At what height?
"""
import pickle, sys, itertools
from fractions import Fraction as F
from math import lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import json
import nice_search as NS

def log(s): print(s, flush=True)
results = pickle.loads((HERE/'subset_exact.pkl').read_bytes())
# rebuild coords exactly as subset_exact.py
oc = json.loads((HERE/'overcomplete_1200.json').read_text())
deps86 = [[F(c) for c in dep['vec'][1:]] for dep in oc['deps']]
rel = oc['target_rels'][0]; c0 = F(rel['vec'][0]); relv = [F(c) for c in rel['vec'][1:]]
kept = oc['kept']
kidx = []
for nm in kept:
    if nm == 'L': kidx.append(70)
    elif nm == 's2L': kidx.append(71)
    elif nm.startswith('s2kv'): kidx.append(35+int(nm[4:]))
    elif nm.startswith('kv'): kidx.append(int(nm[2:]))
    elif nm.startswith('s2u'): kidx.append(79+int(nm[3:]))
    else: kidx.append(72+int(nm[1:]))
Bcols = [[F(int(i==j)) for i in range(86)] for j in kidx]
M = [[Bcols[j][i] for j in range(72)] + [deps86[d][i] for d in range(14)] for i in range(86)]
rhs = [[F(int(i==j)) for j in range(86)] for i in range(86)]
mat, rank = flint.fmpq_mat([[str(c) for c in M[i]+rhs[i]] for i in range(86)]).rref()
inv = [[F(str(mat[i,86+j])) for j in range(86)] for i in range(86)]
coords = [[inv[i][c] for i in range(72)] for c in range(86)]
tbar = [F(-relv[i], c0) for i in range(72)]
# find positions of L, s2L in the 72-basis
iL = kidx.index(70); is2L = kidx.index(71)
log(f'tbar coords on L, s2L: {tbar[iL]}, {tbar[is2L]}')
# allowed: all kv (0..69) + u0..u6,s2u0..s2u6 (72..85); NOT L/s2L
allow = list(range(70)) + list(range(72,86))
Mx = [[str(coords[c][i]) for c in allow] + [str(tbar[i])] for i in range(72)]
mat, rank = flint.fmpq_mat(Mx).rref()
pivots = [next((j for j in range(mat.ncols()) if mat[i,j]), None) for i in range(rank)]
ncol = len(allow)
ok = all(p is None or p < ncol for p in pivots)
log(f'target in span(kv, u0..u6) WITHOUT L: {ok}')
if ok:
    x0 = [F(0)]*ncol
    for i,p in enumerate(pivots):
        if p is not None and p < ncol: x0[p] = F(str(mat[i,ncol]))
    # height of the u-coefficients and kv coefficients (rough)
    ucoef = [(allow[j], x0[j]) for j in range(ncol) if x0[j] and allow[j]>=72]
    log(f'u-coefficients: {ucoef}')
    mx = max(abs(x.numerator) for x in x0 if x)
    dn = lcm(*[x.denominator for x in x0 if x])
    log(f'solution-vector max|num|={mx} lcm={dn} (atom expansion would be much bigger)')
