"""M4f: height of the overcomplete-LLL relation (S={0,2,3}), expanded to atoms,
plus optimizer machinery test. Reuses subset_exact.pkl.
"""
import json, pickle, sys
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nice_search as NS

results = pickle.loads((HERE/'subset_exact.pkl').read_bytes())
oc = json.loads((HERE/'overcomplete_1200.json').read_text())
kept = oc['kept']
rel = oc['target_rels'][0]
c0 = F(rel['vec'][0]); relv = [F(c) for c in rel['vec'][1:]]

# atom map identical to subset_exact.py
sol = NS.solve(list(NS.S0), 'h')
names = sol['names']; xp0, xq0 = sol['xp'], sol['xq']
K = json.loads((HERE/'nice_ksat_s0.json').read_text())['K']
pool = json.loads((HERE/'pool.json').read_text())['pool']
s0 = set(NS.S0)
extras = [e for e in pool if NS.W.Q(F(e['p']),F(e['q'])) not in s0]
NATOM = 146
def acol(c):
    v = [F(0)]*NATOM
    if c < 35:
        for i,a in enumerate(K[c]): v[i] = F(a)
    elif c < 70:
        for i,a in enumerate(K[c-35]): v[65+i] = F(a)
    elif c == 70: v[130] = F(1)
    elif c == 71: v[131] = F(1)
    elif c < 79:
        i = c-72; e = extras[i]
        v[132+i] = F(1)
        for j,a in enumerate(e['xp']): v[j] -= F(a)
        for j,a in enumerate(e['xq']): v[65+j] -= F(a)
    else:
        i = c-79; e = extras[i]
        v[139+i] = F(1)
        for j,a in enumerate(e['xq']): v[j] -= 2*F(a)
        for j,a in enumerate(e['xp']): v[65+j] -= F(a)
    return v
def kidx_of(nm):
    if nm == 'L': return 70
    if nm == 's2L': return 71
    if nm.startswith('s2kv'): return 35+int(nm[4:])
    if nm.startswith('kv'): return int(nm[2:])
    if nm.startswith('s2u'): return 79+int(nm[3:])
    return 72+int(nm[1:])
ib2 = names.index('Li5((3-2*sqrt2))'); ibp = names.index('Li5((-1+1*sqrt2))')
ibm = names.index('Li5((1-1*sqrt2))'); im1 = names.index('Li5(-1)'); iz5 = names.index('z5')
def reduced(v):
    v = list(v)
    for off in (0,65):
        a = v[ib2+off]
        if a: v[ibp+off] += 16*a; v[ibm+off] += 16*a; v[ib2+off] = F(0)
        a = v[im1+off]
        if a: v[iz5+off] += F(-15,16)*a; v[im1+off] = F(0)
    return v
def height(v):
    v = reduced(v)
    nz = [a for a in v if a]
    return max(abs(a.numerator) for a in nz), lcm(*[a.denominator for a in nz]), len(nz)

v = [F(0)]*NATOM
for i,a in enumerate(xp0): v[i] = a
for i,a in enumerate(xq0): v[65+i] = a
for coef, nm in zip(relv, kept):
    if not coef: continue
    a = F(-coef, c0)
    A = acol(kidx_of(nm))
    for s in range(NATOM):
        if A[s]: v[s] += a*A[s]
print('overcomplete LLL relation expanded (S={0,2,3}):')
print('  raw:', height(v)[0], height(v)[1], 'terms', height(v)[2])
# compare: the rref particular for S=(0,2,3)
r0 = [r for r in results if r['S']==(0,2,3)][0]
print('  rref particular:', height(r0['v0']))
# nullity-6 space: check the LLL point is in the affine space: v_lll - v0 in span(NSimg)?
import flint
diff = [a-b for a,b in zip(v, r0['v0'])]
NSimg = r0['NSimg']
M = [[str(NSimg[j][s]) for j in range(len(NSimg))]+[str(diff[s])] for s in range(NATOM)]
mat, rank = flint.fmpq_mat(M).rref()
pivots = [next((j for j in range(mat.ncols()) if mat[i,j]), None) for i in range(rank)]
print('  LLL point in rref affine space:', all(p is None or p < len(NSimg) for p in pivots))
