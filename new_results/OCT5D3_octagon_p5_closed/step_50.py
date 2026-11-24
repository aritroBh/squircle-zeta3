import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
src = '''"""q2_lll.py — integer-relation search over Q(sqrt2) via flint LLL.

fit(target, atoms, dps, cexp=None): atoms = list of (name, mpf value). Returns list of
relations (height, residual, {name: (p,q)}) with target = sum (p+q sqrt2) atom, sorted by height.
The lattice: vals = [target] + [v, r2*v for each atom]; rows = I | C*vals (C = 10^cexp).
"""
import mpmath as mp
from flint import fmpz_mat
from fractions import Fraction as F

def fit(target, atoms, dps=None, cexp=None, nbest=6):
    if dps is None: dps = mp.mp.dps
    if cexp is None: cexp = int(dps*0.90)
    r2 = mp.sqrt(2)
    vals = [target]
    for nm, v in atoms:
        vals += [v, r2*v]
    n = len(vals)
    C = mp.mpf(10)**cexp
    rows = []
    for i in range(n):
        row = [0]*(n+1); row[i] = 1
        row[n] = int(mp.nint(C*vals[i]))
        rows.append(row)
    M = fmpz_mat(rows).lll()
    cands = []
    seen = set()
    for i in range(n):
        v = [int(M[i,j]) for j in range(n)]
        if v[0] == 0: continue
        res = abs(sum(mp.mpf(v[j])*vals[j] for j in range(n)))
        maxc = max(abs(xx) for xx in v)
        key = tuple((x//2 if all(y % 2 == 0 for y in v) else x) for x in v)  # crude dedupe aid
        rel = {}
        for j, (nm, _) in enumerate(atoms):
            p, q = v[1+2*j], v[2+2*j]
            if p or q: rel[nm] = (F(-p, v[0]), F(-q, v[0]))
        cands.append((maxc, res, rel, v[0]))
    cands.sort(key=lambda z: z[0])
    return cands[:nbest]

def report(label, cands, digits_resid=None):
    print("== %s: %d candidates ==" % (label, len(cands)))
    for maxc, res, rel, v0 in cands:
        print("  height %.3e  residual %s  denom %d" % (maxc, mp.nstr(res, 3), abs(v0)))
'''
open("q2_lll.py","w").write(src)
print("written")

# quick test on W1, W2
import sys, pickle
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 300
from oct5d3_atoms import base_atoms, W1_BASIS, W2_BASIS, prod
from q2_lll import fit, report
A = base_atoms()
W = pickle.load(open("wparts_1650.pkl","rb"))
W1 = mp.mpf(W['W1']); W2 = mp.mpf(W['W2'])
at1 = [(nm, prod(A, sp)) for nm, sp in W1_BASIS]
c1 = fit(W1, at1, cexp=280)
report("W1", c1)
for h, r, rel, v0 in c1[:2]:
    print("   ", rel)
at2 = [(nm, prod(A, sp)) for nm, sp in W2_BASIS]
c2 = fit(W2, at2, cexp=280)
report("W2", c2)
for h, r, rel, v0 in c2[:2]:
    print("   ", rel)
