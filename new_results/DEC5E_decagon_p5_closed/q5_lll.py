"""q5_lll.py - integer-relation search over Q(sqrt5) via flint LLL (adapted from OCT5D3 q2_lll.py).

fit(target, atoms, dps, cexp): atoms = list of (name, mpf value). Returns list of
(height, residual, {name: (p,q)}, v0) with target = sum (p+q sqrt5) atom, sorted by height.
Lattice: vals = [target] + [v, sqrt5*v for each atom]; rows = I | C*vals (C = 10^cexp).
"""
import mpmath as mp
from flint import fmpz_mat
from fractions import Fraction as F

def fit(target, atoms, dps=None, cexp=None, nbest=6):
    if dps is None: dps = mp.mp.dps
    if cexp is None: cexp = int(dps*0.90)
    r5 = mp.sqrt(5)
    vals = [target]
    for nm, v in atoms:
        vals += [v, r5*v]
    n = len(vals)
    C = mp.mpf(10)**cexp
    rows = []
    for i in range(n):
        row = [0]*(n+1); row[i] = 1
        row[n] = int(mp.nint(C*vals[i]))
        rows.append(row)
    M = fmpz_mat(rows).lll()
    cands = []
    for i in range(n):
        v = [int(M[i,j]) for j in range(n)]
        if v[0] == 0: continue
        res = abs(sum(mp.mpf(v[j])*vals[j] for j in range(n)))
        maxc = max(abs(xx) for xx in v)
        rel = {}
        for j, (nm, _) in enumerate(atoms):
            p, q = v[1+2*j], v[2+2*j]
            if p or q: rel[nm] = (F(-p, v[0]), F(-q, v[0]))
        cands.append((maxc, res, rel, v[0]))
    cands.sort(key=lambda z: z[0])
    return cands[:nbest]

def report(label, cands):
    print("== %s: %d candidates ==" % (label, len(cands)))
    for maxc, res, rel, v0 in cands:
        print("  height %.3e  residual %s  denom %d" % (maxc, mp.nstr(mp.mpf(res), 3), abs(v0)))
