src = '''"""q5_lll.py - integer-relation search over Q(sqrt5) via flint LLL (adapted from OCT5D3 q2_lll.py).

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
'''
open('q5_lll.py','w').write(src)
import sys
sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 800
from q5_lll import fit, report
phi = (1+mp.sqrt(5))/2
lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi2 = mp.pi**2; l2 = mp.log(2)
# Probe: any Q(sqrt5)-relation among ring4 + Li4(-1/phi), Li4(1/phi^2)?
ring4 = [('pi4', mp.pi**4), ('z3lf', z3*lf), ('z3l2', z3*l2), ('pi2lf2', pi2*lf**2), ('pi2lfl2', pi2*lf*l2),
         ('pi2l22', pi2*l2**2), ('lf4', lf**4), ('lf3l2', lf**3*l2), ('lf2l22', lf**2*l2**2), ('lfl23', lf*l2**3), ('l24', l2**4)]
# use zero-target trick: fit target=Li4(-1/phi) against ring4 + [Li4(1/phi^2)]
c = fit(mp.polylog(4,-1/phi), ring4 + [('Li4q', mp.polylog(4,1/phi**2))], cexp=int(800*0.90), nbest=3)
report("Li4(-1/phi) vs ring4 + Li4(1/phi^2)", c)
ring5 = [('z5', z5), ('pi2z3', pi2*z3), ('z3lf2', z3*lf**2), ('z3lfl2', z3*lf*l2), ('z3l22', z3*l2**2),
         ('pi4lf', mp.pi**4*lf), ('pi4l2', mp.pi**4*l2)]
c5 = fit(mp.polylog(5,-1/phi), ring5 + [('Li5q', mp.polylog(5,1/phi**2))], cexp=int(800*0.90), nbest=3)
report("Li5(-1/phi) vs ring5 + Li5(1/phi^2)", c5)
