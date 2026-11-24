import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
from math import gcd, lcm
import gold5_lib as G
Fraction = Fr
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
nz = [z for z, cl in li5cl.items() if any(any(q) for q in cl)]
def q12(v):
    out = []
    for q in v: out += [q[0], q[1]]
    return tuple(out)
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
SQ5v = (Fr(0), Fr(1))
def cols_z(z):
    cl = li5cl[z]
    return [q12(cl), q12(tuple(nfmul(SQ5v, q) for q in cl))]
targets = {
 'gol1': ((Fr(84), Fr(-3)), (Fr(0), Fr(0)), (Fr(60), Fr(21)), (Fr(0), Fr(0)),
          (Fr(513), Fr(-63,2)), (Fr(1689), Fr(-129))),
 'gol3': ((Fr(-6), Fr(0)), (Fr(0), Fr(0)), (Fr(-21), Fr(3)), (Fr(0), Fr(0)),
          (Fr(-84), Fr(18)), (Fr(-357), Fr(87))),
}
from flint import fmpz_mat, fmpz
for tag, T in targets.items():
    cols = []
    for z in nz: cols += cols_z(z)
    colsT = [q12(T), q12(tuple(nfmul(SQ5v, q) for q in T))]
    allcols = cols + colsT
    # scale rows to integers
    M = [[c[r] for c in allcols] for r in range(12)]
    den = 1
    for row in M:
        for x in row:
            den = lcm(den, x.denominator)
    Mi = [[int(x * den) for x in row] for row in M]
    # kernel of 12x(2m+2) integer matrix via flint: use fmpz_mat.nullspace? compute via HNF/LLL.
    A = fmpz_mat(12, len(allcols))
    for r in range(12):
        for c in range(len(allcols)):
            A[r, c] = fmpz(Mi[r][c])
    # kernel basis: use sympy for exact nullspace
    import sympy as sp
    As = sp.Matrix(Mi)
    ns = As.nullspace()
    print(f"{tag}: kernel dim {len(ns)}")
    # build integer kernel lattice, LLL-reduce
    K = fmpz_mat(len(ns), len(allcols))
    for i, v in enumerate(ns):
        den_v = 1
        for x in v: den_v = lcm(den_v, x.q)
        for j in range(len(allcols)):
            K[i, j] = fmpz(int(v[j] * den_v))
    R = K.lll(delta=0.99, eta=0.51)
    # find short vectors with nonzero target component (last 2 coords)
    cands = []
    for i in range(len(ns)):
        row = [int(R[i, j]) for j in range(len(allcols))]
        t1, t2 = row[-2], row[-1]
        if t1 or t2:
            h = max(abs(x) for x in row)
            cands.append((h, i, row))
    cands.sort(key=lambda t: t[0])
    print(f"  target-supporting kernel vectors, best heights: {[c[0] for c in cands[:8]]}")
    h, i, row = cands[0]
    print(f"  best: h={h}, target comps ({row[-2]}, {row[-1]})")
