import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
from math import lcm
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
def solve_lin_full(cols, T):
    n, m = len(T), len(cols)
    A = [[cols[c][r] for c in range(m)] + [T[r]] for r in range(n)]
    R = 0; piv = []
    for c in range(m):
        p = next((r for r in range(R, n) if A[r][c]), None)
        if p is None: continue
        A[R], A[p] = A[p], A[R]
        q = A[R][c]; A[R] = [x / q for x in A[R]]
        for r in range(n):
            if r != R and A[r][c]:
                A[r] = [a - b * A[R][c] for a, b in zip(A[r], A[R])]
        piv.append(c); R += 1
        if R == n: break
    for r in range(R, n):
        if A[r][m]: return None
    free = [c for c in range(m) if c not in piv]
    part = [Fr(0)] * m
    for i, c in enumerate(piv): part[c] = A[i][m]
    ker = []
    for f in free:
        v = [Fr(0)] * m; v[f] = Fr(1)
        for i, c in enumerate(piv): v[c] = -A[i][f]
        ker.append(v)
    return part, ker
from flint import fmpz_mat, fmpz
for tag, T in targets.items():
    cols = []
    for z in nz: cols += cols_z(z)
    res = solve_lin_full(cols, q12(T))
    part, ker = res
    # scale to integers
    den = 1
    for x in part: den = lcm(den, x.denominator)
    for v in ker:
        for x in v: den = lcm(den, x.denominator)
    m = len(part)
    # lattice rows: (ker_i*den, 0), and target vector (part*den, den) with marker
    rows = []
    for v in ker:
        rows.append([int(x * den) for x in v] + [0])
    L = fmpz_mat(len(rows), m + 1)
    for i, row in enumerate(rows):
        for j, x in enumerate(row): L[i, j] = fmpz(x)
    R = L.lll(delta=0.99, eta=0.51)
    # reduced kernel basis; now greedy reduce part
    pv = [int(x * den) for x in part] + [den]
    # CVP approx: use reduced basis, Babai rounding in exact arithmetic over Q
    import sympy as sp
    B = [[int(R[i, j]) for j in range(m + 1)] for i in range(len(rows))]
    # sort by norm
    def nrm(v): return sum(x*x for x in v)
    B.sort(key=nrm)
    cur = list(pv)
    for it in range(3):
        for b in B:
            # find integer t minimizing |cur - t*b|
            bb = nrm(b)
            if bb == 0: continue
            t = round(sum(c*x for c, x in zip(cur, b)) / bb)
            if t:
                cur = [c - t*x for c, x in zip(cur, b)]
    # cur[m] should still be den; coefficients = cur[:m]/den
    coeffs = [Fr(c, den) for c in cur[:m]]
    h = max(max(abs(x.numerator), x.denominator) for x in coeffs if x)
    print(f"{tag}: reduced representative height {h} (marker {cur[m]}, should be {den})")
    terms = [(nz[j//2], 'a' if j % 2 == 0 else 'b', coeffs[j]) for j in range(m) if coeffs[j]]
    print(f"  {len(terms)} terms:")
    for z, ab, c in terms:
        coef = (c, Fr(0)) if ab == 'a' else (Fr(0), c)
        print(f"    ({G.nf_str(coef)})*Li5({z})")
