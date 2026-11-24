import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
from cor6_sym import rref

li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
def q12(v):
    out = []
    for q in v: out += [q[0], q[1]]
    return tuple(out)
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
SQ5v = (Fr(0), Fr(1))
zs = list(li5cl.keys())
cols, names = [], []
for z in zs:
    cl = li5cl[z]
    cols.append(q12(cl)); names.append((z, 0))
    cols.append(q12(tuple(nfmul(SQ5v, q) for q in cl))); names.append((z, 1))

# target classes (from earlier computation)
targets = {
 'gol1': ((Fr(84), Fr(-3)), (Fr(0), Fr(0)), (Fr(60), Fr(21)), (Fr(0), Fr(0)),
          (Fr(513), Fr(-63,2)), (Fr(1689), Fr(-129))),
 'gol3': ((Fr(-6), Fr(0)), (Fr(0), Fr(0)), (Fr(-21), Fr(3)), (Fr(0), Fr(0)),
          (Fr(-84), Fr(18)), (Fr(-357), Fr(87))),
}

def solve_lin(cols, T):
    """Solve A x = T exactly; returns (particular, kernel_basis, pivots) or None."""
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
        if A[r][m]: return None  # inconsistent
    free = [c for c in range(m) if c not in piv]
    part = [Fr(0)] * m
    for i, c in enumerate(piv):
        part[c] = A[i][m]
    ker = []
    for f in free:
        v = [Fr(0)] * m; v[f] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -A[i][f]
        ker.append(v)
    return part, ker, piv

for tag, T in targets.items():
    res = solve_lin(cols, q12(T))
    if res is None:
        print(f"{tag}: class(Omega5 - N5) NOT in Li5 span!")
        continue
    part, ker, piv = res
    print(f"{tag}: class in Li5 span: YES. kernel dim {len(ker)}")
    terms = []
    for (z, ab), c in zip(names, part):
        if c:
            terms.append((z, ab, c))
    print(f"  particular representative P5 ({len(terms)} nonzero coeffs, free=0):")
    for z, ab, c in terms:
        coef = (c, Fr(0)) if ab == 0 else (Fr(0), c)
        print(f"    ({G.nf_str(coef)})*Li5({z})   ~ {float(mp.mpf(z[0].numerator)/z[0].denominator + mp.mpf(z[1].numerator)/z[1].denominator*mp.sqrt(5)):+.6f}")
import mpmath as mp
mp.mp.dps = 20
