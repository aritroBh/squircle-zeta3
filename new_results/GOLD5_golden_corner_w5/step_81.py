import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
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
def solve_lin(cols, T):
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
    part = [Fr(0)] * m
    for i, c in enumerate(piv):
        part[c] = A[i][m]
    return part
def height(sol):
    return max(max(abs(x.numerator), x.denominator) for x in sol if x)
best = {}
for tag, T in targets.items():
    T12 = q12(T)
    cands = []
    for r in (5, 6):
        for combo in itertools.combinations(nz, r):
            cols = []
            for z in combo: cols += cols_z(z)
            s = solve_lin(cols, T12)
            if s is None: continue
            cands.append((height(s), combo, s))
    cands.sort(key=lambda t: t[0])
    best[tag] = cands
    print(f"==== {tag}: best of 5/6-arg:")
    for h, combo, s in cands[:3]:
        print(f"  h={h}:")
        for i, z in enumerate(combo):
            d = (s[2*i], s[2*i+1])
            if any(d):
                print(f"    ({G.nf_str(d)})*Li5({z})")
import pickle
pickle.dump({'cands_%s' % t: [(str(h), [str(z) for z in combo], [str(x) for x in s]) for h, combo, s in cands[:20]] for t, cands in best.items()},
            open(GOLD + 'p5_cands.pkl', 'wb'))
