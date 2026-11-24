import sys, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
Fraction = Fr
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def q12(v):
    out = []
    for q in v: out += [q[0], q[1]]
    return tuple(out)
SQ5v = (Fr(0), Fr(1))
def cols_z(z):
    cl = li5cl[z]
    return [q12(cl), q12(tuple(nfmul(SQ5v, q) for q in cl))]
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
    for i, c in enumerate(piv): part[c] = A[i][m]
    return part
combo = ((Fr(3,4),Fr(-1,4)), (Fr(1,4),Fr(1,4)), (Fr(-2),Fr(1)), (Fr(3),Fr(-1)))
cols = []
for z in combo: cols += cols_z(z)
T = ((Fr(192), Fr(-24)), (Fr(0), Fr(0)), (Fr(-108), Fr(0)), (Fr(0), Fr(0)),
     (Fr(-147), Fr(-3)), (Fr(636), Fr(-69)))
T12 = q12(T)
sp = solve_lin(cols, T12)
print("solve:", [str(x) for x in sp])
# verify
acc = [Fr(0)]*12
for j, c in enumerate(sp):
    if c:
        for r in range(12):
            acc[r] += c * cols[j][r]
print("residual:", [str(x) for x in acc], "== 0?", all(x == 0 for x in acc))
for i, z in enumerate(combo):
    print(f"  ({G.nf_str((sp[2*i], sp[2*i+1]))})*Li5({z})")
