import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
import mpmath as mp
mp.mp.dps = 30
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
def solve_lin_dbg(cols, T):
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
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[R])]
        piv.append(c); R += 1
        if R == n: break
    for r in range(R, n):
        if A[r][m]: return None
    part = [Fr(0)] * m
    for i, c in enumerate(piv): part[c] = A[i][m]
    return part
combo = ((Fr(1,2),Fr(-1,4)), (Fr(3,4),Fr(-1,4)), (Fr(1,4),Fr(-1,4)), (Fr(1,4),Fr(1,4)))
cols = []
for z in combo: cols += cols_z(z)
T = ((Fr(192), Fr(-24)), (Fr(0), Fr(0)), (Fr(-108), Fr(0)), (Fr(0), Fr(0)),
     (Fr(-147), Fr(-3)), (Fr(636), Fr(-69)))
T12 = q12(T)
sp = solve_lin_dbg(cols, T12)
print("sp:", [str(x) for x in sp])
acc = [Fr(0)]*12
for j, c in enumerate(sp):
    if c:
        for rr in range(12): acc[rr] += c*cols[j][rr]
print("cols*sp:", [str(x) for x in acc])
print("T12   :", [str(x) for x in T12])
print("residual:", [str(a-b) for a, b in zip(acc, T12)])
