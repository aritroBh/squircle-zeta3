from fractions import Fraction as Fr
def solve_lin(cols, T):
    """Solve A x = T exactly (over Q). Returns particular solution or None."""
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
cols = [(Fr(1), Fr(0), Fr(0)), (Fr(1), Fr(1), Fr(0)), (Fr(0), Fr(1), Fr(1))]
sp = solve_lin(cols, (Fr(1), Fr(2), Fr(3)))
print("fixed solve:", sp, "(expect 2, -1, 3)")
acc = [Fr(0)]*3
for j, c in enumerate(sp):
    for r in range(3): acc[r] += c*cols[j][r]
print("check:", acc)
cols2 = [(Fr(1), Fr(2), Fr(3), Fr(4)), (Fr(1), Fr(-1), Fr(1), Fr(0))]
sp2 = solve_lin(cols2, (Fr(3), Fr(0), Fr(5), Fr(4)))
acc = [Fr(0)]*4
for j, c in enumerate(sp2):
    for r in range(4): acc[r] += c*cols2[j][r]
print("check2:", acc, "(target 3,0,5,4)")
