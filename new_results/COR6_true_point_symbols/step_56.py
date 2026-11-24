from fractions import Fraction as Fr
import importlib
import cor6_nf, cor6_sunit
importlib.reload(cor6_nf); importlib.reload(cor6_sunit)
from cor6_nf import NF
from cor6_sunit import *
import mpmath as mp

F = NF([-3, 0, 1])
S0 = [2, 3]
basis = [(Fr(2), Fr(0)), (Fr(3), Fr(0)), (Fr(-3), Fr(-3))]
c = (Fr(0), Fr(1, 2))
dps = 120
mp.mp.dps = dps
rts = F.roots(dps)
bx = LV(F, c, rts, S0)
B = [LV(F, g, rts, S0) for g in basis]
Bt = [list(col) for col in zip(*B)]
r = len(basis)
G = [[sum(Bt[k][i] * Bt[k][j] for k in range(len(Bt))) for j in range(r)] for i in range(r)]
b = [sum(Bt[k][i] * bx[k] for k in range(len(Bt))) for i in range(r)]
A = [row[:] + [b[i]] for i, row in enumerate(G)]
for col in range(r):
    p = max(range(col, r), key=lambda rr: abs(A[rr][col]))
    A[col], A[p] = A[p], A[col]
    for rr in range(col + 1, r):
        A[rr] = [x - y * A[rr][col] / A[col][col] for x, y in zip(A[rr], A[col])]
e = [mp.mpf(0)] * r
for i in range(r - 1, -1, -1):
    e[i] = (A[i][r] - sum(A[i][j] * e[j] for j in range(i + 1, r))) / A[i][i]
print("least squares e:", [mp.nstr(x, 25) for x in e])
resid = [bx[k] - sum(Bt[k][i] * e[i] for i in range(r)) for k in range(len(Bt))]
print("residual:", [mp.nstr(x, 5) for x in resid])
