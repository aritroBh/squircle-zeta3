import sys
sys.path.insert(0, '<project>/research/08_theory/COR6_true_point_symbols')
from fractions import Fraction as Fr
import mpmath as mp
from cor6_nf import NF
from cor6_sunit import LV, exponent_vector, log_matrix, num_rank

F = NF([-17,0,1])
basis = [(Fr(2), Fr(0)), (Fr(-4), Fr(-4)), (Fr(-4), Fr(-1))]
x = (Fr(-5), Fr(1))   # -5 + sqrt17 = -(5 - sqrt17)
rts = F.roots(150)
bx = LV(F, x, rts, [2])
B = [LV(F, g, rts, [2]) for g in basis]
print("LV(x) =", [mp.nstr(v, 10) for v in bx])
for g, b in zip(basis, B):
    print(f"LV{g} =", [mp.nstr(v, 10) for v in b])
# least squares by hand
import itertools
Bt = [list(col) for col in zip(*B)]
r = 3
G = [[sum(Bt[k][i]*Bt[k][j] for k in range(len(Bt))) for j in range(r)] for i in range(r)]
bb = [sum(Bt[k][i]*bx[k] for k in range(len(Bt))) for i in range(r)]
A = [row[:] + [bb[i]] for i, row in enumerate(G)]
for col in range(r):
    p = max(range(col, r), key=lambda rr: abs(A[rr][col]))
    A[col], A[p] = A[p], A[col]
    for rr in range(col+1, r):
        A[rr] = [x_ - y_*A[rr][col]/A[col][col] for x_, y_ in zip(A[rr], A[col])]
e = [mp.mpf(0)]*r
for i in range(r-1, -1, -1):
    e[i] = (A[i][r] - sum(A[i][j]*e[j] for j in range(i+1, r)))/A[i][i]
print("raw exponents:", [mp.nstr(v, 30) for v in e])
# exact check: is -5+sqrt17 == +- 2^e0 * (-4-4s17)^e1 * (-4-s17)^e2 for nice e's?
# try: 5-sqrt17 vs basis elements: N(5-s17) = 8. guess (5-s17) = 2^a * (-4-4s17)^b * (-4-s17)^c?
g0 = F.from_int(2); g1 = (Fr(-4), Fr(-4)); g2 = (Fr(-4), Fr(-1))
for a, b, c in itertools.product(range(-4, 5), range(-4, 5), range(-4, 5)):
    v = F.mul(F.pow(g0, a), F.mul(F.pow(g1, b), F.pow(g2, c)))
    if F.eq(v, x) or F.eq(v, F.neg(x)):
        print("EXACT:", (a, b, c))
