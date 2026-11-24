import mpmath
from mpmath import mp, mpf, sqrt, pi, zeta, matrix, lu_solve
mp.dps = 75
# test matrix construction
Q = [mpf(1)+i for i in range(6)]
b = matrix(6,1)
for i in range(6): b[i] = Q[i]
A = matrix(6,6)
for i in range(6):
    for j in range(6):
        A[i,j] = mpf(i+1)**j
c = lu_solve(A,b)
print("test solve (expect [1,1,0,0,0,0]):", [mp.nstr(ci,6) for ci in c])
# residual
print("resid:", mp.nstr((A*c-b).norm(),4))
