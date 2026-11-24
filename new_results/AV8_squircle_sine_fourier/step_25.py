from mpmath import mp, mpf, pi, sin, cos, matrix, lu_solve
mp.dps = 80
import sys; sys.path.insert(0,'.')
from squircle import Psi3, Phi_p
phi = mpf('0.9')
h = mpf('0.003'); offs=list(range(-4,5))
A = matrix(9,9); b2 = matrix(9,1)
for i,j in enumerate(offs):
    for k in range(9): A[i,k]=(j*h)**k
    b2[i] = Phi_p(sin(phi), 2+j*h)
c2 = lu_solve(A,b2)
print("Psi3:", mp.nstr(Psi3(phi),30))
print("FD  :", mp.nstr(c2[3]*6,30))
print("diff:", mp.nstr(Psi3(phi)-c2[3]*6,5))
