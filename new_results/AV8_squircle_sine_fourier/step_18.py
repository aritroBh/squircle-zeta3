from mpmath import mp, mpf, pi, matrix, lu_solve, factorial
mp.dps = 60
import sys; sys.path.insert(0,'.')
from squircle import Psi1, Psi2, Psi3, Phi_p

phi = mpf('0.9')
h = mpf('0.005')
offs = list(range(-4,5))
A = matrix(9,9); b = matrix(9,1)
for i,j in enumerate(offs):
    for k in range(9):
        A[i,k] = (j*h)**k
    b[i] = Phi_p(__import__('squircle').sin(phi), 2+j*h)
c = lu_solve(A,b)
print("Psi1:", mp.nstr(Psi1(phi),25), " FD:", mp.nstr(c[1]*1,25))
print("Psi2:", mp.nstr(Psi2(phi),25), " FD:", mp.nstr(c[2]*2,25))
print("Psi3:", mp.nstr(Psi3(phi),25), " FD:", mp.nstr(c[3]*6,25))
