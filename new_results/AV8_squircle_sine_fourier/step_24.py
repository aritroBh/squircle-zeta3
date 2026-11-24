from mpmath import mp, mpf, pi, sin, cos, matrix, lu_solve
mp.dps = 60
import sys; sys.path.insert(0,'.')
from squircle import H1u, H2u, H3u, Psi1, Psi2, Psi3, Phi_p

u = mpf('0.9'); y = sin(u); cu = cos(u)
h = mpf('0.005'); offs=list(range(-4,5))
A = matrix(9,9); b = matrix(9,1)
for i,j in enumerate(offs):
    for k in range(9): A[i,k]=(j*h)**k
    b[i] = (1-y**(2+j*h))**(-1/(2+j*h))
c = lu_solve(A,b)
print("H1:", mp.nstr(H1u(u)/cu - c[1],6))
print("H2:", mp.nstr(H2u(u)/cu - c[2]*2,6))
print("H3:", mp.nstr(H3u(u)/cu - c[3]*6,6))
# now Psi3
phi = u
b2 = matrix(9,1)
for i,j in enumerate(offs):
    b2[i] = Phi_p(sin(phi), 2+j*h)
c2 = lu_solve(A,b2)
print("Psi1:", mp.nstr(Psi1(phi)-c2[1],6))
print("Psi2:", mp.nstr(Psi2(phi)-c2[2]*2,6))
print("Psi3:", mp.nstr(Psi3(phi)-c2[3]*6,6))
