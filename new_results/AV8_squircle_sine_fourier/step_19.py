from mpmath import mp, mpf, pi, sin, cos, matrix, lu_solve
mp.dps = 60
import sys; sys.path.insert(0,'.')
from squircle import H1u, H2u, H3u

u = mpf('0.9'); y = sin(u); cu = cos(u)
h = mpf('0.005'); offs=list(range(-4,5))
A = matrix(9,9); b = matrix(9,1)
for i,j in enumerate(offs):
    for k in range(9): A[i,k]=(j*h)**k
    b[i] = (1-y**(2+j*h))**(-1/(2+j*h))
c = lu_solve(A,b)
print("H1:", mp.nstr(H1u(u)/cu,25), mp.nstr(c[1],25))
print("H2:", mp.nstr(H2u(u)/cu,25), mp.nstr(c[2]*2,25))
print("H3:", mp.nstr(H3u(u)/cu,25), mp.nstr(c[3]*6,25))
