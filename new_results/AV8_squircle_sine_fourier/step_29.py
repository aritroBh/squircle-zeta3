from mpmath import mp, mpf, pi, sin, cos, matrix, lu_solve
mp.dps = 60
import sys; sys.path.insert(0,'.')
from squircle import H1u, H2u, H3u, Psi3, Phi_p
u = mpf('0.9'); y=sin(u); cu=cos(u)
h = mpf('0.005'); offs=list(range(-4,5))
A = matrix(9,9); b = matrix(9,1); b2 = matrix(9,1)
for i,j in enumerate(offs):
    for k in range(9): A[i,k]=(j*h)**k
    b[i] = (1-y**(2+j*h))**(-1/(2+j*h))
    b2[i] = Phi_p(y, 2+j*h)
c = lu_solve(A,b); c2 = lu_solve(A,b2)
print("H1err:", mp.nstr(H1u(u)/cu - c[1],5), "H2err:", mp.nstr(H2u(u)/cu - 2*c[2],5), "H3err:", mp.nstr(H3u(u)/cu - 6*c[3],5))
print("Psi3err:", mp.nstr(Psi3(u)-6*c2[3],5))
