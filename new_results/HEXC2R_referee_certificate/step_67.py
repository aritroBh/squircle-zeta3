
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 40
from fractions import Fraction as Fr
import reflib as R

def qadd(x, y): return (x[0] + y[0], x[1] + y[1])
def qmul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c + b*d)
def qdiv(x, y):
    a, b = y
    n = a*a + a*b + b*b
    return qmul(x, (Fr(a + b, n), Fr(-b, n)))
def qpow(x, n):
    out = (Fr(1), Fr(0))
    for _ in range(n): out = qmul(out, x)
    return out
def qsub(x, y): return (x[0]-y[0], x[1]-y[1])
ZETA = (Fr(0), Fr(1))
MU6 = [qpow(ZETA, k) for k in range(6)]

def mobius_from_3(f0, f1, finf):
    A = finf; C = (Fr(1), Fr(0))
    D = qdiv(qsub(finf, f1), qsub(f1, f0))
    B = qmul(f0, D)
    return (A, B, C, D)
def inv_map(M):
    A, B, C, D = M
    return (D, (-B[0], -B[1]), (-C[0], -C[1]), A)
def ev(x):
    return mp.mpf(str(x[0])) + mp.mpc(str(x[1])) * (mp.e**(mp.pi*mp.j/3))

E = R.MAPS6[0]
f0, f1, finf = MU6[E[0]], MU6[E[2]], MU6[E[4]]
M = mobius_from_3(f0, f1, finf)
# sanity: does M give the right images?
def qmob(M, x):
    A, B, C, D = M
    n = qadd(qmul(A, (x, Fr(0))), B); d = qadd(qmul(C, (x, Fr(0))), D)
    return qdiv(n, d)
print("M(0) =", ev(qmob(M, Fr(0))), "expect 1")
print("M(1) =", ev(qmob(M, Fr(1))), "expect zeta^2", mp.e**(2*mp.pi*mp.j/3))
print("M(1/2) =", ev(qmob(M, Fr(1,2))), "expect zeta", mp.e**(mp.pi*mp.j/3))
A, B, C, D = inv_map(M)
g = lambda u: (ev(A)*u + ev(B))/(ev(C)*u + ev(D))
print("g(f(1/2)) =", mp.nstr(g(mp.e**(mp.pi*mp.j/3)), 10))
# sig=1 direct
A2 = qsub(A, C); B2 = qsub(B, D)
Rnum = [qmul((Fr(-1), Fr(0)), qsub(qmul(A2, D), qmul(B2, C)))]
Rden = [B2, A2], [D, C]
def poly_ev(P, u):
    out = mp.mpc(0)
    for a in reversed(P): out = out * u + ev(a)
    return out
u = mp.mpf(1)/3
print("R(u) =", mp.nstr(poly_ev(Rnum, u) / (poly_ev(Rden[0], u) * poly_ev(Rden[1], u)), 15))
print("expect -1.19911209754768j")
