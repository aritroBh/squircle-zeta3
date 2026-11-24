from fractions import Fraction as Fr
import importlib, json
import cor6_sym
importlib.reload(cor6_sym)
from cor6_sym import *
from cor6_ft import *
import sympy as sp

BASE = "<project>/research/"
D = json.load(open(BASE + "10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
W4TERMS = [(int(k.split(' ', 1)[0][2:]), eval(k.split(' ', 1)[1][1:]), Fr(v))
           for k, v in D if len(eval(k.split(' ', 1)[1][1:])) == 4]

def vecQ(x):
    x = abs(Fr(x)); a = b = 0
    n, d = x.numerator, x.denominator
    while n % 2 == 0: n //= 2; a += 1
    while d % 2 == 0: d //= 2; a -= 1
    while n % 3 == 0: n //= 3; b += 1
    while d % 3 == 0: d //= 3; b -= 1
    assert n == 1 and d == 1, x
    return (Fr(a), Fr(b))

def omega_rho(c0):
    SOm = {}
    for e, w, val in W4TERMS:
        SOm = sadd(SOm, sscale(val * c0**e, sym_H(w, c0, vecQ)))
    SE = expand_sym(SOm, 2)
    T = {}
    for w, c in SE.items():
        iw = tuple(0 if l == (Fr(1), Fr(0)) else 1 for l in w)
        T[iw] = T.get(iw, 0) + c
    return rho(T), delta22(SOm, 2)

def li4_rho(x, co):
    vx, vomx = vecQ(x), vecQ(1 - x)
    S = {(lneg(vomx), vx, vx, vx): co}
    SE = expand_sym(S, 2)
    T = {}
    for w, c in SE.items():
        iw = tuple(0 if l == (Fr(1), Fr(0)) else 1 for l in w)
        T[iw] = T.get(iw, 0) + c
    return T

# (b) paper's 3-Li4 combination at c=1/2
rOm, d22h = omega_rho(Fr(1, 2))
T = {}
for x, co in ((Fr(-2), Fr(-21)), (Fr(1, 3), Fr(13, 2)), (Fr(2, 3), Fr(13, 2))):
    for k, v in li4_rho(x, co).items(): T[k] = T.get(k, 0) + v
r3 = rho(T)
print("c=1/2: rho(Omega) == rho(-21 Li4(-2) + 13/2 Li4(1/3) + 13/2 Li4(2/3)):", rOm == r3)
print("delta22(Omega(1/2)):", d22h)

# (a) c=1/3
rOm3, d223 = omega_rho(Fr(1, 3))
print("\nc=1/3: rho(Omega(1/3)):", {k: str(v) for k, v in sorted(rOm3.items())})
print("delta22(Omega(1/3)):", d223)
# Li4 span test: find combination of Li4 at {2,3}-unit args matching the class
sols23 = [Fr(-8), Fr(-3), Fr(-2), Fr(-1), Fr(-1, 2), Fr(-1, 3), Fr(-1, 8), Fr(1, 9), Fr(1, 4),
          Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(3, 4), Fr(8, 9), Fr(9, 8), Fr(4, 3), Fr(3, 2),
          Fr(2), Fr(3), Fr(4), Fr(9)]
# classes of Li4 at the 21 args span Q^4 (dim 3); find a combination equal to rOm3
Q2 = W4Quotient(2)
def li4_sym(x):
    vx, vomx = vecQ(x), vecQ(1 - x)
    return {(lneg(vomx), vx, vx, vx): Fr(1)}
clsOm3 = Q2.coords({k: v for k, v in expand_sym({}, 2).items()} | {})  # zero
SOm3 = {}
for e, w, val in W4TERMS:
    SOm3 = sadd(SOm3, sscale(val * Fr(1, 3)**e, sym_H(w, Fr(1, 3), vecQ)))
clsOm3 = Q2.coords(SOm3)
print("class(Omega(1/3)) in Q4 coords:", clsOm3)
import itertools
vecs = [Q2.coords(li4_sym(x)) for x in sols23]
# solve for a small combination
from flint import fmpz_mat, fmpz
M = fmpz_mat(3, len(sols23))
for j, v in enumerate(vecs):
    for i in range(3):
        M[i, j] = fmpz(v[i].numerator * 10**0)  # denominators: handle below
# instead do exact rational solve over subsets of size 3
found = None
for combo in itertools.combinations(range(len(sols23)), 3):
    A = [[vecs[j][i] for j in combo] for i in range(3)]
    det = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
    if det == 0: continue
    # solve A e = clsOm3
    import numpy as np
    b = clsOm3
    e = []
    ok = True
    # Cramer
    def det3(M):
        return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])
    for col in range(3):
        Mc = [row[:] for row in A]
        for i in range(3): Mc[i][col] = b[i]
        e.append(det3(Mc)/det)
    found = (combo, e)
    break
if found:
    combo, e = found
    print("Li4 combination with same class:", [(str(sols23[j]), str(c)) for j, c in zip(combo, e)])
