import sys, json
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *

sols = json.load(open("sunit23_solutions.json"))
def vec23(fac):
    s, a, b = fac
    return (Fr(a), Fr(b))
# recompute the li4 coords with a fresh quotient
Q = W4Quotient(2)
li4 = []
for xs, fx, fy in sols:
    x = Fr(xs)
    v = Q.coords(li_sym(4, vec23(fx), vec23(fy)))
    li4.append((x, v))
M = [list(v) for _, v in li4]
Dm, piv = rref(M)
print("rank Li4:", len(Dm), "pivots:", piv)
vF = (Fr(5,2), Fr(1,2), Fr(0))
aug = rref([list(v) for _, v in li4] + [list(vF)])
print("rank augmented:", len(aug))
# solve directly: find coefficients of vF in the 21 vectors using rational linear algebra
# use sympy
import sympy as sp
A = sp.Matrix([[sp.Rational(v.numerator, v.denominator) for v in vv] for _, vv in li4]).T
b = sp.Matrix([sp.Rational(v.numerator, v.denominator) for v in vF])
print("A shape:", A.shape, "rank A:", A.rank(), "rank [A|b]:", A.row_join(b).rank())
sol = A.gauss_jordan_solve(b) if A.rank() == A.row_join(b).rank() else None
print("solution exists:", sol is not None)
