import sympy as sp
from fractions import Fraction as Fr
# debug the type returned by gauss_jordan_solve
import os, sys, itertools, pickle
base = "<project>"
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT8_symbol_test"))
from oct8 import Quotient, S_Li3, vec_of_rational as vf
Q = Quotient((2,3,5))
g = Q.coord(S_Li3(Fr(1,3), vf))
M = sp.Matrix.hstack(g.T, *[d.T for d in Q.D])
t = g.T
sol = M.gauss_jordan_solve(t)[0]
x = sol[0,0]
print(type(x), repr(x))
print(sp.nsimplify(x) if not isinstance(x, sp.Rational) else "rational")
