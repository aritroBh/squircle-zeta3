import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp
sqrt2 = sp.sqrt(2)
t = sp.symbols('t', real=True)
normals = [(sqrt2/2,sqrt2/2),(0,1),(-sqrt2/2,sqrt2/2),(-1,0)]
V0 = sp.Matrix([0, sqrt2-1]); V1 = sp.Matrix([1-sqrt2/2, 1-sqrt2/2])
y = V0 + t*(V1-V0)
print("y(t) =", [sp.expand(c) for c in y])
# term values (signed) s_j(t) = <n_j, y(t)>
print("\nSigned term values s_j(t) = <n_j, y(t)>:")
for j,n in enumerate(normals):
    s = sp.expand(sp.Matrix(n).dot(y))
    print(f"  s_{j+1}(t) =", sp.factor(s), " = ", sp.nsimplify(s))
# Determine sign on (0,1): evaluate at t=1/2
print("\nSigns at t=1/2:")
for j,n in enumerate(normals):
    s = sp.Matrix(n).dot(y).subs(t, sp.Rational(1,2))
    print(f"  s_{j+1}(1/2) =", sp.nsimplify(s), "=", sp.N(s))
