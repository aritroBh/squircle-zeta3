import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp
sqrt2 = sp.sqrt(2)
t = sp.symbols('t')
normals = [(sqrt2/2,sqrt2/2),(0,1),(-sqrt2/2,sqrt2/2),(-1,0)]
V0 = sp.Matrix([0, sqrt2-1]); V1 = sp.Matrix([1-sqrt2/2, 1-sqrt2/2])
y = V0 + t*(V1-V0); yp = V1-V0

# zeros of a_j (signed term values) 
svals = [sp.expand(sp.Matrix(n).dot(y)) for n in normals]
print("Zeros of term values a_j(t):")
zeros=[]
for j,s in enumerate(svals):
    z = sp.solve(sp.Eq(s,0), t)
    zeros += z
    print(f"  a_{j+1}(t)=0 at t =", [sp.simplify(zz) for zz in z], "=", [sp.N(zz,6) for zz in z])

# geometric quadratics
yysq = sp.expand(y.dot(y))          # |y(t)|^2
ypyp = sp.expand(yp.dot(yp))        # |y'|^2 (constant)
yyp  = sp.expand(y.dot(yp))         # y.y'
gram = sp.expand(yysq*ypyp - yyp**2)  # Gram det = |y x y'|^2 (should be const in 2D for linear y)
print("\n|y(t)|^2 =", sp.factor(yysq))
print("|y'|^2 =", sp.nsimplify(ypyp), "=", sp.N(ypyp))
print("y.y' =", sp.expand(yyp))
print("Gram =", sp.nsimplify(gram), "=", sp.N(gram))

print("\nRoots of |y(t)|^2 = 0 (complex):")
ry = sp.solve(sp.Eq(yysq,0), t)
for r in ry:
    print("   t =", sp.simplify(r), "=", sp.N(r, 10))
log("step: a_j zeros at t=0,1,-(1+sqrt2),2+sqrt2; computed |y|^2 roots (complex)")
