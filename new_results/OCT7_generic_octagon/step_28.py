import sympy as sp
from sympy import Matrix, Rational
nm = [(1,0),(0,1),(1,1),(1,2)]
ms = [Matrix([Rational(v[0]), Rational(v[1])]) for v in nm]
s = (1,1,-1,-1)
M = sum((s[j]*ms[j] for j in range(4)), Matrix([Rational(0),Rational(0)]))
print("M:", M.T)
for j in range(4):
    A = Matrix([[M[0], M[1]],[ms[j][0], ms[j][1]]])
    print(f"j={j}: det={A.det()}", end=" ")
    if A.det() != 0:
        v = A.solve(Matrix([1, 0]))
        checks = [s[i]*(ms[i].dot(v)) for i in range(4)]
        print("v=", v.T, "checks:", checks, "all>=0:", all(c>=0 for c in checks))
    else:
        print()
