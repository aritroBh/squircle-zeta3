import sympy as sp, itertools
from sympy import Matrix, Rational

nm = [(1,0),(0,1),(1,1),(1,2)]
ms = [Matrix([Rational(v[0]), Rational(v[1])]) for v in nm]
for bits in itertools.product([1,-1], repeat=3):
    s = (1,) + bits
    M = sum((s[j]*ms[j] for j in range(4)), Matrix([Rational(0),Rational(0)]))
    if M == Matrix([0,0]): 
        print(s, "M=0"); continue
    verts = []
    for j in range(4):
        A = Matrix([[M[0], M[1]],[ms[j][0], ms[j][1]]])
        if A.det() == 0: continue
        v = A.solve(Matrix([1, 0]))
        checks = [s[i]*(ms[i].dot(v)) for i in range(4)]
        if all(c >= 0 for c in checks):
            verts.append((j, v))
    print(s, "M=", M.T, "verts:", [(j, v.T) for j,v in verts])
