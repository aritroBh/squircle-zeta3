import os, sympy as sp, itertools
from sympy import Matrix, Rational
WORK = "research/08_theory/OCT7_generic_octagon"
os.chdir(WORK)

def build_octagon(normals):
    ms = [Matrix([Rational(v[0]), Rational(v[1])]) for v in normals]
    edges = []
    for bits in itertools.product([1,-1], repeat=3):
        s = (1,) + bits
        M = sum((s[j]*ms[j] for j in range(4)), Matrix([Rational(0),Rational(0)]))
        if M == Matrix([0,0]): continue
        verts = []
        for j in range(4):
            A = Matrix([[M[0], M[1]],[ms[j][0], ms[j][1]]])
            if A.det() == 0: continue
            v = A.solve(Matrix([1, 0]))
            if all(s[i]*(ms[i].dot(v)) >= 0 for i in range(4)):
                if all((v - w).norm() != 0 for _, w in verts):
                    verts.append((j, v))
        if len(verts) == 2:
            edges.append({'s': s, 'V0': verts[0][1], 'V1': verts[1][1],
                          'j0': verts[0][0], 'j1': verts[1][0]})
        elif len(verts) > 2:
            return None  # degenerate
    return edges if len(edges)==8 else None

def edge_terms(e, normals):
    ms = [Matrix([Rational(v[0]), Rational(v[1])]) for v in normals]
    s, V0, V1 = e['s'], e['V0'], e['V1']
    terms = [(s[j]*ms[j].dot(V0), s[j]*ms[j].dot(V1)) for j in range(4)]
    return terms

def ratios(terms):
    nl = [(a0,a1) for (a0,a1) in terms if a0!=0 and a1!=0]
    rs = [a1/a0 for a0,a1 in nl]
    return rs

candidates = [
    [(1,0),(0,1),(1,1),(1,2)],
    [(1,0),(0,1),(1,1),(2,1)],
    [(1,0),(0,1),(1,1),(1,3)],
    [(1,0),(0,1),(2,1),(1,3)],
    [(1,0),(0,1),(1,2),(2,3)],
    [(1,0),(0,1),(1,1),(3,1)],
    [(1,0),(0,1),(3,1),(1,2)],
    [(1,0),(0,1),(1,2),(3,4)],
    [(2,1),(1,2),(1,0),(0,1)],
    [(1,0),(1,2),(0,1),(2,1)],
]
for nm in candidates:
    edges = build_octagon(nm)
    if edges is None:
        print("normals", nm, "-> DEGENERATE"); continue
    print("normals", nm, "-> 8 edges OK")
    seen = set()
    for e in edges[:4]:  # half boundary
        terms = edge_terms(e, nm)
        rs = ratios(terms)
        print("   edge", e['V0'].T, "->", e['V1'].T, "terms:", terms, "ratios:", [str(r) for r in rs])
