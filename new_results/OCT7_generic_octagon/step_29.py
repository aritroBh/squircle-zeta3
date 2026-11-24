import sympy as sp, itertools
from sympy import Matrix, Rational

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
                if not any((v-w).norm()==0 for _,w in verts):
                    verts.append((j,v))
        if len(verts) == 2:
            edges.append({'s': s, 'V0': verts[0][1], 'V1': verts[1][1],
                          'j0': verts[0][0], 'j1': verts[1][0], 'M': M})
        elif len(verts) > 2:
            return None
    return edges if len(edges)==4 else None

def edge_terms(e, normals):
    ms = [Matrix([Rational(v[0]), Rational(v[1])]) for v in normals]
    s, V0, V1 = e['s'], e['V0'], e['V1']
    return [(s[j]*ms[j].dot(V0), s[j]*ms[j].dot(V1)) for j in range(4)]

def analyze(nm):
    edges = build_octagon(nm)
    if edges is None: return None
    out = []
    for e in edges:
        terms = edge_terms(e, nm)
        nl = [(j,a0,a1) for j,(a0,a1) in enumerate(terms) if a0!=0 and a1!=0]
        rs = [a1/a0 for j,a0,a1 in nl]
        out.append((e, terms, nl, rs))
    return out

def props(rs):
    if len(rs)!=2: return "not 2 nonlocal"
    r1,r2 = rs
    x1,x2 = 1-r1, 1-r2
    flags = []
    if r1==1 or r2==1: flags.append("r=1")
    if r1==r2: flags.append("r1=r2")
    if r1*r2==1: flags.append("LANDEN(r1r2=1)")
    if x2==x1/(x1-1): flags.append("x-landen")
    if x2==1/x1: flags.append("INVERSION")
    if x2==x1: flags.append("x1=x2")
    return ",".join(flags) if flags else "GENERIC"

candidates = [
    [(1,0),(0,1),(1,1),(1,2)],
    [(1,0),(0,1),(1,1),(2,1)],
    [(1,0),(0,1),(1,1),(1,3)],
    [(1,0),(0,1),(2,1),(1,3)],
    [(1,0),(0,1),(1,2),(2,3)],
    [(1,0),(0,1),(1,1),(3,1)],
    [(1,0),(0,1),(3,1),(1,2)],
    [(1,0),(0,1),(1,2),(3,4)],
    [(1,0),(0,1),(1,3),(2,5)],
    [(1,0),(0,1),(2,3),(3,5)],
    [(1,0),(0,1),(1,2),(2,5)],
    [(1,0),(0,1),(1,4),(2,7)],
]
for nm in candidates:
    res = analyze(nm)
    if res is None:
        print("normals", nm, "DEGEN"); continue
    print("normals", nm)
    for e, terms, nl, rs in res:
        print("   ", e['V0'].T, "->", e['V1'].T, "r:", [str(r) for r in rs], "|", props(rs))
