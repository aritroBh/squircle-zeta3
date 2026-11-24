import os, sympy as sp
from sympy import Matrix, Rational, sqrt

WORK = "research/08_theory/OCT7_generic_octagon"
os.chdir(WORK)
open("progress.log","a").write("STEP 01: geometry construction for rational centrally symmetric octagon\n")

def build_octagon(normals):
    """normals: list of 4 integer 2-vectors. K = {x: sum |<m_j,x>| <= 1}.
    Returns list of edges; each edge = dict(wedge signs s, V0, V1, terms)."""
    import itertools
    ms = [Matrix([Rational(v[0]), Rational(v[1])]) for v in normals]
    # all sign patterns mod global sign: fix s0 = +1
    pats = []
    for bits in itertools.product([1,-1], repeat=3):
        s = (1,) + bits
        pats.append(s)
    edges = []
    for s in pats:
        M = sum((s[j]*ms[j] for j in range(4)), Matrix([Rational(0),Rational(0)]))
        if M == Matrix([0,0]):  # degenerate wedge
            return None
        # candidate endpoints: intersect <M,x>=1 with <m_j,x>=0
        verts = []
        for j in range(4):
            A = Matrix([[M[0], M[1]],[ms[j][0], ms[j][1]]])
            b = Matrix([1, 0])
            if A.det() == 0: continue
            v = A.solve(b)
            # check wedge sign conditions: s_i * <m_i, v> >= 0 for all i
            ok = True
            for i in range(4):
                val = s[i]*(ms[i].dot(v))
                if val < 0: ok=False; break
            if ok:
                # check it's a vertex of the segment (some other sign strictly... )
                verts.append((j, v))
        if len(verts) != 2:
            return None  # not a proper octagon edge
        (j0, V0), (j1, V1) = verts
        edges.append({'s': s, 'M': M, 'V0': V0, 'V1': V1, 'j0': j0, 'j1': j1})
    return edges

def edge_data(e, normals):
    ms = [Matrix([Rational(v[0]), Rational(v[1])]) for v in normals]
    s, V0, V1 = e['s'], e['V0'], e['V1']
    d = V1 - V0
    terms = []
    for j in range(4):
        a0 = s[j]*ms[j].dot(V0)   # value at t=0
        a1 = s[j]*ms[j].dot(V1)   # value at t=1
        terms.append((a0, a1))
    # identify local (vanish at endpoints) and non-local
    loc0 = [j for j in range(4) if terms[j][0]==0]  # vanish at t=0 -> A0*t
    loc1 = [j for j in range(4) if terms[j][1]==0]  # vanish at t=1 -> A1*(1-t)
    nl   = [j for j in range(4) if terms[j][0]!=0 and terms[j][1]!=0]
    return terms, loc0, loc1, nl

# Try the brief-style family first, then scan
candidates = [
    [(1,0),(0,1),(1,1),(1,2)],
    [(1,0),(0,1),(1,1),(2,1)],
    [(1,0),(0,1),(1,1),(1,3)],
    [(1,0),(0,1),(2,1),(1,3)],
    [(1,0),(0,1),(1,2),(2,3)],
    [(1,0),(0,1),(1,1),(3,1)],
    [(1,0),(0,1),(3,1),(1,2)],
    [(1,0),(0,1),(1,2),(3,4)],
]
for nm in candidates:
    edges = build_octagon(nm)
    print("normals", nm, "->", "DEGENERATE" if edges is None else f"{len(edges)} edges OK")
