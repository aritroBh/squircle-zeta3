code = r'''
import os, pickle, itertools
import sympy as sp
from sympy import Matrix, Rational

WORK = "research/08_theory/OCT7_generic_octagon"
os.chdir(WORK)

NORMALS = [(1,0),(0,1),(2,1),(1,3)]
ms = [Matrix([Rational(v[0]), Rational(v[1])]) for v in NORMALS]

def build_half_boundary(ms):
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
        if len(verts)==2:
            edges.append({'s': s, 'V0': verts[0][1], 'V1': verts[1][1], 'j0': verts[0][0], 'j1': verts[1][0], 'M': M})
        elif len(verts)>2:
            raise ValueError("degenerate")
    return edges

edges = build_half_boundary(ms)
assert len(edges)==4
# order edges into a chain (V1 of one = V0 of next, or j-chain)
def canon(V): return (V[0], V[1])
# sort by angle of midpoint
import math
def ang(e):
    m = (e['V0']+e['V1'])/2
    return math.atan2(float(m[1]), float(m[0]))
edges.sort(key=ang, reverse=True)
print("OCTAGON K = { |x| + |y| + |2x+y| + |x+3y| <= 1 }")
print("Half boundary (4 edges, CCW):")
all_edge_data = []
for idx, e in enumerate(edges):
    s, V0, V1 = e['s'], e['V0'], e['V1']
    terms = []
    for j in range(4):
        a0 = s[j]*ms[j].dot(V0); a1 = s[j]*ms[j].dot(V1)
        terms.append((a0,a1))
    # local terms
    jA0 = [j for j in range(4) if terms[j][0]==0][0]  # vanishes at t=0
    jA1 = [j for j in range(4) if terms[j][1]==0][0]  # vanishes at t=1
    nl  = [j for j in range(4) if terms[j][0]!=0 and terms[j][1]!=0]
    A0 = terms[jA0][1]           # a(t) = A0 * t
    A1 = terms[jA1][0]           # a(t) = A1 * (1-t)
    jB1, jB2 = nl[0], nl[1]
    B1, r1 = terms[jB1][0], terms[jB1][1]/terms[jB1][0]
    B2, r2 = terms[jB2][0], terms[jB2][1]/terms[jB2][0]
    kap1, kap2 = r1-1, r2-1
    # order so that kap1 > kap2
    if kap2 > kap1:
        (jB1,B1,r1,kap1),(jB2,B2,r2,kap2) = (jB2,B2,r2,kap2),(jB1,B1,r1,kap1)
    # term values as linear functions
    a = [None]*4
    t = sp.symbols('t')
    a[jA0] = A0*t
    a[jA1] = A1*(1-t)
    a[jB1] = B1*(1+kap1*t)
    a[jB2] = B2*(1+kap2*t)
    # checks: sum a_j = 1 identically; sign pattern consistency
    assert sp.expand(sum(a)-1)==0, ("sum rule fails", e)
    # interior positivity
    for j in range(4):
        assert sp.simplify(a[j].subs(t, Rational(1,2))) > 0
    d = V1-V0
    L2 = d.dot(d)
    print(f"\nedge {idx}: V0={V0.T} V1={V1.T}  |y'|^2={L2}  |y'|={sp.sqrt(L2)}")
    print(f"  A0={A0} (j{jA0})  A1={A1} (j{jA1})   B1={B1} kap1={kap1} r1={r1} (j{jB1})   B2={B2} kap2={kap2} r2={r2} (j{jB2})")
    x1, x2 = -kap1, -kap2
    print(f"  x1={x1} x2={x2}  beta1={kap1/(1+kap1)} beta2={kap2/(1+kap2)}   r1*r2={r1*r2}  Landen? r1*r2==1: {r1*r2==1}; x2==x1/(x1-1): {x2==sp.simplify(x1/(x1-1))}; x2==1/x1: {x2==1/x1}")
    all_edge_data.append({'V0':V0,'V1':V1,'A0':A0,'A1':A1,'B1':B1,'kap1':kap1,'B2':B2,'kap2':kap2,
                          'r1':r1,'r2':r2,'L2':L2,'terms':a,'signs':e['s'],'jA0':jA0,'jA1':jA1,'jB1':jB1,'jB2':jB2})

# convexity check: consecutive edge directions cross > 0 (CCW) & 8 distinct vertices
verts = [e['V0'] for e in edges]
print("\nvertices (half):", [v.T for v in verts])
crosses = []
for i in range(4):
    d1 = edges[i]['V1']-edges[i]['V0']
    d2 = edges[(i+1)%4]['V1']-edges[(i+1)%4]['V0'] if i<3 else None
    if d2 is not None:
        cr = d1[0]*d2[1]-d1[1]*d2[0]
        crosses.append(cr)
print("consecutive cross products:", crosses, " (all same sign => convex chain)")

pickle.dump({'normals': NORMALS, 'edges': all_edge_data}, open('geometry.pkl','wb'))
with open("progress.log","a") as f:
    f.write("STEP 01: built octagon |x|+|y|+|2x+y|+|x+3y|<=1; 4 half-edges; edge data exact; working edge = r=(4/3,10/9)\n")
print("\nsaved geometry.pkl")
'''
with open("step_01.py","w") as f: f.write(code)
exec(code)
