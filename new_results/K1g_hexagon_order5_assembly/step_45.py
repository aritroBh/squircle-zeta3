import sympy as sp, mpmath as mp
mp.mp.dps = 60
# Golden hexagon (from planner_check_irrational.py): V=(1,0),(phi,1),(0,1), w=(2,3,5)
phi = (1+sp.sqrt(5))/2
V = [sp.Matrix([1,0]), sp.Matrix([phi,1]), sp.Matrix([0,1])]
w = [2,3,5]
cross = lambda P,Q: P[0]*Q[1]-P[1]*Q[0]
c = sp.symbols('c0:3')
sol = sp.solve([sp.Eq(sum(c[j]*sp.Abs(cross(V[j],V[q])) for j in range(3)),1) for q in range(3)], c, dict=True)[0]
k = [sp.simplify(sol[c[j]]/w[j]) for j in range(3)]
print("k =", [sp.nsimplify(kj) for kj in k])
print("c =", [sp.simplify(sol[c[j]]) for j in range(3)])
# Edge data for 3 edges: verts[0]->verts[1], [1]->[2], [2]->[3]=-V0
verts = [V[0],V[1],V[2],-V[0]]
edges=[]
for e in range(3):
    V0,V1 = verts[e],verts[e+1]
    j0=[j for j in range(3) if sp.simplify(cross(V[j],V0))==0][0]
    j1=[j for j in range(3) if sp.simplify(cross(V[j],V1))==0][0]
    jo=3-j0-j1
    A0=sp.simplify(w[j0]*sp.Abs(k[j0]*cross(V[j0],V1)))
    A1=sp.simplify(w[j1]*sp.Abs(k[j1]*cross(V[j1],V0)))
    B=sp.simplify(w[jo]*sp.Abs(k[jo]*cross(V[jo],V0)))
    r=sp.simplify(cross(V[jo],V1)/cross(V[jo],V0))
    L2=sp.simplify((V1-V0).dot(V1-V0)); q0=sp.simplify(V0.dot(V1-V0)); V02=sp.simplify(V0.dot(V0))
    assert sp.simplify(A1+B-1)==0 and sp.simplify(A0-(1-B*r))==0
    edges.append(dict(e=e+1,V0=V0,V1=V1,j0=j0,j1=j1,jo=jo,A0=A0,A1=A1,B=B,r=r,kap=sp.simplify(r-1),
                      L2=L2,q0=q0,V02=V02,w0=sp.Integer(w[j0]),w1=sp.Integer(w[j1]),wo=sp.Integer(w[jo])))
    print(f"edge {e+1}: A0={sp.nsimplify(A0)}, A1={sp.nsimplify(A1)}, B={sp.nsimplify(B)}, r={sp.nsimplify(r)}, kap={sp.nsimplify(sp.simplify(r-1))}, L2={L2}, q0={q0}, V02={V02}, w=({w[j0]},{w[j1]},{w[jo]})")
import pickle
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/golden_edges.pkl','wb') as f: pickle.dump({'edges':edges,'k':k,'V':V,'w':w},f)
print("saved golden_edges.pkl")
