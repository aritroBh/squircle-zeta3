
import sympy as sp, json, time
WORK='research/08_theory/K1e_hexagon_closed'
def cross(P,Q): return P[0]*Q[1]-P[1]*Q[0]
V = [sp.Matrix([1,0]), sp.Matrix([3,3]), sp.Matrix([0,1])]
kk = [sp.Rational(1,18)]*3; ww = [sp.Rational(3), sp.Rational(5), sp.Rational(3)]
# norm check
Nf = lambda x: sum(ww[j]*abs(kk[j]*cross(V[j],x)) for j in range(3))
print("N(V_q):", [Nf(Vq) for Vq in V])
# convexity
verts = V + [-v for v in V]
for q in range(6):
    e1 = verts[(q+1)%6]-verts[q]; e2 = verts[(q+2)%6]-verts[(q+1)%6]
    assert cross(e1,e2) > 0, ("non-convex at", q)
print("convex: OK (all consecutive edge crosses = 5)")
out=[]
for q in range(3):
    V0, V1 = verts[q], verts[(q+1)%6]
    D = V1-V0; L2 = D.dot(D); L = sp.sqrt(L2); tau = D/L
    h = sp.simplify(abs(cross(V0,tau))); s0 = V0.dot(tau); u = sp.simplify(-s0/L); v = sp.simplify(1-u)
    j0 = [j for j in range(3) if cross(V[j],V0)==0][0]
    j1 = [j for j in range(3) if cross(V[j],V1)==0][0]
    al = sp.simplify(ww[j0]*abs(kk[j0]*cross(V[j0],V1)))
    be = sp.simplify(ww[j1]*abs(kk[j1]*cross(V[j1],V0)))
    other = [j for j in range(3) if j not in (j0,j1)]
    assert len(other)==1
    jo = other[0]
    r = sp.simplify(cross(V[jo],V1)/cross(V[jo],V0))
    B = sp.simplify(ww[jo]*abs(kk[jo]*cross(V[jo],V0)))
    Z = sp.simplify(6*al*be*h**2/L*(3*(al*u+be*v)-2))
    q0 = V0.dot(D); V02 = V0.dot(V0); kap = sp.simplify(r-1)
    assert sp.simplify(al-(be-B*kap))==0 and sp.simplify(be+B-1)==0
    assert r > 0
    out.append(dict(edge=q+1, A0=str(al), A1=str(be), B=str(B), r=str(r), kap=str(kap),
                    L2=str(L2), q0=str(q0), V02=str(V02), Z=str(Z),
                    w0=str(ww[j0]), w1=str(ww[j1]), wo=str(ww[jo]),
                    h=str(h), u=str(u), v=str(v), V0=[str(c) for c in V0], V1=[str(c) for c in V1]))
    print(f"edge {q+1}: V0={list(V0)} V1={list(V1)} A0={al} A1={be} B={B} r={r} kap={kap} L2={L2} q0={q0} V02={V02}")
    print(f"         h={h} u={u} v={v} Z={Z} = {sp.N(Z,15)}")
json.dump(out, open(f"{WORK}/geometry_fresh.json","w"), indent=1)
print("saved geometry_fresh.json")
