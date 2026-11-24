
import sympy as sp, json, os
WORK = os.path.dirname(os.path.abspath(__file__))
def cross(P,Q): return P[0]*Q[1]-P[1]*Q[0]
def edge_data(b):
    V = [sp.Matrix([1,0]), sp.Matrix([2+2*b,2]), sp.Matrix([b,1])]
    kk = [sp.Rational(1,4), sp.Rational(1,8), sp.Rational(1,4)]; ww = [1,3,1]
    verts = V + [-v for v in V]; out = []
    for q in range(3):
        V0, V1 = verts[q], verts[(q+1)%6]
        D = V1-V0; L2 = D.dot(D); L = sp.sqrt(L2); tau = D/L
        h = sp.simplify(abs(cross(V0,tau))); s0 = V0.dot(tau); u = sp.simplify(-s0/L); v = sp.simplify(1-u)
        j0 = [j for j in range(3) if cross(V[j],V0)==0][0]
        j1 = [j for j in range(3) if cross(V[j],V1)==0][0]
        al = sp.simplify(ww[j0]*abs(kk[j0]*cross(V[j0],V1)))
        be = sp.simplify(ww[j1]*abs(kk[j1]*cross(V[j1],V0)))
        other = [j for j in range(3) if j not in (j0,j1)]
        assert len(other)==1, "edge has != 1 non-local ratio"
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
    return out
res = {}
for b in [sp.Rational(0), sp.Rational(1)]:
    eds = edge_data(b)
    res[str(b)] = eds
    print(f"===== b = {b} : 3 edges, each with exactly ONE non-local ratio (asserted) =====")
    for e in eds:
        print(f"edge {e['edge']}: V0={e['V0']} V1={e['V1']}  A0={e['A0']} A1={e['A1']} B={e['B']} "
              f"r={e['r']} kap={e['kap']} L2={e['L2']} q0={e['q0']} V02={e['V02']} Z={e['Z']}")
json.dump(res, open(os.path.join(WORK,"geometry.json"),"w"), indent=1)
print("saved geometry.json")
