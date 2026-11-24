import pickle, sys
sys.path.insert(0, ".")
from fractions import Fraction as F
import obstruction32 as OB
D = pickle.load(open("dod5d_symbols.pkl","rb"))
Psi = OB.psi_tensor(D['total'])
print("Psi entries:", len(Psi))
# Dynkin projections
def D2(w):
    (a, b) = w
    return {(a, b): F(1,2), (b, a): F(-1,2)}
def D3(w):
    (a, b, c) = w
    out = {}
    for u, cf in [((a,b,c), 1), ((a,c,b), -1), ((b,c,a), -1), ((c,b,a), 1)]:
        out[u] = out.get(u, F(0)) + F(cf,3)
    return {k: v for k, v in out.items() if v}
# check D3 is a projection onto Lie_3 (idempotent on basis elts): quick test on [X0,[X0,X1]] tensor form
t = D3((0,0,1))
t2 = {}
for w, cf in t.items():
    for w2, cf2 in D3(w).items():
        t2[w2] = t2.get(w2, F(0)) + cf*cf2
print("D3 idempotent on its image:", {k: v for k, v in t2.items() if v} == t)
proj = {}
for (u, v), c in Psi.items():
    du = D3(u); dv = D2(v)
    for u2, cu in du.items():
        for v2, cv in dv.items():
            key = (u2, v2)
            proj[key] = proj.get(key, 0) + c*(cu*cv)
diff = {}
for k in set(Psi) | set(proj):
    d = proj.get(k, 0) - Psi.get(k, 0)
    if d != 0: diff[k] = d
print("Psi in Lie_3(x)Lie_2 (projected == itself):", not diff, " (ndiff:", len(diff), ")")
# same check for the octagon total (control)
bs = pickle.load(open("<project>/research/08_theory/OCT5D2_octagon_depth_test/block_symbols.pkl","rb"))
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from qsqrt2 import QS2 as Q
So = {tuple(map(int, w)): Q(F(p), F(q)) for w, (p, q) in bs['total'].items()}
Psio = OB.psi_tensor(So)
projo = {}
for (u, v), c in Psio.items():
    for u2, cu in D3(u).items():
        for v2, cv in D2(v).items():
            key = (u2, v2); projo[key] = projo.get(key, 0) + c*(cu*cv)
diffo = {k: projo.get(k,0)-Psio.get(k,0) for k in set(Psio)|set(projo)}
diffo = {k: v for k, v in diffo.items() if v != 0}
print("octagon Psi in Lie_3(x)Lie_2:", not diffo)
