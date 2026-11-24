import pickle, sys
sys.path.insert(0, ".")
from fractions import Fraction as F
import obstruction32 as OB
D = pickle.load(open("dod5d_symbols.pkl","rb"))
def D2(w):
    (a, b) = w
    return {(a, b): F(1,2), (b, a): F(-1,2)}
def D3(w):
    (a, b, c) = w
    out = {}
    for u, cf in [((a,b,c), 1), ((a,c,b), -1), ((b,c,a), -1), ((c,b,a), 1)]:
        out[u] = out.get(u, F(0)) + F(cf,3)
    return {k: v for k, v in out.items() if v}
L3 = OB.lie3_basis(3); L2 = OB.lie2_basis(3)
# pairing via Omega = (rho3 x rho2)(Psi) using D3,D2 as tensors and rho = D/n
Psi = OB.psi_tensor(D['total'])
def pair_tensor(T, A, B):
    tot = 0
    for (u, v), c in T.items():
        if u in A and v in B:
            tot = tot + c*A[u]*B[v]
    return tot
# build Omega explicitly
Omega = {}
for (u, v), c in Psi.items():
    for u2, cu in D3(u).items():
        for v2, cv in D2(v).items():
            key = (u2, v2)
            Omega[key] = Omega.get(key, 0) + c*(cu*cv)
print("Omega entries:", len(Omega))
ok = True
for (_, A) in L3:
    for (_, B) in L2:
        p1 = OB.pair(Psi, A, B)
        p2 = pair_tensor(Omega, A, B)
        if p1 != p2:
            ok = False
            print("MISMATCH:", p1, p2)
print("pairings via Psi == pairings via Dynkin-projected Omega:", ok)
