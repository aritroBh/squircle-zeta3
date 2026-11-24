import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 450
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
d = pickle.load(open("symbol_identities_B.pkl","rb"))
G_names, D_names, sol_out = d["G_names"], d["D_names"], d["sol_out"]
Lam = {k: mp.mpf(v) for k, v in pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"].items()}
logv = {p: mp.log(p) for p in (2,3,5)}
z3, z2 = mp.zeta(3), mp.zeta(2)

def RHS(name):
    c = [Fr(x) for x in sol_out[name]]
    tot = mp.mpf(0)
    for i, gn in enumerate(G_names):
        if c[i]: tot += M(c[i])*mp.polylog(3, M(Fr(gn[4:-1])))
    for j, dn in enumerate(D_names):
        cf = c[len(G_names)+j]
        if not cf: continue
        m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
        w, p = Fr(m.group(1)), int(m.group(2))
        tot += M(cf)*(-mp.polylog(2, M(w)))*logv[p]
    return tot

results = {}
for name in ["Lam1","Lam2","Lam3","Lam4"]:
    rem = Lam[name] - RHS(name)
    vec = [rem, z3, z2*logv[2], z2*logv[3], z2*logv[5]]
    rel = mp.pslq(vec, tol=mp.mpf('1e-390'), maxcoeff=10**40, maxsteps=1000000)
    print(name, "|rem| =", mp.nstr(abs(rem),3))
    print("   PSLQ:", rel)
    results[name] = rel
    if rel:
        q = [Fr(int(r), int(rel[0])) for r in rel[1:]]
        print("   const =", -q[0], "*z3 +", -q[1], "*z2l2 +", -q[2], "*z2l3 +", -q[3], "*z2l5")
pickle.dump(results, open("verify_B_full.pkl","wb"))
