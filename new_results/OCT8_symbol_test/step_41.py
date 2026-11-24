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

for name in ["Lam1","Lam3"]:
    rem = Lam[name] - RHS(name)
    basis = [z3, z2*logv[2], z2*logv[3], z2*logv[5], mp.pi**3]
    rel = mp.pslq([rem]+basis, tol=mp.mpf('1e-400'), maxcoeff=10**60, maxsteps=1000000)
    print(name, "|rem| =", mp.nstr(abs(rem),3), " PSLQ(60):", rel)
    # also single-term checks
    print("   rem/z3 =", mp.nstr(rem/z3, 15), " rem/pi^3 =", mp.nstr(rem/mp.pi**3, 15))
