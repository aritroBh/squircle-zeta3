import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 420
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator

d = pickle.load(open("symbol_identities_B.pkl","rb"))
G_names, D_names, sol_out = d["G_names"], d["D_names"], d["sol_out"]
Lam = {k: mp.mpf(v) for k, v in pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"].items()}

z3 = mp.zeta(3); z2 = mp.zeta(2)
logv = {p: mp.log(p) for p in (2,3,5)}

def RHS(name):
    c = [Fr(x) for x in sol_out[name]]
    tot = mp.mpf(0)
    for i, gn in enumerate(G_names):
        if c[i] == 0: continue
        z = Fr(gn[4:-1])
        tot += M(c[i]) * mp.polylog(3, M(z))
    for j, dn in enumerate(D_names):
        cf = c[len(G_names)+j]
        if cf == 0: continue
        if dn.startswith("-Li2("):
            m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
            w, p = Fr(m.group(1)), int(m.group(2))
            tot += M(cf) * (-mp.polylog(2, M(w))) * logv[p]
        else:
            ps = [int(x[3:]) for x in dn.split("*")]
            pr = mp.mpf(1)
            for p in ps: pr *= logv[p]
            tot += M(cf) * pr
    return tot

for name in ["Lam2","Lam4","Lam1","Lam3"]:
    r = RHS(name)
    rem = Lam[name] - r
    # PSLQ remainder vs {z3, z2*log2, z2*log3, z2*log5, 1}
    vec = [rem, z3, z2*logv[2], z2*logv[3], z2*logv[5]]
    rel = mp.pslq(vec, tol=mp.mpf('1e-380'), maxcoeff=10**12, maxsteps=10000)
    print(name, "|remainder| =", mp.nstr(abs(rem),3), " PSLQ:", rel)
    if rel:
        print("   => remainder =", -rel[1]/rel[0], "*z3 +", -rel[2]/rel[0], "*z2l2 +",
              -rel[3]/rel[0], "*z2l3 +", -rel[4]/rel[0], "*z2l5")
