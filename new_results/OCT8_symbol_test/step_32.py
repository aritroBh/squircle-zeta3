import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 100
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
d = pickle.load(open("symbol_identities_B.pkl","rb"))
G_names, D_names, sol_out = d["G_names"], d["D_names"], d["sol_out"]
print("len G,D, coeffs:", len(G_names), len(D_names), len(sol_out["Lam1"]))
logv = {p: mp.log(p) for p in (2,3,5)}

name = "Lam1"
c = [Fr(x) for x in sol_out[name]]
sG = mp.mpf(0)
for i, gn in enumerate(G_names):
    if c[i] == 0: continue
    z = Fr(gn[4:-1])
    sG += M(c[i]) * mp.polylog(3, M(z))
sD = mp.mpf(0)
for j, dn in enumerate(D_names):
    cf = c[len(G_names)+j]
    if cf == 0: continue
    m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
    w, p = Fr(m.group(1)), int(m.group(2))
    sD += M(cf) * (-mp.polylog(2, M(w))) * logv[p]
Lam = {k: mp.mpf(v) for k, v in pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"].items()}
print("sum Li3 part:", mp.nstr(sG, 15))
print("sum prod part:", mp.nstr(sD, 15))
print("RHS total:", mp.nstr(sG+sD, 15))
print("LHS:", mp.nstr(Lam[name], 15))
print("diff:", mp.nstr(abs(sG+sD-Lam[name]),3))
