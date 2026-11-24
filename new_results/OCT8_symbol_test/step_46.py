import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 800
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
dB = pickle.load(open("symbol_identities_B.pkl","rb"))
G_names, D_names, solB = dB["G_names"], dB["D_names"], dB["sol_out"]
Lam = {k: mp.mpf(v) for k, v in pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"].items()}
logv = {p: mp.log(p) for p in (2,3,5)}
z3, z2 = mp.zeta(3), mp.zeta(2)

name = "Lam1"
c = [Fr(x) for x in solB[name]]
terms = []
tot = mp.mpf(0)
for i, gn in enumerate(G_names):
    if c[i]:
        v = M(c[i])*mp.polylog(3, M(Fr(gn[4:-1])))
        terms.append((gn, c[i], v)); tot += v
for j, dn in enumerate(D_names):
    cf = c[len(G_names)+j]
    if cf:
        m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
        w, p = Fr(m.group(1)), int(m.group(2))
        v = M(cf)*(-mp.polylog(2, M(w)))*logv[p]
        terms.append((dn, cf, v)); tot += v
remB = Lam[name] - tot
print("RHS =", mp.nstr(tot, 20))
print("remB =", mp.nstr(remB, 20))
# any complex terms?
print("any complex:", [t[0] for t in terms if mp.im(t[2]) != 0])
# PSLQ with maximal settings
rel = mp.pslq([remB, z3, z2*logv[2], z2*logv[3], z2*logv[5]], tol=mp.mpf('1e-750'), maxcoeff=10**300, maxsteps=10**7)
print("PSLQ(max):", rel)
# identify remB against the same Li3's in its RHS + constants
li3s = sorted(set(Fr(g[4:-1]) for g,_,_ in terms if g.startswith("Li3")))
vals = [remB] + [mp.polylog(3, M(z)) for z in li3s] + [z3, z2*logv[2], z2*logv[3], z2*logv[5]]
rel2 = mp.pslq(vals, tol=mp.mpf('1e-750'), maxcoeff=10**200, maxsteps=10**7)
print("PSLQ vs own Li3s+const:", rel2)
