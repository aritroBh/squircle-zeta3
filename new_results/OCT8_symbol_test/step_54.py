import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 350
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
NS = pickle.load(open("nullspace.pkl","rb"))["NS"]
dB = pickle.load(open("symbol_identities_B.pkl","rb"))
G_names, D_names = dB["G_names"], dB["D_names"]
logv = {p: mp.log(p) for p in (2,3,5)}
Li3v = [mp.polylog(3, M(Fr(gn[4:-1]))) for gn in G_names]
Dv = []
for dn in D_names:
    m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
    w, p = Fr(m.group(1)), int(m.group(2))
    Dv.append(-mp.polylog(2, M(w))*logv[p])
z3, z2 = mp.zeta(3), mp.zeta(2)
inb = [z3, z2*logv[2], z2*logv[3], z2*logv[5]]
rev = []
for nv in NS:
    tot = mp.mpc(0)
    for i in range(len(G_names)):
        if nv[i] != '0': tot += M(Fr(nv[i]))*Li3v[i]
    for j in range(len(D_names)):
        if nv[len(G_names)+j] != '0': tot += M(Fr(nv[len(G_names)+j]))*Dv[j]
    rev.append(tot)
re_nonzero = [mp.re(v) for v in rev if abs(mp.re(v)) > mp.mpf('1e-300')]
print("#nonzero real values:", len(re_nonzero))
outside = []
in_count = 0
for v in re_nonzero:
    rel = mp.pslq([v]+inb, tol=mp.mpf('1e-320'), maxcoeff=10**80, maxsteps=10**6)
    if rel is None: outside.append(v)
    else: in_count += 1
print("in {z3,z2log}:", in_count, " outside:", len(outside))
for v in outside[:10]:
    print("   OUT:", mp.nstr(v, 18))
# check remB against outside vals + inb
Lam = {k: mp.mpf(v) for k, v in pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"].items()}
solB = dB["sol_out"] if "sol_out" in dB else dB["sol"]
c = [Fr(x) for x in solB["Lam1"]]
tot = mp.mpf(0)
for i, gn in enumerate(G_names):
    if c[i]: tot += M(c[i])*mp.re(Li3v[i])
for j, dn in enumerate(D_names):
    cf = c[len(G_names)+j]
    if cf: tot += M(cf)*Dv[j]
remB = Lam["Lam1"] - tot
print("remB =", mp.nstr(remB, 15))
rel = mp.pslq([remB]+inb+outside, tol=mp.mpf('1e-320'), maxcoeff=10**80, maxsteps=10**6)
print("remB PSLQ vs {z3,z2log}+outside:", rel)
pickle.dump({"outside":[mp.nstr(v,350) for v in outside]}, open("outside_vals.pkl","wb"))
