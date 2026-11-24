import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 500
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
dA = pickle.load(open("symbol_identities.pkl","rb"))
dB = pickle.load(open("symbol_identities_B.pkl","rb"))
logv = {p: mp.log(p) for p in (2,3,5)}
z3, z2 = mp.zeta(3), mp.zeta(2)
Lam = {k: mp.mpf(v) for k, v in pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"].items()}

def RHS(d, name):
    G_names, D_names, sol_out = d["G_names"], d["D_names"], d["sol_out"]
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

def Li21A(x, y, N=3000):
    x = M(x); y = M(y); s = mp.mpf(0); yn = y
    for n in range(2, N):
        yn *= y
        t = mp.mpf(0); xm = x
        for m in range(1, n):
            t += xm/m**2; xm *= x
        s += t*yn/n
        if abs(yn) < mp.mpf(10)**(-mp.mp.dps-10): break
    return s

lhsA = Li21A(Fr(-1,3), Fr(1,3))
lhsB = Lam["Lam1"]
remA = lhsA - RHS(dA, "Lam1")
remB = lhsB - RHS(dB, "Lam1")
print("remA =", mp.nstr(remA, 15))
print("remB =", mp.nstr(remB, 15))
print("remA+remB =", mp.nstr(remA+remB, 15))
basis = [z3, z2*logv[2], z2*logv[3], z2*logv[5]]
for nm, r in [("remA",remA), ("remB",remB), ("remA+remB",remA+remB)]:
    rel = mp.pslq([r]+basis, tol=mp.mpf('1e-450'), maxcoeff=10**80, maxsteps=1000000)
    print(nm, "PSLQ:", rel)
