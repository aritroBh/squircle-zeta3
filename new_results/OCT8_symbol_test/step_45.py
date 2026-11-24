import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 500
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
dB = pickle.load(open("symbol_identities_B.pkl","rb"))
Lam = {k: mp.mpf(v) for k, v in pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"].items()}
logv = {p: mp.log(p) for p in (2,3,5)}
z3, z2 = mp.zeta(3), mp.zeta(2)

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
remB = Lam["Lam1"] - RHS(dB, "Lam1")
print("remB =", mp.nstr(remB, 25))

# LLL integer relation with python-flint
import flint
from flint import fmpz, fmpz_mat
def int_rel(vals, digits, maxcoeff_exp=None):
    scale = 10**digits
    n = len(vals)
    A = fmpz_mat(n, n+1)
    for i, v in enumerate(vals):
        A[i,0] = fmpz(int(mp.nint(v*scale)))
        A[i,i+1] = fmpz(1)
    B = A.lll()
    # first few rows are short vectors
    out = []
    for r in range(min(3,n)):
        row = [B[r,0]] + [B[r,i+1] for i in range(n)]
        if row[1:] and all(abs(x) < 10**10 for x in row[1:]):
            out.append(row)
    return out

# identify remB among Li3 values + constants
sols = pickle.load(open("sols235.pkl","rb"))
vals = [remB] + [mp.polylog(3, M(z)) for z in sols] + [z3, z2*logv[2], z2*logv[3], z2*logv[5]]
rels = int_rel(vals, 480)
print("LLL relations (remB, Li3's, z3, z2log p):")
for r in rels:
    print(r)
