import mpmath as mp, pickle
from fractions import Fraction as Fr
mp.mp.dps = 300
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
dB = pickle.load(open("symbol_identities_B.pkl","rb"))
G_names, D_names, solB = dB["G_names"], dB["D_names"], dB["sol_out"]
c = [Fr(x) for x in solB["Lam1"]]
print("nonzero coeffs:", sum(1 for x in c if x != 0))
tot = mp.mpf(0)
print("--- Li3 terms ---")
for i, gn in enumerate(G_names):
    if c[i]:
        z = Fr(gn[4:-1])
        v = M(c[i])*mp.polylog(3, M(z))
        tot += v
        print(f"{gn:>14} coeff {str(c[i]):>12}  term {mp.nstr(v,12)}")
print("--- prod terms ---")
import re
for j, dn in enumerate(D_names):
    cf = c[len(G_names)+j]
    if cf:
        m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
        w, p = Fr(m.group(1)), int(m.group(2))
        v = M(cf)*(-mp.polylog(2, M(w)))*mp.log(p)
        tot += v
        print(f"{dn:>20} coeff {str(cf):>10}  term {mp.nstr(v,12)}")
print("TOTAL RHS:", mp.nstr(tot, 20))
Lam = {k: mp.mpf(v) for k, v in pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"].items()}
print("LHS:", mp.nstr(Lam["Lam1"], 20))
print("rem:", mp.nstr(Lam["Lam1"]-tot, 20))
