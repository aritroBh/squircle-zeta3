"""w6 identification step 22: joint LLL of the 12 basic integrals N (dps-900 run) with the large template family, then an exact test
whether the exact coefficient vector rho^(m) of c_m (m=4,3) lies in the Q-span of the N-parts of the found relations (i.e. c_m in the family).
usage: s29_joint_N.py GROUPS D   (groups as in s26: G,L2,A,T3,T4,T5,T6,AA)       Printed numbers from this script."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from w6lib import *
from w6lll import lll_relations
import sympy
groups = sys.argv[1].split(","); D = int(sys.argv[2])
src = open(os.path.join(W, "s26_c3_bigfamily.py")).read()
# reuse the family builder of s26 (everything between the argv parsing and the targets)
pre = src.split('groups = sys.argv[1].split(",")')[1]
exec("RUN = 900\n" + pre.split("cs = [mp.mpf(coef(RUN, j))")[0].split("\n", 1)[1].replace("mp.mp.dps = D + 30", "mp.mp.dps = D + 30"))
names_fam = [n for n, _ in B]
N = json.load(open(os.path.join(W, "Nints_900.json")))["N"]
mp.mp.dps = D + 30
Nnames = list(N.keys())
vals = [mp.mpf(N[k]) for k in Nnames] + [v for _, v in B]
names = [f"N[{k}]" for k in Nnames] + names_fam
print(f"{len(vals)} numbers ({len(Nnames)} N + {len(B)} family), D={D}, mp.mp.dps={mp.mp.dps}, spurious size ~10^{D/len(vals):.1f}", flush=True)
rows = lll_relations(vals, D)
good = [(rel, res) for rel, res, last in rows if abs(res) < mp.mpf(10) ** (-(D - 25))]
print("rows with residual < 1e-%d: %d" % (D - 25, len(good)))
cred = [(rel, res) for rel, res in good if max(abs(a) for a in rel) < 10 ** (D / len(vals) - 2)]
print("credible relations (max coeff < 10^(spurious-2)):", len(cred))
for rel, res in cred:
    print("  max|coeff|=%d  residual %s" % (max(abs(a) for a in rel), mp.nstr(res, 3)))
    print("    N-part:", [(Nnames[i], a) for i, a in enumerate(rel[:12]) if a != 0])
    print("    family-part:", [(names_fam[i - 12], a) for i, a in enumerate(rel) if i >= 12 and a != 0][:40])
# exact test of rho^(m)
rho = json.load(open(os.path.join(W, "rho_exact.json")))
def rho_vec(m):
    v = [Fr(0)] * 12
    for name, co in rho.items():
        key, let, comp = name.split("|")
        if comp == "1": continue            # sqrt2-components only enter m <= 2
        v[Nnames.index(f"{key}|{let}")] += Fr(co[m])
    return v
Mrows = [[Fr(a) for a in rel[:12]] for rel, res in cred]
if Mrows:
    M = sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in r] for r in Mrows])
    r0 = M.rank()
    print("rank of N-parts of credible relations:", r0, "(so dim of N-span modulo family <=", 12 - r0, ")")
    for m in (4, 3):
        v = rho_vec(m)
        M2 = M.col_join(sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in v]]))
        print(f"rho^({m}) in the span of the relation N-parts: {M2.rank() == r0}   (rank with rho^({m}) appended: {M2.rank()})")
json.dump({"groups": groups, "D": D, "relations": [{"rel": [int(a) for a in rel], "max": int(max(abs(a) for a in rel))} for rel, res in cred],
           "names": names}, open(os.path.join(W, "joint_relations_%s.json" % "_".join(groups)), "w"))
