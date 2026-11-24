
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json
import mpmath as mp
import w5lib as W
import pipe
D=5
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
RM = json.load(open(pipe.cand_dir(D)+"/rel_min_atom.json"))
names = RM["names"]; c = [Fraction(x) for x in RM["c"]]
cmap = dict(zip(names, c))
used = [(nm, cf) for nm, cf in cmap.items() if nm != "T" and cf]
print("terms:")
for nm, cf in used: print(f"  {nm:>26s}  {cf}")
# verify numerically fresh at 1200 digits
P = 1200
mp.mp.dps = P + 30
lphi = mp.log((1+mp.sqrt(5))/2); ls5 = mp.log(mp.sqrt(5)); l2 = mp.log(2); pi2 = mp.pi**2
logmap = {"lphi": lphi, "ls5": ls5, "l2": l2}
tot = W.Lchi_hurwitz(D, P)/mp.sqrt(5)   # T = L/sqrt(5)
for nm, cf in cmap.items():
    if nm == "T" or not cf: continue
    if nm == "z5": v = mp.zeta(5)
    elif nm.startswith("Li5("):
        i = [j for j, cd in enumerate(cands) if f"Li5({pipe.fmt_qd(cd['x'])})" == nm][0]
        v = W.li5_series(cands[i]['x'].emb(P+30), P)
    else:
        v = mp.mpf(1)
        for f in nm.split("*"): v *= {"lphi": lphi, "ls5": ls5, "l2": l2, "pi2": pi2}[f]
    tot -= cf*v
print("fresh 1200-digit residual T - RHS:", mp.nstr(tot, 4))
# also evaluate with the CACHED 1150-digit values
vals = pipe.load_vals(D, 1150)
mp.mp.dps = 1200
tot2 = vals["T"]
prods = pipe.prod_values(C, vals)
vmap = {"z5": vals["z5"]}
for i, cd in enumerate(cands): vmap[f"Li5({pipe.fmt_qd(cd['x'])})"] = vals[f"li5:{i}"]
for kind, msk, nm2, pv in prods: vmap[nm2] = pv
for nm, cf in cmap.items():
    if nm == "T" or not cf: continue
    tot2 -= cf*vmap[nm]
print("cached 1150-digit residual:", mp.nstr(tot2, 4))
