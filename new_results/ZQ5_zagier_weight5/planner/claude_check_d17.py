# Planner (main session) check of the W4 session's L(chi_17,5) relation, with own code: mpmath polylog and Hurwitz zeta, dps 80.
import json, re, mpmath as mp
from sympy import legendre_symbol
mp.mp.dps = 80
d = json.load(open("<project>/research/08_theory/ZQ5_zagier_weight5/planner/d17_relation_planner.json"))
s17 = mp.sqrt(17)
def arg(e):
    e = e.replace("+-", "-").replace("sqrt17", "*s17")
    return eval(re.sub(r"(\d+)/(\d+)", r"mp.mpf(\1)/\2", e).replace("**s17", "*s17"), {"mp": mp, "s17": s17})
chi = lambda a: 0 if a % 17 == 0 else legendre_symbol(a % 17, 17)
L = mp.mpf(17)**-5 * mp.fsum(chi(a) * mp.zeta(5, mp.mpf(a)/17) for a in range(1, 17))
T = L / s17
var = {"lE17": mp.log(4 + s17), "lp": mp.log((3 + s17)/2), "lp'": mp.log(abs((3 - s17)/2)), "pi2": mp.pi**2}
big = 0
tot = int(d["T_coef"]) * T + int(d["z5_coef"]) * mp.zeta(5)
for x, c in d["li5"].items():
    v = arg(x)
    if abs(v) >= 1: big += 1
    tot += int(c) * mp.re(mp.polylog(5, v))
for m, c in d["products"].items():
    tot += int(c) * mp.fprod(var[f] for f in m.split("*"))
print("L(chi_17,5) =", mp.nstr(L, 30), "(W4 session: 1.027634038145162923246265...)")
print("arguments with |x| >= 1:", big, "of", len(d["li5"]))
print("relation residual at dps 80:", mp.nstr(tot, 5), "| size of the T term:", mp.nstr(int(d["T_coef"]) * T, 5))
