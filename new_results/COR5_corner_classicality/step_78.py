import sys, json
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
sys.path.insert(0, ".")
from fractions import Fraction as Fr
import cor5_sym as cs
import cor4_sym as c4

base = "<project>/research"
# rebuild F_sym from R3 with cor4 conventions
R3 = json.load(open(base + "/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
import re
def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1))
    w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, Fr(int(num), int(den))
terms = [parse(n, c) for n, c in R3]
w4 = [(e, w, c) for e, w, c in terms if len(w) == 4]
Z = (Fr(0),)*4
def vec(j, val=1):
    v = list(Z); v[j] = Fr(val); return tuple(v)
G_c = {0: vec(1), 1: c4.lneg(vec(2)), -1: vec(3)}
F_sym = {}
for e, w, cf in w4:
    for word, c0 in c4.hpl_sym(w, G_c).items():
        F_sym[(e, word)] = F_sym.get((e, word), Fr(0)) + cf * c0
F_sym = {k: v for k, v in F_sym.items() if v}
# delta_22 as a function of c
d22 = c4.delta22L2_func(F_sym, 4)
BASIS = ["2", "c", "1-c", "1+c"]
for (p1, p2), poly in sorted(d22.items()):
    print(f"  ({BASIS[p1[0]]}^{BASIS[p1[1]]}) x ({BASIS[p2[0]]}^{BASIS[p2[1]]}): "
          + " + ".join(f"{v} c^{e}" for e, v in sorted(poly.items())))
