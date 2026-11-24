"""step11_F_symbol.py -- Task 2: the symbol of F (weight-4 part of R3) as a function of c,
and its antisymmetrised (2,2) component (Goncharov delta_{2,2}, symbol level)."""
import sys, json, re
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
R3 = json.load(open(G1 + "hpl_fit_k3_(0,0,0).json"))

def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1))
    w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, Fr(int(num), int(den))

terms = [parse(n, c) for n, c in R3]
w4 = [(e, w, c) for e, w, c in terms if len(w) == 4]
w23 = [(e, w, c) for e, w, c in terms if len(w) < 4]
print(f"R3: {len(terms)} terms total; weight-4: {len(w4)}; weight<4: {len(w23)}")
from collections import Counter
print("word multiplicities at weight 4:", len(set(w for _, w, _ in w4)))

# letter basis: index 0=[2], 1=[c], 2=[1-c], 3=[1+c]
BASIS = ["2", "c", "1-c", "1+c"]
Z = (Fr(0),)*4
def vec(j, val=1):
    v = list(Z); v[j] = Fr(val); return tuple(v)
G_c  = {0: vec(1), 1: lneg(vec(2)), -1: vec(3)}        # argument  c
G_mc = {0: vec(1), 1: lneg(vec(3)), -1: vec(2)}        # argument -c

# F as a function symbol: dict (e, word) -> Fr
F_sym = {}
for e, w, cf in w4:
    s = hpl_sym(w, G_c)
    for word, c0 in s.items():
        F_sym[(e, word)] = F_sym.get((e, word), 0) + cf * c0
F_sym = {k: v for k, v in F_sym.items() if v}
print(f"function symbol of F: {len(F_sym)} terms")
for (e, w), cf in sorted(F_sym.items()):
    def lstr(v):
        parts = []
        for j, b in enumerate(BASIS):
            if v[j]:
                parts.append(f"{b}^{v[j]}" if v[j] != 1 else b)
        return "*".join(parts) if parts else "1"
    print(f"  c^{e:<2} * {cf!s:>4} * {tuple(lstr(v) for v in w)}")

json.dump({f"c^{e} " + "|".join("".join(str(int(x)) for x in v) for v in w): str(cf)
           for (e, w), cf in F_sym.items()},
          open("<project>/research/08_theory/COR4_M3alpha_classical/F_symbol.json", "w"), indent=1)

# ---- the (2,2) component as a Laurent polynomial in c for each wedge pair ----
d22 = delta22L2_func(F_sym, 4)
print()
print("delta_{2,2}(F) -- components ((i,j),(k,l)) with letters", BASIS)
if not d22:
    print("  == 0 : no (2,2) obstruction; F is a classical candidate")
for (p1, p2), poly in sorted(d22.items()):
    print(f"  ({BASIS[p1[0]]}^{BASIS[p1[1]]}) x ({BASIS[p2[0]]}^{BASIS[p2[1]]}): "
          + " + ".join(f"{v} c^{e}" for e, v in sorted(poly.items())))
print()
print("number of nonzero wedge components:", len(d22))
nonzero = len(d22) > 0
print("VERDICT:", "delta22(F) != 0 -> F is NOT classical as a function of c (symbol level)"
      if nonzero else "delta22(F) = 0")
# also the (3,1) and (1,3) components
for ee in (-1, 0, 1):
    Se = {w: cf for (e, w), cf in F_sym.items() if e == ee}
    d31 = delta31(Se, 4)
    d13 = delta13(Se, 4)
    print(f"(3,1) c^{ee}: {'nonzero, %d terms' % len(d31) if d31 else '0'} ; (1,3) c^{ee}: {'nonzero, %d terms' % len(d13) if d13 else '0'}")
