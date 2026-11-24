import sys, json
sys.path.insert(0, ".")
from fractions import Fraction as Fr
import cor5_sym as cs

base = "<project>/research"
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
print("weight-4 terms:", len(w4))
Z = (Fr(0),)*4
def vec(j, val=1):
    v = list(Z); v[j] = Fr(val); return tuple(v)
G_c = {0: vec(1), 1: cs.lneg(vec(2)), -1: vec(3)}
# build function symbol, then specialize
letter_map = [(Fr(1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(1))]
F_half = {}
for e, w, cf in w4:
    for word, c0 in cs.hpl_sym(w, G_c).items():
        nw = tuple(cs.lincomb(l, letter_map) for l in word)
        cval = Fr(2) ** (-e)
        F_half[nw] = F_half.get(nw, Fr(0)) + cf * c0 * cval
F_half = {w: c for w, c in F_half.items() if c}
Q = cs.W4Quotient(2)
print("vF =", Q.coords(F_half))
# per-term contributions
print("\nper-term quotient contributions:")
tot = [Fr(0)]*3
for e, w, cf in sorted(w4):
    cval = Fr(2) ** (-e)
    for word, c0 in cs.hpl_sym(w, G_c).items():
        nw = tuple(cs.lincomb(l, letter_map) for l in word)
        qc = Q.coords({nw: cf * c0 * cval})
        if any(qc):
            print(f"  c^{e} {w} cf={cf}: {qc}")
        for j in range(3):
            tot[j] += qc[j]
print("sum:", tuple(tot))
