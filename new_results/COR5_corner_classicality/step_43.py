import sys, json
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
import cor5_sym as cs
import cor4_sym as c4

# ---------- Independent rebuild of S(F(1/2)) from R3 directly ----------
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
# letters over basis (2, c, 1-c, 1+c) as in cor4
Z = (Fr(0),)*4
def vec(j, val=1):
    v = list(Z); v[j] = Fr(val); return tuple(v)
G_c = {0: vec(1), 1: cs.lneg(vec(2)), -1: vec(3)}
F_sym = {}
for e, w, cf in w4:
    for word, c0 in cs.hpl_sym(w, G_c).items():
        F_sym[(e, word)] = F_sym.get((e, word), Fr(0)) + cf * c0
F_sym = {k: v for k, v in F_sym.items() if v}
print("rebuilt F_sym terms:", len(F_sym))
# compare with F_symbol.json
Fjson = cs.parse_F_symbol(cor4 := base + "/08_theory/COR4_M3alpha_classical" + "/F_symbol.json")
same = {(e, w): cf for e, w, cf in Fjson} == F_sym
print("matches F_symbol.json:", same)

# specialize to c = 1/2
letter_map = [(Fr(1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(1))]
F_half = {}
for (e, word), cf in F_sym.items():
    cval = Fr(2) ** (-e)
    nw = tuple(cs.lincomb(l, letter_map) for l in word)
    F_half[nw] = F_half.get(nw, Fr(0)) + cf * cval
F_half = {w: c for w, c in F_half.items() if c}

# delta_{2,2} at c=1/2 via cor4 machinery (Lambda^2 of Lambda^2 V), dim=2
d22 = c4.delta22L2_spec(F_half, 2)
print("delta_{2,2}(F(1/2)) =", d22 if d22 else "0  (automatic: dim Lambda^2 V = 1)")

# quotient coords
Q = cs.W4Quotient(2)
vF = Q.coords(F_half)
print("pi(S(F(1/2))):", vF)
# the four named combinations
sols = json.load(open("sunit23_solutions.json"))
V23 = {Fr(xs): (Fr(fx[1]), Fr(fx[2])) for xs, fx, fy in sols}
OM23 = {Fr(xs): (Fr(fy[1]), Fr(fy[2])) for xs, fx, fy in sols}
def li4v(x): return Q.coords(cs.li_sym(4, V23[x], OM23[x]))
combos = [
    ([Fr(-2), Fr(1,3), Fr(2,3)], [Fr(-2), Fr(1,2), Fr(1,2)]),
    ([Fr(-1,2), Fr(1,3), Fr(2,3)], [Fr(2), Fr(1,2), Fr(1,2)]),
    ([Fr(1,3), Fr(2,3), Fr(3,4)], [Fr(-1,2), Fr(-3,2), Fr(-1,2)]),
    ([Fr(1,4), Fr(1,3), Fr(2,3)], [Fr(1,4), Fr(1,2), Fr(1,2)]),
]
for args, cfs in combos:
    tot = tuple(sum(cf*v[i] for (cf, v) in zip(cfs, [li4v(a) for a in args])) for i in range(3))
    print(f"  {args} {cfs}: combo == vF: {tot == vF}")
