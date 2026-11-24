import sys, json, re
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor5_sym import *
import cor4_sym as c4

base = "<project>/research"
R3 = json.load(open(base + "/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, Fr(int(num), int(den))
terms = [parse(n, c) for n, c in R3]
w4 = [(e, w, c) for e, w, c in terms if len(w) == 4]
Z = (Fr(0),)*4
def vec(j, val=1):
    v = list(Z); v[j] = Fr(val); return tuple(v)
G_c = {0: vec(1), 1: lneg(vec(2)), -1: vec(3)}
letter_map = [(Fr(1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(1))]
F_half = {}
for e, w, cf in w4:
    for word, c0 in c4.hpl_sym(w, G_c).items():
        nw = tuple(lincomb(l, letter_map) for l in word)
        F_half[nw] = F_half.get(nw, Fr(0)) + cf * c0 * (Fr(2) ** (-e))
F_half = {w: c for w, c in F_half.items() if c}
d31_F = c4.delta31(F_half, 2)
print("delta_31(F(1/2)) [from F_symbol.json pipeline]: terms", len(d31_F))

# delta_31 of rel2's RHS symbol
e1 = (Fr(1), Fr(0)); e2 = (Fr(0), Fr(1))
def li4s(rv, om):  # S(Li4(r)) = -(om)[rv]^3 with signs inside vectors
    return {(lneg(om), rv, rv, rv): Fr(1)}
S = {}
def add(S0, q):
    for w, c in S0.items():
        S[w] = S.get(w, 0) + q*c
add(li4s(e1, e2), Fr(21))            # -21*Li4(-2): S(Li4(-2)) = -(e2,e1,e1,e1), coeff -21 -> +21(e2...)
add(li4s(lneg(e2), ladd(e1, lneg(e2))), Fr(65,8))   # Li4(1/3): rv=[1/3]=-e2, om=[2/3]=e1-e2: S=-(e1-e2)(-e2)^3 = +(e1-e2)(e2)^3
add(li4s(ladd(e1, lneg(e2)), lneg(e2)), Fr(-13,2))  # Li4(2/3): rv=e1-e2, om=-e2: S=-(-e2)(e1-e2)^3, coeff 13/2 -> -13/2
add(li4s(lneg(e1), lneg(e1)), Fr(-137,6))   # Li4(1/2): S=-(-e1)(-e1)^3=(e1)^4, coeff 137/6 -> wait S(Li4(1/2)) = -[1/2][1/2]^3 = -(-e1)^4 = -(e1^4); coeff in F4 is +137/6
add(li4s(lscale(2, lneg(e2)), ladd(lscale(3, e1), lscale(2, lneg(e2)))), Fr(-39,64))  # Li4(1/9): rv=[1/9]=-2e2, om=[8/9]=3e1-2e2
# products with nonzero delta31: Li3(1/3)*log2, Li3(-1/2)*log2, Li3(-1/2)*log3, zeta3*log2, zeta3*log3
# S(Li3(1/3)) = -[2/3][1/3]^2 = -((e1-e2),(-e2),(-e2)) 
S_Li313 = {(ladd(e1, lneg(e2)), lneg(e2), lneg(e2)): Fr(-1)}
# S(Li3(-1/2)) = -[3/2][1/2]^2 = -((-e1+e2),(-e1),(-e1))
S_Li3m12 = {(ladd(lneg(e1), e2), lneg(e1), lneg(e1)): Fr(-1)}
prods = [(S_Li313, e1, Fr(15,2)), (S_Li3m12, e1, Fr(37,2)), (S_Li3m12, e2, Fr(-5,4))]
for S3, l1, q in prods:
    for w, c in sshuffle(S3, {(l1,): Fr(1)}).items():
        S[w] = S.get(w, 0) + q*c
S = {w: c for w, c in S.items() if c}
d31_R = c4.delta31(S, 2)
print("delta_31(rel2 RHS): terms", len(d31_R))
# compare
keys = set(d31_F) | set(d31_R)
diff = {k: d31_F.get(k, 0) - d31_R.get(k, 0) for k in keys}
diff = {k: v for k, v in diff.items() if v}
print("delta_31(F) - delta_31(rel2):", "ZERO (consistent)" if not diff else f"{len(diff)} differing comps")
for k, v in sorted(diff.items())[:10]:
    print("   ", k, v)
