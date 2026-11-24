import sys, json, re, itertools
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
R2 = json.load(open(G1 + "hpl_fit_k2_(0,0,0).json"))
def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, Fr(int(num), int(den))
w3 = [parse(n, c) for n, c in R2 if len(parse(n, c)[1]) == 3]
print("weight-3 terms of R2:", w3)

h = Fr(1,2)
g = {0: (-h, Fr(0)), 1: (h, Fr(1)), -1: (-h, Fr(1))}
def cmul2(a, b): return (a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
c0pow = {-1: (Fr(0), Fr(1)), 0: (Fr(1), Fr(0)), 1: (Fr(0), h), 2: (Fr(1,2), Fr(0)), -2: (Fr(2), Fr(0))}
SF2 = {}
for e, w, cf in w3:
    word = tuple(g[a] for a in reversed(w))
    cc = cmul2(c0pow[e], (cf, Fr(0)))
    cur = SF2.get(word, (Fr(0),Fr(0)))
    SF2[word] = (cur[0]+cc[0], cur[1]+cc[1])
SF2 = {w: c for w, c in SF2.items() if c != (Fr(0),Fr(0))}
print("S(F2)(c0) terms:")
for w, c in sorted(SF2.items()): print("   ", w, c)

# classical Li3 library at c0: li_sym(3, rv, omrv) = -[1-r] [r]^2
Z2 = (Fr(0), Fr(0))
lib = {
 'Li3(b)':    ((Fr(0),Fr(-1)), (h,Fr(-1))),
 'Li3(-b)':   ((Fr(0),Fr(-1)), (h,Fr(0))),
 'Li3(1/sq2)':((-h,Fr(0)), (-h,Fr(-1))),
 'Li3(1/2)':  ((Fr(-1),Fr(0)), (Fr(-1),Fr(0))),
 'Li3(1-1/sq2)':((-h,Fr(-1)), (-h,Fr(0))),
 'Li3(2-sq2)':((h,Fr(-1)), (Fr(0),Fr(-1))),
 'Li3(-b^2)': ((Fr(0),Fr(-2)), (Fr(3,2),Fr(-1))),
 'Li3(b^2)':  ((Fr(0),Fr(-2)), (-h,Fr(1))),   # 1-b^2 = 1-(3-2sqrt2) = 2sqrt2-2 = 2(√2-1) = 2/(1+√2): [1-b^2]=(1,-1)
}
# wait check 1-b^2: b^2 = (√2-1)^2 = 3-2√2; 1-b^2 = 2√2-2 = 2(√2-1): [2(√2-1)] = (1,-1): so omrv=(1,-1), and -[1-b^2] = (-1,1). Let me just use li_sym directly which negates.
