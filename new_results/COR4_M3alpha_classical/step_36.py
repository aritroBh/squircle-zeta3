import sys, json, re
from fractions import Fraction as Fr
from itertools import product as iproduct

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
R3 = json.load(open(G1 + "hpl_fit_k3_(0,0,0).json"))
def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, Fr(int(num), int(den))

# letters at c0 = 1/sqrt2, basis (2, 1+sqrt2): index 0 = [2], 1 = [1+sqrt2]
# g(0) = [c] = (-1/2, 0); g(1) = -[1-c] = (1/2, 1); g(-1) = [1+c] = (-1/2, 1)
h = Fr(1, 2)
g = {0: (-h, Fr(0)), 1: (h, Fr(1)), -1: (-h, Fr(1))}

# coefficient scalar in Q(sqrt2): (p,q) = p + q*sqrt2
def cmul(a, b): return (a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
c0pow = {-1: (Fr(0), Fr(1)), 0: (Fr(1), Fr(0)), 1: (Fr(0), h)}

# build specialized symbol: word (4 letters) -> (p,q)
S = {}
for n, c in R3:
    e, w, cf = parse(n, c)
    if len(w) != 4: continue
    word = tuple(g[a] for a in reversed(w))   # S(H_w) = g(a_n)...g(a_1)
    coef = cmul(c0pow[e], (cf, Fr(0)))
    cur = S.get(word, (Fr(0), Fr(0)))
    S[word] = (cur[0] + coef[0], cur[1] + coef[1])
S = {w: c for w, c in S.items() if c != (Fr(0), Fr(0))}
print("specialized symbol terms:", len(S))
for w, c in sorted(S.items()):
    print("  ", w, c)

# delta22: (a,b,c,d) -> (a^b)x(c^d)
d = {}
for w, cf in S.items():
    a, b_, c_, d_ = w
    w1 = a[0]*b_[1] - a[1]*b_[0]
    w2 = c_[0]*d_[1] - c_[1]*d_[0]
    v = w1 * w2
    d[(0,1)] = (d.get((0,1), (Fr(0),Fr(0)))[0] + cf[0]*v, d.get((0,1), (Fr(0),Fr(0)))[1] + cf[1]*v)
print("delta22 (from scratch):", d)
