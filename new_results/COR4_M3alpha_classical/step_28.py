import sys
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *
import json, re

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
R3 = json.load(open(G1 + "hpl_fit_k3_(0,0,0).json"))
def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, Fr(int(num), int(den))
w4 = [parse(n, c) for n, c in R3 if len(parse(n, c)[1]) == 4]

def vec(dim, j, val=1):
    v = [Fr(0)]*dim; v[j] = Fr(val); return tuple(v)
G_c = {0: vec(4,1), 1: lneg(vec(4,2)), -1: vec(4,3)}
F_sym = {}
for e, w, cf in w4:
    for word, c0 in hpl_sym(w, G_c).items():
        F_sym[(e, word)] = F_sym.get((e, word), 0) + cf * c0
F_sym = {k: v for k, v in F_sym.items() if v}

# --- specialize to c0 = 1/sqrt(2): basis (2, 1+sqrt2) ---
# [2] -> (1,0); [c] -> (-1/2, 0); [1-c] -> (-1/2, -1); [1+c] -> (-1/2, 1)
h = Fr(1,2)
images = [(Fr(1),Fr(0)), (-h,Fr(0)), (-h,Fr(-1)), (-h,Fr(1))]
lmap = make_letter_map(images)
# c0^e as (p,q) = p + q sqrt2: c0 = sqrt2/2 -> (0, 1/2); c0^-1 = sqrt2 -> (0,1)
cpow = {-1: (Fr(0),Fr(1)), 0: (Fr(1),Fr(0)), 1: (Fr(0),Fr(1,2))}
S_spec = specialize_func(F_sym, lmap, cpow)
print("specialized symbol terms:", len(S_spec))
d = delta22_spec(S_spec, 2)
print("delta22 at c=1/sqrt2:", d if d else "= 0  ✓ (square corner classical, as required)")

# --- specialize to c0 = 1/2: basis (2,3) ---
# [2]->(1,0); [c]=[1/2]->(-1,0); [1-c]=[1/2]->(-1,0); [1+c]=[3/2]->(-1,1)
images = [(Fr(1),Fr(0)), (Fr(-1),Fr(0)), (Fr(-1),Fr(0)), (Fr(-1),Fr(1))]
lmap = make_letter_map(images)
cpow = {-1: Fr(2), 0: Fr(1), 1: Fr(1,2)}
S2 = specialize_func(F_sym, lmap, cpow)
d2 = delta22_spec(S2, 2)
print("delta22 at c=1/2 (alpha=2pi/3):", d2 if d2 else "= 0")

# --- rational point c0 = 1/3: basis (2,3) ---
# [c]=[1/3]->(0,-1); [1-c]=[2/3]->(1,-1); [1+c]=[4/3]->(2,-1)
images = [(Fr(1),Fr(0)), (Fr(0),Fr(-1)), (Fr(1),Fr(-1)), (Fr(2),Fr(-1))]
lmap = make_letter_map(images)
cpow = {-1: Fr(3), 0: Fr(1), 1: Fr(1,3)}
S3 = specialize_func(F_sym, lmap, cpow)
d3 = delta22_spec(S3, 2)
print("delta22 at c=1/3:", d3 if d3 else "= 0")

# --- c0 = 2/5: primes {2,5}: [c]=( -1? ) 2/5 -> (1,-1); [1-c]=[3/5] -> prime 3! need basis {2,3,5}
images = [(Fr(1),Fr(0),Fr(0)), (Fr(1),Fr(0),Fr(-1)), (Fr(0),Fr(1),Fr(-1)), (Fr(0),Fr(1),Fr(0))]
# [2]->(1,0,0); [c]=[2/5]->(1,0,-1); [1-c]=[3/5]->(0,1,-1); [1+c]=[7/5]-> prime 7! 
print("c0=2/5 needs primes {2,3,5,7}: recompute")
images = [(Fr(1),Fr(0),Fr(0),Fr(0)), (Fr(1),Fr(0),Fr(-1),Fr(0)), (Fr(0),Fr(1),Fr(-1),Fr(0)), (Fr(0),Fr(0),Fr(-1),Fr(1))]
# basis (2,3,5,7): [2]; [c]=2/5 -> (1,0,-1,0); [1-c]=3/5 -> (0,1,-1,0); [1+c]=7/5 -> (0,0,-1,1)
lmap = make_letter_map(images)
cpow = {-1: Fr(5,2), 0: Fr(1), 1: Fr(2,5)}
S25 = specialize_func(F_sym, lmap, cpow)
d25 = delta22_spec(S25, 4)
print("delta22 at c=2/5:", d25 if d25 else "= 0")
