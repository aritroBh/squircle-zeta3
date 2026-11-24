"""a01_obstruction_check.py -- viability check for Route A (function-level identity in c).

Reads G1's fitted B3 (M3 = -z3*Phi + z2*Psi + R3, HPLs of c with letters {0,1,-1}),
extracts the weight-4 part, splits by the prefactor c^e, and forms the two pure
functions whose values at c = 1/sqrt2 are the rational and sqrt2 parts:
    U(c) = F_0(c) + 2 F_{-2}(c),     V(c) = F_{-1}(c) + F_1(c)/2.
Computes the (function) symbol in the alphabet {c, 1-c, 1+c} and the
Lambda^2 B_2 obstruction delta_{2,2} (Goncharov / Duhr-Gangl-Rhodes) to a
classical (Li_4 + products) representation.  Exact rational arithmetic.
NOTE: B3 itself is (at this point) only EXACT through c^640 + numerical; this
script only decides whether the function-level route can work.
"""
import json, sys
from fractions import Fraction as Fr
from collections import defaultdict

G1 = "research/10_frontier/G1_zeta3_any_angle/"
fit = json.load(open(G1 + "hpl_fit_k3_(0,0,0).json"))

def parse(name):
    # 'c^-1 H[0, 1, 0, 1]'
    pre, w = name.split(' H')
    e = int(pre.replace('c^', ''))
    word = tuple(int(x) for x in w.strip('[]').split(','))
    return e, word

F = defaultdict(lambda: defaultdict(Fr))   # F[e][word] = coeff  (weight-4 only)
for name, v in fit:
    if name in ('1', 'c'):
        continue
    e, w = parse(name)
    if len(w) == 4:
        F[e][w] += Fr(v)

print("weight-4 prefactor exponents present:", sorted(F))
U = defaultdict(Fr); V = defaultdict(Fr)
for w, v in F.get(0, {}).items(): U[w] += v
for w, v in F.get(-2, {}).items(): U[w] += 2*v
for w, v in F.get(-1, {}).items(): V[w] += v
for w, v in F.get(1, {}).items(): V[w] += v/2

# symbol: S(H_{a1..an}) = s(an) x ... x s(a1); s(0)=[c], s(1)=-[1-c], s(-1)=[1+c]
IDX = {0: 0, 1: 1, -1: 2}
SGN = {0: 1, 1: -1, -1: 1}
def symbol(lin):
    out = defaultdict(Fr)
    for w, v in lin.items():
        sg = 1
        for a in w: sg *= SGN[a]
        key = tuple(IDX[a] for a in reversed(w))
        out[key] += v*sg
    return {k: x for k, x in out.items() if x}

def wedge(a, b):
    if a == b: return None, 0
    return ((a, b), 1) if a < b else ((b, a), -1)

def delta22(S):
    out = defaultdict(Fr)
    for (a, b, c, d), v in S.items():
        p, s1 = wedge(a, b); q, s2 = wedge(c, d)
        if s1 and s2:
            key, s3 = wedge(p, q)
            if s3: out[key] += v*s1*s2*s3
    return {k: x for k, x in out.items() if x}

for name, lin in (("U", U), ("V", V)) + tuple((f"F_{e}", F[e]) for e in sorted(F)):
    S = symbol(lin)
    print(f"{name}: {len(lin)} HPL terms, {len(S)} symbol terms, delta22 = {delta22(S)}")
print("alphabet: 0=c, 1=1-c, 2=1+c")
