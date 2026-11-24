import sys
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *

h = Fr(1,2)
# basis e0 = [2], e1 = [1+sqrt2]
# symbol of paper Li4part at c0: li_sym(4, rv, omrv) = -[1-r] [r]^3, coefficients (p,q)=p+q sqrt2
def cmul2(a, b): return (a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def pair(x): return (Fr(x[0]), Fr(x[1]))
terms = [
    ((1,-573), (Fr(0),Fr(-1)), (h,Fr(-1))),      # Li4(b): [b]=(0,-1), [1-b]=(1/2,-1)
    ((-1,-573), (Fr(0),Fr(-1)), (h,Fr(0))),      # Li4(-b): [1+b]=[sqrt2]=(1/2,0)
    ((0,56), (Fr(0),Fr(-2)), (Fr(3,2),Fr(-1))),  # Li4(-b^2): [b^2]=(0,-2), [1+b^2]=(3/2,-1)
    ((0,17), (-h,Fr(0)), (-h,Fr(-1))),           # Li4(1/sqrt2)
    ((-Fr(37,4),-Fr(17,16)), (Fr(-1),Fr(0)), (Fr(-1),Fr(0))),  # Li4(1/2)
    ((10,86), (-h,Fr(-1)), (-h,Fr(0))),          # Li4(1-1/sqrt2)
    ((-10,86), (h,Fr(-1)), (Fr(0),Fr(-1))),      # Li4(2-sqrt2)
]
SL = {}
for cf, rv, omrv in terms:
    s = li_sym(4, rv, omrv)
    for w, c in s.items():
        cc = cmul2(pair(cf), (c, Fr(0)))
        cur = SL.get(w, (Fr(0),Fr(0)))
        SL[w] = (cur[0]+cc[0], cur[1]+cc[1])
SL = {w: c for w, c in SL.items() if c != (Fr(0),Fr(0))}
print("Li4part symbol terms:", len(SL))
dL = delta22_spec(SL, 2)
print("delta22(Li4part):", dL if dL else "= 0 (as expected)")

# rebuild S(F) at c0 (from previous cell, via module)
import json, re
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
R3 = json.load(open(G1 + "hpl_fit_k3_(0,0,0).json"))
def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, Fr(int(num), int(den))
g = {0: (-h, Fr(0)), 1: (h, Fr(1)), -1: (-h, Fr(1))}
c0pow = {-1: (Fr(0), Fr(1)), 0: (Fr(1), Fr(0)), 1: (Fr(0), h)}
SF = {}
for n, c in R3:
    e, w, cf = parse(n, c)
    if len(w) != 4: continue
    word = tuple(g[a] for a in reversed(w))
    cc = cmul2(c0pow[e], (cf, Fr(0)))
    cur = SF.get(word, (Fr(0), Fr(0)))
    SF[word] = (cur[0] + cc[0], cur[1] + cc[1])
SF = {w: c for w, c in SF.items() if c != (Fr(0), Fr(0))}

# difference
DIFF = {}
for w in set(SF) | set(SL):
    a = SF.get(w, (Fr(0),Fr(0))); b_ = SL.get(w, (Fr(0),Fr(0)))
    dv = (a[0]-b_[0], a[1]-b_[1])
    if dv != (Fr(0),Fr(0)): DIFF[w] = dv
print("difference terms:", len(DIFF))
dD = delta22_spec(DIFF, 2)
print("delta22(DIFF):", dD if dD else "= 0")

# is DIFF a shuffle product? build shuffle subspace over 2 letters and test membership via linear algebra
# shuffles: (1,3): e_i ш word3 ; (2,2): word2 ш word2 ; also (1,1,2) etc are in their span
import itertools
words4 = list(itertools.product(range(2), repeat=4))
def shuf_tuples(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],)+w for w in shuf_tuples(u[1:],v)] + [(v[0],)+w for w in shuf_tuples(u,v[1:])]
rows = []  # each row: dict word4 -> coeff (over Q x Q? use two separate Q problems: p-part and q-part)
for i in range(2):
    for w3 in itertools.product(range(2), repeat=3):
        rows.append({w: 1 for w in shuf_tuples((i,), w3)})
for w2a in itertools.product(range(2), repeat=2):
    for w2b in itertools.product(range(2), repeat=2):
        rows.append({w: 1 for w in shuf_tuples(w2a, w2b)})
# DIFF as vector over words in basis letters: need to EXPAND letter tensors into basis words
def expand(S):
    out = {}
    for w, cf in S.items():
        terms = [((), cf)]
        for v in w:
            terms = [(pref+(j,), cmul2(c, ((v[j]), Fr(0))) if False else None) for pref, c in terms for j in range(2) if v[j]]
        # redo properly below
    return None
def expand2(S):
    out = {}
    for w, cf in S.items():
        terms = [((), cf)]
        for v in w:
            newt = []
            for pref, c in terms:
                for j in range(2):
                    if v[j]:
                        newt.append((pref+(j,), (c[0]*v[j], c[1]*v[j])))
            terms = newt
        for pref, c in terms:
            cur = out.get(pref, (Fr(0),Fr(0)))
            out[pref] = (cur[0]+c[0], cur[1]+c[1])
    return {k: v for k, v in out.items() if v != (Fr(0),Fr(0))}
DE = expand2(DIFF)
print("expanded DIFF terms:", len(DE))
# membership: for p-part and q-part separately, is the vector in row span?
import sympy as sp
def in_span(vec_d):
    allw = sorted(set(list(vec_d.keys()) + [w for r in rows for w in r]))
    idx = {w: i for i, w in enumerate(allw)}
    M = sp.Matrix([[r.get(w, 0) for w in allw] for r in rows])
    v = sp.Matrix([vec_d.get(w, 0) for w in allw])
    sol = M.T.gauss_jordan_solve(v)[0] if M.rank() == sp.Matrix.hstack(M.T, v).rank() else None
    return sol
vp = {w: c[0] for w, c in DE.items()}
vq = {w: c[1] for w, c in DE.items()}
sp_p = in_span(vp); sp_q = in_span(vq)
print("p-part in shuffle span:", sp_p is not None)
print("q-part in shuffle span:", sp_q is not None)
