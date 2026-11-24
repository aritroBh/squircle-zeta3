# Planner check (independent of cor4_sym.py): recompute the symbol of F = weight-4 part of R3(c)
# directly from hpl_fit_k3_(0,0,0).json, apply delta_{2,2}: a(x)b(x)c(x)d -> (a^b)^(c^d) in L2(L2 V),
# and specialise at c = 1/sqrt2, sqrt3/2, 2/5.
import json, itertools, sympy as sp
from fractions import Fraction
D = json.load(open('<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
c = sp.symbols('c')
# letters: index 0=[c], 1=[1-c], 2=[1+c]; g(0)=+[c], g(1)=-[1-c], g(-1)=+[1+c]
G = {0: (0, 1), 1: (1, -1), -1: (2, 1)}
sym = {}   # (i,j,k,l) -> coefficient (sympy expr in c)
for key, val in D:
    e, w = key.split(' ', 1)
    e = int(e[2:]); w = eval(w[1:])
    if len(w) != 4: continue
    word = list(reversed(w))       # S(H_{a1..an}) = g(a_n) x ... x g(a_1)
    idx = tuple(G[a][0] for a in word); sg = 1
    for a in word: sg *= G[a][1]
    sym[idx] = sym.get(idx, 0) + sg * sp.Rational(val) * c**e
nW = sum(1 for k, v in D if len(eval(k.split(' ',1)[1][1:])) == 4)
print('weight-4 HPL terms:', nW)
# delta22 in basis of L2(L2 V): pairs of pairs (p<q) with p,q in {(0,1),(0,2),(1,2)}
def wedge(i, j):   # i^j in L2 V as dict
    if i == j: return {}
    return {(i, j): 1} if i < j else {(j, i): -1}
def wedge2(A, B):
    out = {}
    for p, x in A.items():
        for q, y in B.items():
            if p == q: continue
            k, s = ((p, q), 1) if p < q else ((q, p), -1)
            out[k] = out.get(k, 0) + s * x * y
    return out
d22 = {}
for (i, j, k, l), co in sym.items():
    for key, v in wedge2(wedge(i, j), wedge(k, l)).items():
        d22[key] = sp.expand(d22.get(key, 0) + v * co)
d22 = {k: v for k, v in d22.items() if v != 0}
print('delta22(F) as function of c:')
for k, v in d22.items(): print('  ', k, sp.factor(v))
# controls: products and Li4 must give 0
def d22_of(symb):
    out = {}
    for (i, j, k, l), co in symb.items():
        for key, v in wedge2(wedge(i, j), wedge(k, l)).items():
            out[key] = out.get(key, 0) + v * co
    return {k: v for k, v in out.items() if v != 0}
def shuffle(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + s for s in shuffle(u[1:], v)] + [(v[0],) + s for s in shuffle(u, v[1:])]
import random; random.seed(1); bad = 0
for t in range(200):
    a = tuple(random.randrange(3) for _ in range(2)); b = tuple(random.randrange(3) for _ in range(2))
    sh = {}
    for s in shuffle(a, b): sh[s] = sh.get(s, 0) + 1
    a1 = (random.randrange(3),); b3 = tuple(random.randrange(3) for _ in range(3))
    sh2 = {}
    for s in shuffle(a1, b3): sh2[s] = sh2.get(s, 0) + 1
    if d22_of(sh) or d22_of(sh2): bad += 1
print('200 random shuffle controls nonzero:', bad, '| Li4(c) control:', d22_of({(1, 0, 0, 0): -1}))
# specialisations: express letters as vectors in a multiplicative basis
def spec(vecs, cval, label):
    n = len(vecs[0])
    def w(i, j):
        out = {}
        for a in range(n):
            for b in range(n):
                x = vecs[i][a] * vecs[j][b]
                if x == 0 or a == b: continue
                k, s = ((a, b), 1) if a < b else ((b, a), -1)
                out[k] = out.get(k, 0) + s * x
        return out
    tot = {}
    for (i, j, k, l), co in sym.items():
        for key, v in wedge2(w(i, j), w(k, l)).items():
            tot[key] = sp.nsimplify(sp.expand(tot.get(key, 0) + v * co.subs(c, cval)))
    tot = {k: sp.simplify(v) for k, v in tot.items() if sp.simplify(v) != 0}
    print(label, '->', tot if tot else 0)
h = sp.Rational(1, 2)
# c = 1/sqrt2, basis ([2],[1+sqrt2]): [c]=-1/2[2], [1-c]=-[1+s2]-1/2[2], [1+c]=[1+s2]-1/2[2]
spec([(-h, 0), (-h, -1), (-h, 1)], 1/sp.sqrt(2), 'c=1/sqrt2 (square corner)')
# c = sqrt3/2, basis ([2],[3],[2+sqrt3]): [c]=1/2[3]-[2], [1-c]=[2-s3]-[2]=-[2+s3]-[2], [1+c]=[2+s3]-[2]
spec([(-1, h, 0), (-1, 0, -1), (-1, 0, 1)], sp.sqrt(3)/2, 'c=sqrt3/2 (equilateral-triangle corner)')
# c = 1/2 (hexagon corner), basis ([2],[3]): [c]=-[2], [1-c]=-[2], [1+c]=[3]-[2]
spec([(-1, 0), (-1, 0), (-1, 1)], h, 'c=1/2 (hexagon corner)')
# c = 2/5, basis ([2],[3],[5],[7])
spec([(1, 0, -1, 0), (0, 1, -1, 0), (0, 0, -1, 1)], sp.Rational(2, 5), 'c=2/5')
