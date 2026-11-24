
import pickle
from fractions import Fraction as F
import itertools
import dec5d_lib as D
from qsqrt2 import QS2
ZERO2 = QS2(0)
kit2 = D.make_kit(QS2, D.logvec2, ('log2','log(1+sqrt2)'))
sh2, li2 = kit2['shuffle'], kit2['li_symbol']
words5 = list(itertools.product([0,1], repeat=5))
widx = {w:i for i,w in enumerate(words5)}
gens2 = []
for i in (0,1):
    for w in itertools.product([0,1], repeat=4):
        gens2.append((f"x{i}ш{w}", sh2({(i,): QS2(1)}, {w: QS2(1)})))
for u in itertools.product([0,1], repeat=2):
    for v in itertools.product([0,1], repeat=3):
        gens2.append((f"{u}ш{v}", sh2({u: QS2(1)}, {v: QS2(1)})))
from qsqrt2 import SQ2, A as A2
cands2 = []
for a in range(-6,7):
    for b in range(-5,6):
        for sgn in (1,-1):
            z = (SQ2**a)*(A2**b)*(QS2(sgn))
            if z in (QS2(0), QS2(1)): continue
            try: D.logvec2(z); D.logvec2(QS2(1)-z)
            except ValueError: continue
            cands2.append(z)
print("octagon S-units z,1-z over <2,1+sqrt2>:", len(set(cands2)))
for z in sorted(set(cands2), key=float):
    gens2.append((f"Li5({z})", li2(5, z)))
def to_vec2(s):
    v = [F(0)]*32
    for w, c in s.items():
        assert c.q == 0
        v[widx[w]] = F(c.p)
    return v
G2 = [to_vec2(g) for _, g in gens2]
def solve2(G, t):
    M = [[G[j][i] for j in range(len(G))] + [t[i]] for i in range(32)]
    R = [r[:] for r in M]
    pivrow = {}; r = 0
    for c in range(len(G)):
        p = next((i for i in range(r, 32) if R[i][c] != 0), None)
        if p is None: continue
        R[r], R[p] = R[p], R[r]
        inv = 1/R[r][c]
        R[r] = [x*inv for x in R[r]]
        for i in range(32):
            if i != r and R[i][c] != 0:
                f = R[i][c]
                R[i] = [R[i][j] - f*R[r][j] for j in range(len(G)+1)]
        pivrow[c] = r; r += 1
        if r == 32: break
    for i in range(r, 32):
        if R[i][len(G)] != 0: return None
    sol = [F(0)]*len(G)
    for c, i in pivrow.items(): sol[c] = R[i][len(G)]
    for i in range(32):
        if sum(sol[j]*G[j][i] for j in range(len(G))) != t[i]: return "BAD"
    return sol
sub_sym = pickle.load(open("octagon_control_syms.pkl","rb"))['sub_sym']
tot2 = kit2['add'](*sub_sym.values())
tp = [F(0)]*32; tq = [F(0)]*32
for w, c in tot2.items():
    tp[widx[w]] = F(c.p); tq[widx[w]] = F(c.q)
sp = solve2(G2, tp); sq = solve2(G2, tq)
print("octagon total in classical span: p:", sp, " q:", sq)
print("(expected: NOT in span -- the octagon needs depth two)")
# rank of octagon classical span
def rank(G):
    M = [list(col) for col in zip(*G)]
    r = 0
    for c in range(len(G)):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        inv = 1/M[r][c]
        M[r] = [x*inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [M[i][j]-f*M[r][j] for j in range(len(G))]
        r += 1
        if r == len(M): break
    return r
print("octagon classical span rank:", rank(G2), " decagon rank was 30")
