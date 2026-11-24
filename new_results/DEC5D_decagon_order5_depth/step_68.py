
import pickle
from fractions import Fraction as F
import itertools
import dec5d_lib as D
QS5 = D.QS5; ZERO = QS5(0)
def logvec_dec(z):
    lv = D.logvec5(z)
    if lv[1] != 0: raise ValueError("has logsqrt5 part")
    return (lv[0], lv[2])
kit_dec = D.make_kit(QS5, logvec_dec, ('log2','logphi'))
sh, li = kit_dec['shuffle'], kit_dec['li_symbol']
words5 = list(itertools.product([0,1], repeat=5))
widx = {w:i for i,w in enumerate(words5)}
gens = []
for i in (0,1):
    for w in itertools.product([0,1], repeat=4):
        gens.append((f"x{i}sh{w}", sh({(i,): QS5(1)}, {w: QS5(1)})))
for u in itertools.product([0,1], repeat=2):
    for v in itertools.product([0,1], repeat=3):
        gens.append((f"{u}sh{v}", sh({u: QS5(1)}, {v: QS5(1)})))
PHI = D.PHI
cands = []
for a in range(-4,5):
    for b in range(-6,7):
        for sgn in (1,-1):
            z = (QS5(2)**a)*(PHI**b)*(QS5(sgn))
            if z in (QS5(0), QS5(1)): continue
            try: logvec_dec(z); logvec_dec(QS5(1)-z)
            except ValueError: continue
            cands.append(z)
for z in sorted(set(cands), key=float):
    gens.append((f"Li5({z})", li(5, z)))
def to_vec_f(s):
    v = [F(0)]*32
    for w, c in s.items():
        assert c.q == 0
        v[widx[w]] = F(c.p)
    return v
G = [to_vec_f(g) for _, g in gens]
data = pickle.load(open("dec_block_symbols.pkl","rb"))
def remap(s):
    return {tuple(0 if x==0 else 1 for x in w): c for w, c in s.items()}
tot = remap(data['tot'])
tp = [F(0)]*32; tq = [F(0)]*32
for w, c in tot.items():
    tp[widx[w]] = F(c.p); tq[widx[w]] = F(c.q)

def solve2(G, t):
    # column-style rref, with verification
    import copy
    A = [list(col) for col in G]  # ngens columns of length 32
    b = list(t)
    # solve A c = b via row reduction on augmented matrix (32 rows x (ngens+1))
    M = [[A[j][i] for j in range(len(A))] + [b[i]] for i in range(32)]
    R = [r[:] for r in M]
    nrow, ncol = 32, len(A)
    pivrow = {}; r = 0
    for c in range(ncol):
        p = next((i for i in range(r, nrow) if R[i][c] != 0), None)
        if p is None: continue
        R[r], R[p] = R[p], R[r]
        inv = 1/R[r][c]
        R[r] = [x*inv for x in R[r]]
        for i in range(nrow):
            if i != r and R[i][c] != 0:
                f = R[i][c]
                R[i] = [R[i][j] - f*R[r][j] for j in range(ncol+1)]
        pivrow[c] = r; r += 1
        if r == nrow: break
    # check consistency
    for i in range(r, nrow):
        if R[i][ncol] != 0: return None
    sol = [F(0)]*ncol
    for c, i in pivrow.items(): sol[c] = R[i][ncol]
    return sol
sp = solve2(G, tp)
print("solve2 p-part:", sp is not None)
if sp:
    rp = [sum(sp[j]*G[j][i] for j in range(len(gens))) for i in range(32)]
    print("  residual p:", max(abs(rp[i]-tp[i]) for i in range(32)))
sq = solve2(G, tq)
print("solve2 q-part:", sq is not None)
if sq:
    rq = [sum(sq[j]*G[j][i] for j in range(len(gens))) for i in range(32)]
    print("  residual q:", max(abs(rq[i]-tq[i]) for i in range(32)))
