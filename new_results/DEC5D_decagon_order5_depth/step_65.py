
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
cands = sorted(set(cands), key=float)
for z in cands:
    gens.append((f"Li5({z})", li(5, z)))
def to_vec_f(s):
    v = [F(0)]*32
    for w, c in s.items():
        assert c.q == 0
        v[widx[w]] = F(c.p)
    return v
G = [to_vec_f(g) for _, g in gens]
def solve(G, t):
    M = [list(row) + [tt] for row, tt in zip(zip(*G), t)]
    nr, nc = len(M), len(G)
    r = 0; piv = []
    for c in range(nc):
        p = next((i for i in range(r, nr) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(nr):
            if i != r and M[i][c]:
                f = M[i][c]/M[r][c]
                for j in range(c, nc+1): M[i][j] -= f*M[r][j]
        piv.append(c); r += 1
    for i in range(r, nr):
        if M[i][nc] != 0: return None
    sol = [F(0)]*nc
    for i, c in enumerate(piv): sol[c] = M[i][nc]
    return sol
def remap(s):
    return {tuple(0 if x==0 else 1 for x in w): c for w, c in s.items()}
data = pickle.load(open("dec_block_symbols.pkl","rb"))
tot = remap(data['tot'])
tp = [F(0)]*32; tq = [F(0)]*32
for w, c in tot.items():
    tp[widx[w]] = F(c.p); tq[widx[w]] = F(c.q)
solp = solve(G, tp); solq = solve(G, tq)
print("decagon total symbol in classical span: p-part:", solp is not None, " q-part:", solq is not None)
if solp is not None and solq is not None:
    used = [(gens[i][0], solp[i], solq[i]) for i in range(len(gens)) if solp[i] or solq[i]]
    print("solution uses %d generators" % len(used))
    for n, a, b in used: print("   %-44s (%s) + (%s)*sqrt5" % (n, a, b))
    pickle.dump({'used': used}, open("dec_classical_fit.pkl","wb"))
