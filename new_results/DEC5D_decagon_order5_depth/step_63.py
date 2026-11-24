
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

# generators
gens = []  # (name, {word: coeff Fraction})
letters = {(0,): {(0,): QS5(1)}, (1,): {(1,): QS5(1)}}
for i in (0,1):
    for w in itertools.product([0,1], repeat=4):
        gens.append((f"x{i}sh{w}", sh({(i,): QS5(1)}, {w: QS5(1)})))
for u in itertools.product([0,1], repeat=2):
    for v in itertools.product([0,1], repeat=3):
        gens.append((f"{u}sh{v}", sh({u: QS5(1)}, {v: QS5(1)})))
# Li_5 at S-units over <2, phi>
PHI = D.PHI
cands = []
for a in range(-4,5):
    for b in range(-6,7):
        for sgn in (1,-1):
            z = (QS5(2)**a)*(PHI**b)*(QS5(sgn))
            if z in (QS5(0), QS5(1)): continue
            try:
                lv1 = logvec_dec(z); lv2 = logvec_dec(QS5(1)-z)
            except ValueError:
                continue
            cands.append(z)
cands = sorted(set(cands), key=float)
print("S-units with z,1-z over <2,phi>:", len(cands), [str(z) for z in cands])
for z in cands:
    gens.append((f"Li5({z})", li(5, z)))

def to_vec(s):
    v = [F(0)]*32
    for w, c in s.items():
        v[widx[w]] = F(c.p)   # rational part; gens all have rational coeffs
        assert c.q == 0
    return v
G = [to_vec(g) for _, g in gens]
print("generators:", len(G))
# rank of generators
def rank(mat):
    M = [row[:] for row in mat if any(row)]
    if not M: return 0
    nr, nc = len(M), len(M[0])
    r = 0; piv = []
    for c in range(nc):
        p = next((i for i in range(r, nr) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(nr):
            if i != r and M[i][c]:
                f = M[i][c]/M[r][c]
                for j in range(c, nc): M[i][j] -= f*M[r][j]
        piv.append(c); r += 1
        if r == nr: break
    return r
rk = rank(G)
print("rank of classical generator matrix:", rk)
pickle.dump({'gens': [(n, {''.join(map(str,w)): str(c) for w,c in g.items()}) for n,g in gens]}, open("classical_gens.pkl","wb"))
