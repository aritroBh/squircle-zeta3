import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
from oct5d3_basis import *
from fractions import Fraction as F

t0 = time.time()
# target symbol S(W5): load validated block_symbols (will recompute from parsed_mons in verifier)
base2 = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
raw = pickle.load(open(base2+"block_symbols.pkl","rb"))
def parse(s): return {tuple(map(int,w)): Q(F(p), F(q)) for w,(p,q) in s.items()}
W5sym = parse(raw['total'])
T = add(W5sym, scale(C_symbol(), -1))     # classical target = S(W5) - S(C)
print("obstruction of T (must be [0,0]):", obstruction(T))

atoms = atom_list()
vis = [a for a in atoms if a[1]=='vis']
words = sorted(set(T) | set().union(*[set(a[2]) for a in vis]))
print("words:", len(words), "visible atoms:", len(vis))
M = [[a[2].get(w, ZERO) for a in vis] for w in words]
bvec = [T.get(w, ZERO) for w in words]

# RREF over QS2
nr, nc = len(M), len(M[0])
A = [row[:] + [bvec[i]] for i, row in enumerate(M)]
pivots = []; r = 0
for c in range(nc):
    piv = next((i for i in range(r, nr) if A[i][c] != ZERO), None)
    if piv is None: continue
    A[r], A[piv] = A[piv], A[r]
    pv = A[r][c]
    A[r] = [x/pv for x in A[r]]
    for i in range(nr):
        if i != r and A[i][c] != ZERO:
            f = A[i][c]
            A[i] = [A[i][j] - f*A[r][j] for j in range(nc+1)]
    pivots.append(c); r += 1
    if r == nr: break
incons = [i for i in range(nr) if all(A[i][c]==ZERO for c in range(nc)) and A[i][nc]!=ZERO]
print("rank:", r, "of", nr, "; inconsistent rows:", len(incons),
      "=> T in span of M4f visible atoms:", not incons, " (%.1fs)" % (time.time()-t0))
x0 = [ZERO]*nc
for i, c in enumerate(pivots): x0[c] = A[i][nc]
resid = [sum((M[i][c]*x0[c] for c in range(nc)), ZERO) - bvec[i] for i in range(nr)]
print("residual tensor zero:", all(v==ZERO for v in resid))
# kernel basis
free = [c for c in range(nc) if c not in pivots]
print("kernel dim:", len(free))
K = []
for f in free:
    v = [ZERO]*nc; v[f] = ONE
    for i, c in enumerate(pivots): v[c] = -A[i][f]
    K.append(v)
pickle.dump({'x0': x0, 'K': K, 'pivots': pivots, 'free': free,
             'vis_labels': [a[0] for a in vis]}, open("symfit_m4f.pkl","wb"))
print("saved symfit_m4f.pkl;  time %.1fs" % (time.time()-t0))
