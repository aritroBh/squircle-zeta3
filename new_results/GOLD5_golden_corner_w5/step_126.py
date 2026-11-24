import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
Fraction = Fr
import gold5_lib as G
def nfadd(a, b): return (a[0]+b[0], a[1]+b[1])
def nfmul(a, b): return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def pad(a): return a if len(a) == 2 else (a[0], Fr(0))
Z2 = (Fr(0), Fr(0))
adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
keys = [k for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
V = {k: (Fr(adata[k]['cvec'][0][0]), Fr(adata[k]['cvec'][1][0])) for k in keys}
P = G.golden_point('gol3')
atomSE = {}
for k in keys:
    x, y = eval(adata[k]['xy'][0]), eval(adata[k]['xy'][1])
    atomSE[k] = G.ft_expand_sym_nf(G.li32_sym_nf(P, x, y), 2)
def img32(SE):
    img = {}
    for w, cf in SE.items():
        cf = pad(cf)
        u, v = w[:3], w[3:]
        qc, _ = G.quotient_coords_nf({u: cf})
        if any(any(q) for q in qc):
            a = img.get(v, [Z2, Z2]); a[0] = nfadd(a[0], qc[0]); a[1] = nfadd(a[1], qc[1]); img[v] = a
    return {k: v for k, v in img.items() if any(any(q) for q in v)}
Som = G.ft_expand_sym_nf(P.omega5(), 2)
i32om = img32(Som)
sup32 = sorted(set(i32om) | {w for k in keys for w in img32(atomSE[k])})
print("sup32:", sup32)
# NF-packed vectors: [cvec1, cvec2] + [img32[(v,fc)] for v in sup32 for fc in (0,1)]
def nfvec_of(k):
    i32 = img32(atomSE[k])
    out = [(V[k][0], Fr(0)), (V[k][1], Fr(0))]
    for w in sup32:
        a = i32.get(w, [Z2, Z2])
        out += [a[0], a[1]]
    return tuple(out)
atomnf = {k: nfvec_of(k) for k in keys}
Tom = ((Fr(0), Fr(0)), (Fr(972), Fr(432))) + tuple()
i32omv = [i32om.get(w, [Z2, Z2]) for w in sup32]
# NOTE: target image is for Omega - N5 == shuffle-free: img(Omega) - img(N5) must have ZERO Lie3 image,
# i.e. img(N5) == img(Omega). And cobracket(N5) == cobracket(Omega) = ((0,0),(972,432)).
T = ((Fr(0), Fr(0)), (Fr(972), Fr(432))) + tuple(q for a in i32omv for q in a)
dimF = len(T)
def qvec(nfv):
    out = []
    for a in nfv: out += [a[0], a[1]]
    return tuple(out)
Tq = qvec(T)
atomq = {k: qvec(v) for k, v in atomnf.items()}
SQ5v = (Fr(0), Fr(1))
def cols_of(k):
    vv = atomnf[k]
    c1 = qvec(vv)
    c2 = qvec(tuple(nfmul(SQ5v, a) for a in vv))
    return [c1, c2]
def solve_lin(cols, T):
    n, m = len(T), len(cols)
    A = [[cols[c][r] for c in range(m)] + [T[r]] for r in range(n)]
    R = 0; piv = []
    for c in range(m):
        p = next((r for r in range(R, n) if A[r][c]), None)
        if p is None: continue
        A[R], A[p] = A[p], A[R]
        q = A[R][c]; A[R] = [x / q for x in A[R]]
        for r in range(n):
            if r != R and A[r][c]:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[R])]
        piv.append(c); R += 1
        if R == n: break
    for r in range(R, n):
        if A[r][m]: return None
    part = [Fr(0)] * m
    for i, c in enumerate(piv): part[c] = A[i][m]
    return part
def height(sol):
    return max(max(abs(x.numerator), x.denominator) for x in sol if x)
best = []
for k1, k2, k3 in itertools.combinations(keys, 3):
    cols = []
    for k in (k1, k2, k3): cols += cols_of(k)
    s = solve_lin(cols, Tq)
    if s is None: continue
    # verify
    acc = [Fr(0)]*len(Tq)
    for j, c in enumerate(s):
        if c:
            for rr in range(len(Tq)): acc[rr] += c*cols[j][rr]
    if not all(a-b == 0 for a, b in zip(acc, Tq)): continue
    best.append((height(s), k1, k2, k3, s))
best.sort(key=lambda t: t[0])
print(f"{len(best)} gol3 TRIPLE solutions with cobracket + zero-Lie3-image:")
for h, k1, k2, k3, s in best[:8]:
    print(f"  h={h}:")
    for i, k in enumerate((k1, k2, k3)):
        d = (s[2*i], s[2*i+1])
        if any(d): print(f"    ({G.nf_str(d)})*Li32{tuple(adata[k]['xy'])}")
