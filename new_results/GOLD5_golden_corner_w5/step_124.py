import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
Fraction = Fr
import gold5_lib as G
def nfadd(a, b): return (a[0]+b[0], a[1]+b[1])
def nfmul(a, b): return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def nfsub(a, b): return (a[0]-b[0], a[1]-b[1])
def pad(a): return a if len(a) == 2 else (a[0], Fr(0))
Z2 = (Fr(0), Fr(0))
adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
keys = [k for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
V = {k: (Fr(adata[k]['cvec'][0][0]), Fr(adata[k]['cvec'][1][0])) for k in keys}
def solve4(cols, T):
    A = [[cols[c][r] for c in range(4)] + [T[r]] for r in range(4)]
    for col in range(4):
        p = next((r for r in range(col, 4) if A[r][col]), None)
        if p is None: return None
        A[col], A[p] = A[p], A[col]
        q = A[col][col]; A[col] = [x / q for x in A[col]]
        for r in range(4):
            if r != col and A[r][col]:
                f = A[r][col]
                A[r] = [a - f * b for a, b in zip(A[r], A[col])]
    return [A[i][4] for i in range(4)]
def cols_of(k):
    v1, v2 = V[k]
    return [(v1, Fr(0), v2, Fr(0)), (Fr(0), v1, Fr(0), v2)]
P = G.golden_point('gol3')
atomSE = {}
for k in keys:
    x, y = eval(adata[k]['xy'][0]), eval(adata[k]['xy'][1])
    atomSE[k] = G.ft_expand_sym_nf(G.li32_sym_nf(P, x, y), 2)
# Omega_5 symbol for gol3
Som = G.ft_expand_sym_nf(P.omega5(), 2)
def lie3_img(SE):
    """project deconcat_32 onto (T^3/D^3) x T^2 and deconcat_23 onto T^2 x (T^3/D^3)."""
    img32 = {}
    img23 = {}
    for w, cf in SE.items():
        for (u, v), D in ((w[:3], w[3:]), None), ((w[:2], w[3:]), 'x'):
            pass
        u, v = w[:3], w[3:]
        qc, _ = G.quotient_coords_nf({u: cf})
        if any(any(q) for q in qc):
            a = img32.get(v, [Z2, Z2]); a[0] = nfadd(a[0], qc[0]); a[1] = nfadd(a[1], qc[1]); img32[v] = a
        u, v = w[:2], w[3:]
        qc, _ = G.quotient_coords_nf({v: cf})
        if any(any(q) for q in qc):
            a = img23.get(u, [Z2, Z2]); a[0] = nfadd(a[0], qc[0]); a[1] = nfadd(a[1], qc[1]); img23[u] = a
    return ({k: v for k, v in img32.items() if any(any(q) for q in v)},
            {k: v for k, v in img23.items() if any(any(q) for q in v)})
im_om = lie3_img(Som)
print("Omega_5 gol3: img32 components:", {str(k): [[G.nf_str(q) for q in vv] for k, vv in v.items()] for k, v in im_om[0].items()} if False else len(im_om[0]), "img23:", len(im_om[1]))
# scan N5 pairs for zero image
T3 = (Fr(0), Fr(0), Fr(972), Fr(432))
best = []
for k1, k2 in itertools.combinations(keys, 2):
    s = solve4(cols_of(k1) + cols_of(k2), T3)
    if s is None: continue
    d1 = (s[0], s[1]); d2 = (s[2], s[3])
    SN = {}
    for d, k in ((d1, k1), (d2, k2)):
        for w, cf in atomSE[k].items():
            SN[w] = nfadd(SN.get(w, Z2), nfmul(d, pad(cf)))
    Sres = {w: nfsub(pad(cf), pad(SN.get(w, Z2))) for w, cf in Som.items()}
    Sres = {w: cf for w, cf in Sres.items() if any(cf)}
    i32, i23 = lie3_img(Sres)
    bad = len(i32) + len(i23)
    hn5 = max(max(abs(x.numerator), x.denominator) for x in s if x)
    if bad == 0:
        best.append((hn5, k1, k2, (d1, d2)))
best.sort()
print(f"{len(best)} gol3 N5-pairs with ZERO Lie3 image in both splits:")
for hn5, k1, k2, (d1, d2) in best[:10]:
    print(f"  N5h={hn5}: ({G.nf_str(d1)})*Li32{tuple(adata[k1]['xy'])} + ({G.nf_str(d2)})*Li32{tuple(adata[k2]['xy'])}")
