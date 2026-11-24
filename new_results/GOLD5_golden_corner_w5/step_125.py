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
P = G.golden_point('gol3')
# per-atom: cvec (2 Fr), img32/img23 as Q-vectors over the union of support words
atomSE = {}
for k in keys:
    x, y = eval(adata[k]['xy'][0]), eval(adata[k]['xy'][1])
    atomSE[k] = G.ft_expand_sym_nf(G.li32_sym_nf(P, x, y), 2)
# support words for images: all (v2 for 32-split) and (u2 for 23-split) over Omega and atoms
def imgs(SE):
    img32 = {}; img23 = {}
    for w, cf in SE.items():
        cf = pad(cf)
        u, v = w[:3], w[3:]
        qc, _ = G.quotient_coords_nf({u: cf})
        if any(any(q) for q in qc):
            a = img32.get(v, [Z2, Z2]); a[0] = nfadd(a[0], qc[0]); a[1] = nfadd(a[1], qc[1]); img32[v] = a
        u, v = w[:2], w[3:]
        qc, _ = G.quotient_coords_nf({v: cf})
        if any(any(q) for q in qc):
            a = img23.get(u, [Z2, Z2]); a[0] = nfadd(a[0], qc[0]); a[1] = nfadd(a[1], qc[1]); img23[u] = a
    return img32, img23
Som = G.ft_expand_sym_nf(P.omega5(), 2)
i32om, i23om = imgs(Som)
sup32 = sorted(set(i32om) | {k2 for k in keys for k2 in imgs(atomSE[k])[0]})
sup23 = sorted(set(i23om) | {k2 for k in keys for k2 in imgs(atomSE[k])[1]})
print(f"support: {len(sup32)} (3,2)-words, {len(sup23)} (2,3)-words")
def vec_of(img, sup):
    out = []
    for w in sup:
        a = img.get(w, [Z2, Z2])
        out += [a[0][0], a[0][1], a[1][0], a[1][1]]
    return tuple(out)
atomvec = {}
for k in keys:
    i32, i23 = imgs(atomSE[k])
    atomvec[k] = (V[k][0], V[k][1]) + vec_of(i32, sup32) + vec_of(i23, sup23)
target_om = (Fr(0), Fr(0), Fr(972), Fr(432)) + (Fr(0),)*(4*len(sup32)) + (Fr(0),)*(4*len(sup23))
dim = len(target_om)
print("constraint dim:", dim)
SQ5v = (Fr(0), Fr(1))
def cols_of(k):
    vv = atomvec[k]
    return [tuple(vv), tuple(nfmul(SQ5v, (x, Fr(0)))[0] * 0 + x for x in ())]
def cols_of(k):
    vv = atomvec[k]
    c1 = tuple(vv)
    c2 = tuple(5 * vv[i+1] if i % 2 == 0 else vv[i-1] for i in range(len(vv)))
    # careful: multiplication by sqrt5 on Q^2 pairs: (a,b) -> (5b, a); but cvec entries are rational (b=0)
    c2 = []
    for i, x in enumerate(vv):
        c2.append(x)  # placeholder, fixed below
    # pairs are (rational, 0): sqrt5*(r,0) = (0, r): swap within pairs of the img part; for cvec part (2 rationals): also pairs
    c2 = []
    # cvec part: entries 0..1 are single rationals treated as (r,0): sqrt5 path -> (0,r): but our vector packs them singly...
    return [c1, c2]
# SIMPLER: pack everything as NF pairs. cvec entries are rational -> pairs (r,0). img entries already pairs.
def nfvec_of(k):
    i32, i23 = imgs(atomvec and atomSE[k] or atomSE[k])
    v1, v2 = V[k]
    return ((v1, Fr(0)), (v2, Fr(0))) + tuple((a[0], a[1]) for w in sup32 for a in [i32.get(w, [Z2, Z2])] for a in a) + \
           tuple((a[0], a[1]) for w in sup23 for a in [i23.get(w, [Z2, Z2])] for a in a)
atomnf = {k: nfvec_of(k) for k in keys}
# target: cobracket pairs rational, images zero
T = ((Fr(0), Fr(0)), (Fr(0), Fr(0)), (Fr(972), Fr(432)), (Fr(0), Fr(0))) + tuple()
T = ((Fr(0), Fr(0)), (Fr(0), Fr(0)), (Fr(972), Fr(432))) + tuple(Z2 for _ in range(2*len(sup32)*2 + 2*len(sup23)*2 - 1))
print("len check:", len(atomnf[keys[0]]), len(T))
