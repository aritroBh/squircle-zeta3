import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
from cor6_sym import sshuffle, rref
import mpmath as mp
mp.mp.dps = 30
P = G.golden_point('gol1')
F = P.F
# weight-4 quotient T^4/D^4 over 2 letters, NF coefficients
B4 = list(itertools.product(range(2), repeat=4))
def d4_rows():
    rows = []
    for lu in range(1, 4):
        for u in itertools.product(range(2), repeat=lu):
            for v in itertools.product(range(2), repeat=4 - lu):
                rows.append(sshuffle({u: Fr(1)}, {v: Fr(1)}))
    return rows
Dmat4, piv4 = rref([[s.get(w, Fr(0)) for w in B4] for s in d4_rows()])
free4 = [j for j in range(len(B4)) if j not in piv4]
print(f"T^4/D^4 dim {len(free4)} (expect 3), free {[B4[j] for j in free4]}")
def quotient4(SE, dimF=2):
    def _pad(a):
        return tuple(a) + (Fr(0),) * (dimF - len(a))
    v = [_pad(SE.get(w, ())) for w in B4]
    for i, pj in enumerate(piv4):
        if any(v[pj]):
            q = v[pj]
            v = [tuple(vv[t] - q[t] * Dmat4[i][j] for t in range(dimF)) for j, vv in enumerate(v)]
    for j in range(len(B4)):
        if j not in free4:
            assert not any(v[j])
    return tuple(v[j] for j in free4)
def li4_sym(z):
    vz = P.vec(z); v1z = P.vec(F.sub(F.one(), z))
    return {(v1z, vz, vz, vz): (Fr(-1), Fr(0))}
li5cl = {eval(k): None for k in json.load(open(GOLD + 'li5_classes_golden.json'))}
def fval(z):
    return float(mp.mpf(z[0].numerator)/z[0].denominator + mp.mpf(z[1].numerator)/z[1].denominator*mp.sqrt(5))
cl4 = {}
for z in li5cl:
    cl = quotient4(G.ft_expand_sym_nf(li4_sym(z), 2))
    cl4[z] = cl
for z in sorted(cl4, key=fval):
    print(f"  z={str(z):24s} ~ {fval(z):+.6f}: {[G.nf_str(q) for q in cl4[z]]}")
# rank over Q
def q6(v):
    out = []
    for q in v: out += [q[0], q[1]]
    return tuple(out)
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
SQ5v = (Fr(0), Fr(1))
cols = []
for z, cl in cl4.items():
    cols.append(q6(cl)); cols.append(q6(tuple(nfmul(SQ5v, q) for q in cl)))
Mrref, piv = rref([list(r) for r in zip(*cols)])
print("Li4 class span rank over Q:", len(piv))
