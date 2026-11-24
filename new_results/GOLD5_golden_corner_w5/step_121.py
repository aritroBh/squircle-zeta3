import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G

def nfadd(a, b): return (a[0]+b[0], a[1]+b[1])
def nfmul(a, b): return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def nfsub(a, b): return (a[0]-b[0], a[1]-b[1])
def pad(a): return a if len(a) == 2 else (a[0], Fr(0))
Z2 = (Fr(0), Fr(0))
MRHO  = (Fr(1,2), Fr(-1,2)); MPH2  = (Fr(1,4), Fr(-1,4))
PH2   = (Fr(1,4), Fr(1,4)); TP2   = (Fr(3), Fr(-1)); TP2O2 = (Fr(3,4), Fr(-1,4))
RHO3 = (Fr(-2), Fr(1)); RHO6 = (Fr(9), Fr(-4)); M1 = (Fr(-1), Fr(0)); ONE = (Fr(1), Fr(0))
N5 = {'gol1': [((Fr(-8, 3), Fr(0)), MRHO, MPH2), ((Fr(8, 27), Fr(2, 3)), M1, RHO3)],
      'gol3': [((Fr(-9, 8), Fr(-1, 2)), RHO6, ONE)]}
P5 = {'gol1': [((Fr(25, 12), Fr(-37, 12)), TP2O2), ((Fr(-17), Fr(9)), MPH2),
               ((Fr(16, 3), Fr(-4, 3)), PH2), ((Fr(-19, 12), Fr(-3, 4)), TP2)],
      'gol3': [((Fr(11, 6), Fr(-1, 2)), MPH2), ((Fr(17, 2), Fr(-5, 2)), PH2),
               ((Fr(20, 9), Fr(-16, 27)), RHO3), ((Fr(-26, 3), Fr(2)), TP2)]}
for tag in ('gol1', 'gol3'):
    P = G.golden_point(tag)
    S = {}
    for w, cf in P.omega5().items():
        S[w] = nfadd(S.get(w, Z2), pad(cf))
    for d, x, y in N5[tag]:
        for w, cf in G.li32_sym_nf(P, x, y).items():
            S[w] = nfsub(S.get(w, Z2), nfmul(d, pad(cf if isinstance(cf, tuple) else (cf,))))
    for c, z in P5[tag]:
        for w, cf in G.li5_sym_nf(P, z).items():
            S[w] = nfsub(S.get(w, Z2), nfmul(c, pad(cf)))
    S = {k: v for k, v in S.items() if any(v)}
    SE = G.ft_expand_sym_nf(S, 2)
    # deconcat 3+2, project T^3 -> T^3/D^3 (quotient coords), keep T^2 raw
    img = {}
    for w, cf in SE.items():
        u, v = w[:3], w[3:]
        qc, _ = G.quotient_coords_nf({u: cf})
        if any(any(q) for q in qc):
            img.setdefault(v, [Z2, Z2])
            img[v][0] = nfadd(img[v][0], qc[0])
            img[v][1] = nfadd(img[v][1], qc[1])
    img = {k: v for k, v in img.items() if any(any(q) for q in v)}
    print(f"{tag}: remainder symbol in D^5 with {len(SE)} words; Lie_3 x T^2 image has {len(img)} components")
    for k, v in sorted(img.items()):
        print(f"   v={k}: ({G.nf_str(v[0])}), ({G.nf_str(v[1])})")
