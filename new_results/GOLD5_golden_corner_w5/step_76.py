import sys
sys.path.insert(0, '<project>/research/08_theory/GOLD5_golden_corner_w5')
from fractions import Fraction as Fr
import gold5_lib as G

def pad(a):
    return a if len(a) == 2 else (a[0], Fr(0))
def nfmul(a, b):
    a, b = pad(a), pad(b)
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def nfadd(a, b):
    a, b = pad(a), pad(b)
    return (a[0]+b[0], a[1]+b[1])
Z2 = (Fr(0), Fr(0))

PHI3 = (Fr(-2), Fr(1)); PHI6 = (Fr(9), Fr(-4))
PH2  = (Fr(1,4), Fr(1,4)); TPHI2 = (Fr(3), Fr(-1))
N5 = {
 'gol1': [((Fr(2), Fr(-1,2)), PHI3, PH2), ((Fr(2), Fr(3,2)), TPHI2, PH2)],
 'gol3': [((Fr(-9,8), Fr(-1,2)), PHI6, (Fr(1), Fr(0)))],
}
obs_target = {'gol1': ((Fr(72), Fr(0)), (Fr(72), Fr(-54))),
              'gol3': ((Fr(0), Fr(0)), (Fr(972), Fr(432)))}
delta_target = {'gol1': ((Fr(-72), Fr(0)), (Fr(72), Fr(-54))),
                'gol3': ((Fr(0), Fr(0)), (Fr(972), Fr(432)))}
for tag in ('gol1', 'gol3'):
    P = G.golden_point(tag)
    F = P.F
    Dtot = {}
    oQ = Z2; oT = Z2
    for d, x, y in N5[tag]:
        dw = G.delta_word_sym(P, G.li32_word(P, x, y), F.zero(), F.one())
        for lw, v in dw.items():
            Dtot[lw] = nfadd(Dtot.get(lw, Z2), nfmul(d, v))
        S_at = G.ft_expand_sym_nf(G.li32_sym_nf(P, x, y), 2)
        obs = G.obs_of_SE(S_at)
        oQ = nfadd(oQ, nfmul(d, obs[0] if obs[0] else Z2))
        oT = nfadd(oT, nfmul(d, obs[1] if obs[1] else Z2))
    qc, free = G.quotient_coords_nf(G.ft_expand_sym_nf(Dtot, 2))
    print(f"{tag}: N5 quotient cobracket = {[G.nf_str(q) for q in qc]} (target {[G.nf_str(q) for q in delta_target[tag]]})")
    print(f"      N5 obs = [{G.nf_str(oQ)}, {G.nf_str(oT)}] (target [{G.nf_str(obs_target[tag][0])}, {G.nf_str(obs_target[tag][1])}])")
    assert tuple(pad(q) for q in qc) == delta_target[tag] and oQ == obs_target[tag][0] and oT == obs_target[tag][1]
    SN5 = {}
    for d, x, y in N5[tag]:
        for w, cf in G.li32_sym_nf(P, x, y).items():
            SN5[w] = nfadd(SN5.get(w, Z2), nfmul(d, cf if isinstance(cf, tuple) else (cf,)))
    clN5 = G.quotient5(G.ft_expand_sym_nf(SN5, 2))
    print(f"      class(S(N5)) = {[G.nf_str(q) for q in clN5]}")
print("N5 verification PASS (both corners)")
