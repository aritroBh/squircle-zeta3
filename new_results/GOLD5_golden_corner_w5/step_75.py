import sys
sys.path.insert(0, '<project>/research/08_theory/GOLD5_golden_corner_w5')
from fractions import Fraction as Fr
import gold5_lib as G

def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def nfadd(a, b):
    return (a[0]+b[0], a[1]+b[1])

PHI3 = (Fr(-2), Fr(1))     # phi^-3 = sqrt5-2
PHI6 = (Fr(9), Fr(-4))     # phi^-6 = 9-4sqrt5
PH2  = (Fr(1,4), Fr(1,4))  # phi/2 = (1+sqrt5)/4
TPHI2 = (Fr(3), Fr(-1))    # 3-sqrt5 = 2 phi^-2

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
    # cobracket of N5
    Dtot = {}
    oQ = (Fr(0), Fr(0)); oT = (Fr(0), Fr(0))
    for d, x, y in N5[tag]:
        dw = G.delta_word_sym(P, G.li32_word(P, x, y), F.zero(), F.one())
        for lw, v in dw.items():
            Dtot[lw] = nfadd(Dtot.get(lw, (Fr(0), Fr(0))), nfmul(d, v))
        S_at = G.ft_expand_sym_nf(G.li32_sym_nf(P, x, y), 2)
        obs = G.obs_of_SE(S_at)
        oQ = nfadd(oQ, nfmul(d, obs[0] if obs[0] else (Fr(0), Fr(0))))
        oT = nfadd(oT, nfmul(d, obs[1] if obs[1] else (Fr(0), Fr(0))))
    qc, free = G.quotient_coords_nf(G.ft_expand_sym_nf(Dtot, 2))
    print(f"{tag}: N5 quotient cobracket = {[G.nf_str(q) for q in qc]} "
          f"(target {[G.nf_str(q) for q in delta_target[tag]]})")
    print(f"      N5 obs = [{G.nf_str(oQ)}, {G.nf_str(oT)}] "
          f"(target [{G.nf_str(obs_target[tag][0])}, {G.nf_str(obs_target[tag][1])}])")
    assert tuple(qc) == delta_target[tag] and oQ == obs_target[tag][0] and oT == obs_target[tag][1]
    # class of S(N5)
    SN5 = {}
    for d, x, y in N5[tag]:
        for w, cf in G.li32_sym_nf(P, x, y).items():
            SN5[w] = nfadd(SN5.get(w, (Fr(0), Fr(0))), nfmul(d, cf if isinstance(cf, tuple) else (cf, Fr(0))))
    clN5 = G.quotient5(G.ft_expand_sym_nf(SN5, 2))
    print(f"      class(S(N5)) = {[G.nf_str(q) for q in clN5]}")
print("N5 verification PASS (both corners)")
