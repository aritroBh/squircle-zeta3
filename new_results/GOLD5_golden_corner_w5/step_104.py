import sys, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G

def nfadd(a, b): return (a[0]+b[0], a[1]+b[1])
def nfmul(a, b): return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def pad(a): return a if len(a) == 2 else (a[0], Fr(0))
Z2 = (Fr(0), Fr(0))

MRHO  = (Fr(1,2), Fr(-1,2))     # -1/phi
MR3O4 = (Fr(1,2), Fr(-1,4))     # (2-sqrt5)/4 = -rho^3/4
M2PHI = (Fr(-1), Fr(-1))        # -1-sqrt5 = -2 phi
RHO3 = (Fr(-2), Fr(1)); RHO6 = (Fr(9), Fr(-4)); M1 = (Fr(-1), Fr(0)); ONE = (Fr(1), Fr(0))
PH2 = (Fr(1,4), Fr(1,4)); MPH2 = (Fr(1,4), Fr(-1,4)); TP2 = (Fr(3), Fr(-1)); TP2O2 = (Fr(3,4), Fr(-1,4))

N5 = {
 'gol1': [((Fr(-8, 3), Fr(0)), MRHO, MPH2), ((Fr(8, 27), Fr(2, 3)), M1, RHO3)],
 'gol3': [((Fr(-9, 8), Fr(-1, 2)), RHO6, ONE)],
}
P5 = {
 'gol1': [((Fr(25, 12), Fr(-37, 12)), TP2O2),
          ((Fr(-17), Fr(9)), MPH2),
          ((Fr(16, 3), Fr(-4, 3)), PH2),
          ((Fr(-19, 12), Fr(-3, 4)), TP2)],
 'gol3': [((Fr(11, 6), Fr(-1, 2)), MPH2),
          ((Fr(17, 2), Fr(-5, 2)), PH2),
          ((Fr(20, 9), Fr(-16, 27)), RHO3),
          ((Fr(-26, 3), Fr(2)), TP2)],
}
for tag in ('gol1', 'gol3'):
    P = G.golden_point(tag)
    F = P.F
    S5 = P.omega5()
    Delta = {}
    for w, coef in G.signed_words(P):
        dw = G.delta_word_sym(P, w, F.zero(), P.c0)
        for lw, v in dw.items():
            Delta[lw] = nfadd(Delta.get(lw, Z2), nfmul(coef, pad(v)))
    qcD, _ = G.quotient_coords_nf(G.ft_expand_sym_nf({k: v for k, v in Delta.items() if any(v)}, 2))
    # N5
    Dtot = {}; oQ = Z2; oT = Z2
    for d, x, y in N5[tag]:
        dw = G.delta_word_sym(P, G.li32_word(P, x, y), F.zero(), F.one())
        for lw, v in dw.items():
            Dtot[lw] = nfadd(Dtot.get(lw, Z2), nfmul(d, pad(v)))
        S_at = G.ft_expand_sym_nf(G.li32_sym_nf(P, x, y), 2)
        obs = G.obs_of_SE(S_at)
        oQ = nfadd(oQ, nfmul(d, pad(obs[0]) if obs[0] else Z2))
        oT = nfadd(oT, nfmul(d, pad(obs[1]) if obs[1] else Z2))
    qcN, _ = G.quotient_coords_nf(G.ft_expand_sym_nf(Dtot, 2))
    obs_om = tuple(pad(o) for o in G.obs_of_SE(G.ft_expand_sym_nf(S5, 2)))
    print(f"{tag}: cobracket match {tuple(pad(q) for q in qcN) == tuple(pad(q) for q in qcD)}, "
          f"obs match {(oQ, oT) == obs_om}")
    assert tuple(pad(q) for q in qcN) == tuple(pad(q) for q in qcD) and (oQ, oT) == obs_om
    # classes
    SN5 = {}
    for d, x, y in N5[tag]:
        for w, cf in G.li32_sym_nf(P, x, y).items():
            SN5[w] = nfadd(SN5.get(w, Z2), nfmul(d, pad(cf if isinstance(cf, tuple) else (cf,))))
    SP5 = {}
    for c, z in P5[tag]:
        for w, cf in G.li5_sym_nf(P, z).items():
            SP5[w] = nfadd(SP5.get(w, Z2), nfmul(c, pad(cf)))
    cl_om = G.quotient5(G.ft_expand_sym_nf(S5, 2))
    cl_N5 = G.quotient5(G.ft_expand_sym_nf(SN5, 2))
    cl_P5 = G.quotient5(G.ft_expand_sym_nf(SP5, 2))
    res = [nfadd(nfadd(a, (-b[0], -b[1])), (-c[0], -c[1])) for a, b, c in zip(cl_om, cl_N5, cl_P5)]
    ok = all(not any(q) for q in res)
    print(f"      class(Omega5) == class(N5) + class(P5): {'PASS' if ok else 'FAIL'}")
    if not ok:
        print("   cl_om :", [G.nf_str(q) for q in cl_om])
        print("   cl_N5 :", [G.nf_str(q) for q in cl_N5])
        print("   cl_P5 :", [G.nf_str(q) for q in cl_P5])
    assert ok
print("ALL EXACT VERIFICATIONS PASS")
