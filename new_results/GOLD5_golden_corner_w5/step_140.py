verify = r'''"""verify_gold5.py -- < 5 min re-verification of the GOLD5 exact results (PROVED parts).

1. CONTROLS: square pairing == (84, -38 sqrt2) (M4c), hexagon == (41, -33/2) (G2b).
2. Golden pairings: c0=1/phi == (72, 72-54 sqrt5); c0=sqrt5-2 == (0, 972+432 sqrt5).
3. CHECK 1: delta_word == -delta^-(Omega_5) mod shuffles (both corners).
4. All 27 golden Li3(z) classes vanish (Delta has no classical Li3 part).
5. N5: gol1 -(8/3)Li32(-1/phi,-1/(2phi)) + (8/27+2sqrt5/3)Li32(-1,phi^-3);
        gol3 -(9/8+sqrt5/2)Li32(phi^-6,1);
   quotient cobracket == Delta's and obs(N5) == obs(Omega_5), exact.
6. class(Omega_5) == class(N5) + class(P5) in T^5/D^5 (both corners), exact.
7. Convergence of all N5 atoms (|x|<=1, |xy|<1).
Prints PASS/FAIL per item; exits nonzero on any FAIL.
"""
import sys, json, time
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
Fraction = Fr
import gold5_lib as G
from gold5_lib import golden_point, s01, W5Point

t0 = time.time()
ok_all = True
def check(name, cond):
    global ok_all
    ok_all = ok_all and bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name} ({time.time()-t0:.0f}s)", flush=True)
def nfadd(a, b): return (a[0]+b[0], a[1]+b[1])
def nfmul(a, b): return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
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
OBS = {'gol1': ((Fr(72), Fr(0)), (Fr(72), Fr(-54))), 'gol3': ((Fr(0), Fr(0)), (Fr(972), Fr(432)))}

print("verify_gold5:")
# 1. controls
P4 = W5Point("c=1/sqrt2", [-2, 0, 1], (Fr(0), Fr(1, 2)), basis=[(Fr(2), Fr(0)), (Fr(1), Fr(1))])
obs4 = s01.obstruction_fixed(s01.ft_expand_sym(P4.omega5(), 2))
check("square control == (84, -38 sqrt2)", obs4[0] == (Fr(84), Fr(0)) and obs4[1] == (Fr(0), Fr(-38)))
P6 = W5Point("c=1/2", [-1, 1], (Fr(1, 2),))
obs6 = s01.obstruction_fixed(s01.ft_expand_sym(P6.omega5(), 2))
check("hexagon control == (41, -33/2)", obs6[0] == (Fr(41),) and obs6[1] == (Fr(-33, 2),))

for tag in ('gol1', 'gol3'):
    P = golden_point(tag)
    S5 = P.omega5()
    obs = tuple(pad(o) for o in s01.obstruction_fixed(s01.ft_expand_sym(S5, 2)))
    check(f"{tag} pairing == {OBS[tag]}", obs == OBS[tag])
    # 3. CHECK 1
    F = P.F
    Delta = {}
    for w, coef in G.signed_words(P):
        for lw, v in G.delta_word_sym(P, w, F.zero(), P.c0).items():
            Delta[lw] = nfadd(Delta.get(lw, Z2), nfmul(coef, pad(v)))
    D5 = G.delta_minus_5sym(S5)
    keys2 = set(Delta) | set(D5)
    diff = {k: nfadd(pad(Delta.get(k, Z2)), pad(D5.get(k, Z2))) for k in keys2}
    diff = {k: v for k, v in diff.items() if any(v)}
    qc, _ = G.quotient_coords_nf(G.ft_expand_sym_nf(diff, 2))
    check(f"{tag} CHECK 1 (delta_word == -delta^- mod shuffles)", all(not any(q) for q in qc))
    qcD, _ = G.quotient_coords_nf(G.ft_expand_sym_nf({k: v for k, v in Delta.items() if any(v)}, 2))
    # 4. Li3 classes vanish
    li5cl = {eval(k): None for k in json.load(open(GOLD + 'li5_classes_golden.json'))}
    nz = 0
    for z in li5cl:
        qc3, _ = G.quotient_coords_nf(G.ft_expand_sym_nf(G.li3_sym_nf(P, z), 2))
        if any(any(q) for q in qc3): nz += 1
    check(f"{tag} all 27 Li3 classes vanish", nz == 0)
    # 5. N5 cobracket + obs
    Dtot = {}; oQ = Z2; oT = Z2
    for d, x, y in N5[tag]:
        for lw, v in G.delta_word_sym(P, G.li32_word(P, x, y), F.zero(), F.one()).items():
            Dtot[lw] = nfadd(Dtot.get(lw, Z2), nfmul(d, pad(v)))
        S_at = G.ft_expand_sym_nf(G.li32_sym_nf(P, x, y), 2)
        ob = G.obs_of_SE(S_at)
        oQ = nfadd(oQ, nfmul(d, pad(ob[0]) if ob[0] else Z2))
        oT = nfadd(oT, nfmul(d, pad(ob[1]) if ob[1] else Z2))
    qcN, _ = G.quotient_coords_nf(G.ft_expand_sym_nf(Dtot, 2))
    check(f"{tag} N5 cobracket == Delta", tuple(pad(q) for q in qcN) == tuple(pad(q) for q in qcD))
    check(f"{tag} obs(N5) == obs(Omega_5)", (oQ, oT) == obs)
    # 6. class match
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
    check(f"{tag} class(Omega_5) == class(N5) + class(P5)", all(not any(q) for q in res))
    # 7. convergence
    import mpmath as mp
    mp.mp.dps = 30
    r5 = mp.sqrt(5)
    def fv(z):
        return abs(mp.mpf(z[0].numerator)/z[0].denominator + mp.mpf(z[1].numerator)/z[1].denominator*r5)
    conv = all(fv(x) <= 1 and fv(F.mul(x, y)) < 1 for d, x, y in N5[tag])
    check(f"{tag} N5 atoms convergent (paper convention)", conv)
print(f"verify_gold5: {'ALL PASS' if ok_all else 'FAILURES'} ({time.time()-t0:.0f}s)")
sys.exit(0 if ok_all else 1)
'''
open('verify_gold5.py', 'w').write(verify)
import subprocess
r = subprocess.run(['nice', '-n', '10', '<project>/.venv/bin/python', 'verify_gold5.py'], capture_output=True, text=True, cwd='.', timeout=300)
print(r.stdout)
print("STDERR:", r.stderr[-800:])
