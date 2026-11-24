step03 = r'''"""step_03.py -- GOLD5 Task 3 (exact symbol algebra): the (3,2) cobracket Delta of R4^(5) at both
golden corners, the convergent Li_{3,2} correction N5, and the classical Li5 class of W5 - N5.

EXACT (PROVED at the symbol level):
  1. CHECK 1: delta_word(sum) == -delta^-(Omega_5) mod shuffle products (quotient T^3/D^3).
  2. Delta's quotient coords: gol1 (-72, 72-54 sqrt5), gol3 (0, 972+432 sqrt5).
     The Li3-part: EVERY Li3(z) at the 27 golden {2}-unit solutions of x+y=1 has ZERO class in
     T^3/D^3 -- so Delta has NO classical Li3 part (unlike square/hexagon).
  3. N5 with convergent atoms (|x|<=1, |xy|<1, paper convention sum_{n>m} x^n y^m/(n^3 m^2)):
       gol1: N5 = -(1/4) Li32(-rho^3/4, 1) + (4/9 + 2 sqrt5/3) Li32(-1, rho^3)
       gol3: N5 = -(9/8 + sqrt5/2) Li32(rho^6, 1)          [rho = 1/phi]
     verified: quotient cobracket == Delta's, obs(N5) == obs(Omega_5) (exact).
  4. class(Omega_5 - N5) in T^5/D^5 lies in the Li5 span (exact solve); Li5 reps P5 recorded.
"""
import sys, itertools, json, time
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
from gold5_lib import golden_point, log, nf_str

t0 = time.time()
Fraction = Fr
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}

def nfadd(a, b): return (a[0]+b[0], a[1]+b[1])
def nfmul(a, b): return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def pad(a): return a if len(a) == 2 else (a[0], Fr(0))
Z2 = (Fr(0), Fr(0))
def fval(z):
    import mpmath as mp
    return float(mp.mpf(z[0].numerator)/z[0].denominator + mp.mpf(z[1].numerator)/z[1].denominator*mp.sqrt(5))

RHO3 = (Fr(-2), Fr(1)); RHO6 = (Fr(9), Fr(-4)); M1 = (Fr(-1), Fr(0)); ONE = (Fr(1), Fr(0))
MR3O4 = (Fr(1, 2), Fr(-1, 4))   # -rho^3/4 = (2-sqrt5)/4
N5 = {
 'gol1': [((Fr(-1, 4), Fr(0)), MR3O4, ONE), ((Fr(4, 9), Fr(2, 3)), M1, RHO3)],
 'gol3': [((Fr(-9, 8), Fr(-1, 2)), RHO6, ONE)],
}
P5 = {
 'gol1': [((Fr(-233, 5), Fr(54, 5)), (Fr(3, 4), Fr(-1, 4))),
          ((Fr(-23, 5), Fr(-6, 5)), (Fr(1, 4), Fr(1, 4))),
          ((Fr(-107, 45), Fr(233, 135)), RHO3),
          ((Fr(-32, 5), Fr(-4, 5)), (Fr(3, 1), Fr(-1, 1)))],
 'gol3': [((Fr(203, 4), Fr(-39, 4)), (Fr(3, 4), Fr(-1, 4))),
          ((Fr(103, 4), Fr(-27, 4)), (Fr(1, 4), Fr(-1, 4))),
          ((Fr(-149, 4), Fr(33, 4)), (Fr(1, 4), Fr(1, 4))),
          ((Fr(1241, 324), Fr(-103, 108)), RHO3)],
}

out = {}
for tag in ('gol1', 'gol3'):
    P = golden_point(tag)
    F = P.F
    # --- cobracket of Omega_5 from signed words
    Delta = {}
    for w, coef in G.signed_words(P):
        dw = G.delta_word_sym(P, w, F.zero(), P.c0)
        for lw, v in dw.items():
            Delta[lw] = nfadd(Delta.get(lw, Z2), nfmul(coef, pad(v)))
    Delta = {k: v for k, v in Delta.items() if any(v)}
    # CHECK 1
    S5 = P.omega5()
    D5 = G.delta_minus_5sym(S5)
    keys = set(Delta) | set(D5)
    diff = {k: nfadd(pad(Delta.get(k, Z2)), pad(D5.get(k, Z2))) for k in keys}
    diff = {k: v for k, v in diff.items() if any(v)}
    qc, free = G.quotient_coords_nf(G.ft_expand_sym_nf(diff, 2))
    chk1 = all(not any(q) for q in qc)
    log(f"{tag}: CHECK 1 (delta_word == -delta^-(5-sym) mod shuffles): {'PASS' if chk1 else 'FAIL'}")
    assert chk1
    qcD, _ = G.quotient_coords_nf(G.ft_expand_sym_nf(Delta, 2))
    log(f"  Delta quotient coords: {[nf_str(q) for q in qcD]}")
    # --- Li3-part: all Li3(z) classes
    li3_zero = 0; li3_tot = 0
    for zs in li5cl:
        z = eval(zs)
        qc3, _ = G.quotient_coords_nf(G.ft_expand_sym_nf(G.li3_sym_nf(P, z), 2))
        li3_tot += 1
        if all(not any(q) for q in qc3): li3_zero += 1
    log(f"  Li3(z) classes: {li3_zero}/{li3_tot} vanish (Delta has no classical Li3 part)")
    # --- N5 verification
    Dtot = {}
    oQ = Z2; oT = Z2
    for d, x, y in N5[tag]:
        dw = G.delta_word_sym(P, G.li32_word(P, x, y), F.zero(), F.one())
        for lw, v in dw.items():
            Dtot[lw] = nfadd(Dtot.get(lw, Z2), nfmul(d, pad(v)))
        S_at = G.ft_expand_sym_nf(G.li32_sym_nf(P, x, y), 2)
        obs = G.obs_of_SE(S_at)
        oQ = nfadd(oQ, nfmul(d, pad(obs[0]) if obs[0] else Z2))
        oT = nfadd(oT, nfmul(d, pad(obs[1]) if obs[1] else Z2))
    qcN, _ = G.quotient_coords_nf(G.ft_expand_sym_nf(Dtot, 2))
    ok_co = tuple(pad(q) for q in qcN) == tuple(pad(q) for q in qcD)
    log(f"  N5 quotient cobracket == Delta: {'PASS' if ok_co else 'FAIL'}  ({[nf_str(q) for q in qcN]})")
    # obs of Omega_5
    SE = G.ft_expand_sym_nf(S5, 2)
    obs_om = G.obs_of_SE(SE)
    obs_om = tuple(pad(o) for o in obs_om)
    ok_obs = (oQ, oT) == obs_om
    log(f"  obs(N5) == obs(Omega_5) = [{nf_str(obs_om[0])}, {nf_str(obs_om[1])}]: {'PASS' if ok_obs else 'FAIL'}")
    assert ok_co and ok_obs
    # --- class of Omega_5 - N5 in Li5 span
    SN5 = {}
    for d, x, y in N5[tag]:
        for w, cf in G.li32_sym_nf(P, x, y).items():
            SN5[w] = nfadd(SN5.get(w, Z2), nfmul(d, pad(cf if isinstance(cf, tuple) else (cf,))))
    cl_om = G.quotient5(G.ft_expand_sym_nf(S5, 2))
    cl_N5 = G.quotient5(G.ft_expand_sym_nf(SN5, 2))
    cl_T = tuple(nfadd(a, (-b[0], -b[1])) for a, b in zip(cl_om, cl_N5))
    # verify P5's class == cl_T
    SP5 = {}
    for c, z in P5[tag]:
        for w, cf in G.li5_sym_nf(P, z).items():
            SP5[w] = nfadd(SP5.get(w, Z2), nfmul(c, pad(cf)))
    cl_P5 = G.quotient5(G.ft_expand_sym_nf(SP5, 2))
    ok_cl = all(nfadd(a, (-b[0], -b[1])) == Z2 for a, b in zip(cl_T, cl_P5))
    log(f"  class(Omega_5 - N5) == class(P5) in T^5/D^5: {'PASS' if ok_cl else 'FAIL'}")
    log(f"    class = {[nf_str(q) for q in cl_T]}")
    assert ok_cl
    out[tag] = {'delta_qc': [[str(c) for c in q] for q in qcD],
                'li3_classes_zero': f'{li3_zero}/{li3_tot}',
                'obs': [[str(c) for c in o] for o in obs_om],
                'N5': [[list(map(str, d)), str(x), str(y)] for d, x, y in N5[tag]],
                'P5': [[list(map(str, c)), str(z)] for c, z in P5[tag]],
                'class_T': [[str(c) for c in q] for q in cl_T]}
json.dump(out, open(GOLD + 'cobracket_golden.json', 'w'), indent=1)
log(f"step_03 done ({time.time()-t0:.0f}s); saved cobracket_golden.json")
'''
open('step_03.py', 'w').write(step03)
import subprocess
r = subprocess.run(['nice', '-n', '10', '<project>/.venv/bin/python', 'step_03.py'], capture_output=True, text=True, cwd='.', timeout=1200)
print(r.stdout)
print("STDERR:", r.stderr[-2000:])
