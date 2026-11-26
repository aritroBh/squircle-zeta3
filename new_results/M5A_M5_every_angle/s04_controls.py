"""s04: controls for s03 (exact).
(1) selftest of the HPL side: h_w(2n+delta) as Z-forms (g8_zsum.hpl_coeff_at) against the exact Taylor coefficients
    from the definition (g8_allN.hpl_coeffs) for every word in the k=5 term lists and its negation, delta=-1..5, n=1..6.
(2) mutants: change one coefficient (or one prefactor exponent) of a term list; the comparison of s03 must fail."""
import sys, json, random
sys.setrecursionlimit(20000)
from fractions import Fraction as Fr
import g8_allN as G, g8_zsum as zs
from s03_prove import parse, check_part
words = set()
for nm in ['1', 'z2', 'z2^2', 'z3', 'z4', 'z5', 'z2*z3']:
    for kind, w, e, q in parse(f"hpl_fit_k5_{nm}.json"):
        if kind != 'const': words.add(w); words.add(G.neg(w))
ok = True; cnt = 0
for w in sorted(words):
    co = G.hpl_coeffs(w, 2*6 + 6)
    for nv in range(1, 7):
        for delta in range(-1, 6):
            cnt += 1
            if G.ev(zs.hpl_coeff_at(w, delta), nv) != co[2*nv + delta]:
                ok = False; print('FAIL', w, nv, delta)
print(f"selftest: {len(words)} words, {cnt} checks of h_w(2n+delta) against the definition: all equal = {ok}", flush=True)
P5 = G.mu_side(5)
Z = {'1': (0,0,0,0), 'z2': (1,0,0,0), 'z2^2': (2,0,0,0)}
random.seed(1)
for nm in ['1', 'z2', 'z2^2']:
    T = parse(f"hpl_fit_k5_{nm}.json")
    i = random.randrange(len(T))
    kind, w, e, q = T[i]
    T1 = list(T); T1[i] = (kind, w, e, q + 1)
    print(f"mutant {nm}: coefficient of term {i} ({kind} e={e} w={w}) +1 ->")
    r1 = check_part(nm + '[mut coef]', T1, P5[Z[nm]])
    j = next(j for j, t in enumerate(T) if t[0] != 'const' and t[2] > -4)
    kind, w, e, q = T[j]
    T2 = list(T); T2[j] = (kind, w, e - 1, q)
    print(f"mutant {nm}: exponent of term {j} lowered by one ->")
    r2 = check_part(nm + '[mut exp]', T2, P5[Z[nm]])
    print(f"   mutants detected: {not r1}, {not r2}", flush=True)
