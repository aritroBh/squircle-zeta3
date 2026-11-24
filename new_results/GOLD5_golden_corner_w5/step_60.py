import sys, time
sys.path.insert(0, '<project>/research/08_theory/GOLD5_golden_corner_w5')
from fractions import Fraction as Fr
import gold5_lib as G

t0 = time.time()
for tag in ('gol1', 'gol3'):
    P = G.golden_point(tag)
    F = P.F
    # assemble Delta from signed words
    Delta = {}
    for w, coef in G.signed_words(P):
        dw = G.delta_word_sym(P, w, F.zero(), P.c0)
        for lw, v in dw.items():
            Delta[lw] = F.add(Delta.get(lw, F.zero()), F.mul(coef, v))
    Delta = {k: v for k, v in Delta.items() if any(v)}
    print(f"{tag}: Delta assembled: {len(Delta)} left words ({time.time()-t0:.0f}s)", flush=True)
    # EXACT CHECK 1: Delta + delta^-(Omega5) is a shuffle product
    S5 = P.omega5()
    D5 = G.delta_minus_5sym(S5)
    keys = set(Delta) | set(D5)
    diff = {k: F.add(Delta.get(k, F.zero()), D5.get(k, F.zero())) for k in keys}
    diff = {k: v for k, v in diff.items() if any(v)}
    # quotient test: coords must vanish
    Dmat, piv = G.d3_rref()
    B = G.t3_basis()
    from cor6_sym import sshuffle
    Drows = G.d3_rows()
    # expand diff words (exponent vectors) to basis indices
    diffE = G.ft_expand_sym_nf(diff, 2)
    qc, free = G.quotient_coords_nf(diffE)
    chk1 = all(not any(q) for q in qc)
    print(f"  EXACT CHECK 1 (delta_word == -delta^-(5-sym) mod shuffles): {'PASS' if chk1 else 'FAIL'}", flush=True)
    # quotient coords of Delta
    DeltaE = G.ft_expand_sym_nf(Delta, 2)
    qcD, freeD = G.quotient_coords_nf(DeltaE)
    print(f"  Delta quotient coords (free words {[B[j] for j in freeD]}):", flush=True)
    for j, q in zip(freeD, qcD):
        print(f"    word {B[j]}: {G.nf_str(q)}   {[str(x) for x in q]}", flush=True)
