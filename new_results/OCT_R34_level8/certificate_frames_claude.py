"""Exact certificate: DIFF(Lam), DIFF(swap) in span_Q(rational base rows + A3h4 weight-3 frame rows), if they are rational.
Run from research/08_theory."""
import sys, pickle
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
exec(open('OCT_R34_level8/membership_ext_claude.py').read().split("F5 = pickle.load")[0])
from qsqrt2 import QS2
FR = pickle.load(open('A3h4_frame_ds/validated_w3.pickle', 'rb'))['rows']
irr = sum(1 for r in FR for w, a in r if isinstance(a, QS2) and a.q != 0)
print(f"frame rows: {len(FR)}, coefficients with a sqrt2 part: {irr}")
Frows = [accum([(tuple(w), (a.p if isinstance(a, QS2) else a)) for w, a in r]) for r in FR]
ROWS = Wrows + Frows
tag = ['base'] * len(Wrows) + ['frame'] * len(Frows)
for key in ('Lam', 'swap'):
    DIFF = DIFFS[key]
    p = 1000000033
    M = [[0] * n for _ in ROWS]
    for i, r in enumerate(ROWS):
        for w, a in r.items(): M[i][idx[w]] = (M[i][idx[w]] + modp(a, p)) % p
    d = [modp(DIFF.get(c, 0), p) for c in cols]
    r0 = flint.nmod_mat(len(M), n, [v for r in M for v in r], p).rank()
    r1 = flint.nmod_mat(len(M) + 1, n, [v for r in M for v in r] + d, p).rank()
    print(f"{key}: rational base + frames rank {r0}, +DIFF {r1}")
    if r0 != r1: continue
    E, rk = flint.nmod_mat(n, len(M), [M[i][j] for j in range(n) for i in range(len(M))], p).rref()
    piv = [next(j for j in range(len(M)) if int(E[i, j]) != 0) for i in range(rk)]
    E2, _ = flint.nmod_mat(rk, n, [v for i in piv for v in M[i]], p).rref()
    pc = [next(j for j in range(n) if int(E2[i, j]) != 0) for i in range(rk)]
    fq = lambda f: flint.fmpq(Fr(f).numerator, Fr(f).denominator)
    sol = flint.fmpq_mat(rk, rk, [fq(ROWS[piv[i]].get(cols[j], 0)) for j in pc for i in range(rk)]).solve(
        flint.fmpq_mat(rk, 1, [fq(DIFF.get(cols[j], 0)) for j in pc]))
    c = [Fr(int(sol[i, 0].p), int(sol[i, 0].q)) for i in range(rk)]
    comb = {}
    for ci, i in zip(c, piv):
        if ci:
            for w, a in ROWS[i].items(): comb[w] = comb.get(w, 0) + ci * a
    comb = {w: a for w, a in comb.items() if a != 0}
    used = [i for ci, i in zip(c, piv) if ci]
    print(f"{key}: EXACT CERTIFICATE {'holds on all 576 words' if comb == DIFF else 'FAILED'}; rows used: "
          f"{sum(tag[i] == 'base' for i in used)} base (max index {max([i for i in used if tag[i]=='base'], default=None)}), "
          f"{sum(tag[i] == 'frame' for i in used)} frame")
    pickle.dump({'rows_used': used, 'multipliers': c, 'piv': piv, 'ok': comb == DIFF}, open(f'OCT_R34_level8/certificate_frames_{key}.pickle', 'wb'))
