base = '<project>'
here = base + '/research/08_theory/M4b2_w4pp_nice'
disp = open(here + '/nice_display.md').read()
NL = chr(10)
res = """# M4b2: W4'' = R''_nice PROVED — NICE closed form with the seven Theorem-W4 Li4 arguments

## Verdict

**W4'' = R''_nice PROVED exactly.** Same template as M4b (which proved W4'' = R'' with awkward
coefficients over a different Li4 set): stage-1 mod-p CRT with hold-out + stage-2 exact solve +
full 5184-column Fraction verification + standalone fresh verifier (`verify_W4pp_nice.py` prints
`VERDICT: W4'' = R''_nice PROVED`).

The atom set is the SAME 22 weight-4 products as M4b plus Li4 at exactly the seven Theorem-W4
arguments **{b, -b, -b^2, 1/sqrt2, 1/2, 1-1/sqrt2, 2-sqrt2}** (b = sqrt2 - 1); all seven convert via
`cmzv.find_cross`/`cmzv.atom` (step02; none fail). Stage 1 found the 57 columns [28 frame lifts | 29
atoms] **independent** modulo the proven relations, so the coefficients c_k are unique; they came
out **EXACTLY equal to the M4a coefficients (w4pp_closed_form.json), 58/58 components** — the
fix-c_k-and-solve-z_j fallback was therefore not needed (it is implemented in step04 and triggers
automatically if rank < 57 or a mismatch).

W4'' = \\int_b^a [l1^2 l2 \\,\\Omega_{{21}} + l1^3\\,\\Omega_{{30}}]\\,dv
= 1.4244536760088729215144224912572705356370998476755244717254677916845719691528061...
agrees with R''_nice to 1.07e-1060 (1059.97 decimal digits; M4a W4pp_1060.txt).

""" + disp + """

## Proof structure (identical to M4b / A3h5)

Exact word-vector identity on the 5184 admissible weight-4 level-8 words, for X in {rat, sq}:

  W4''_X = sum_k c_k^X B_k + sum_j z_j^X F_j + sum_i y_i^X T0_{S_i}

- T0 = standard weight-4 level-8 relations (regularised double shuffle + Zhao distribution
  d in {2,4,8} + branch identities + shuffle lifts): evaluate to 0.
- F_j = the 28 function-level double-shuffle frame lifts used by A3h5 (provenance re-derived and
  asserted by the verifier): evaluate to 0.
- B_k = 29 classical weight-4 constants (22 products + the 7 NICE Li4 values), exact vectors via
  cmzv.atom/mul.

Hence W4'' = R''_nice = sum_k (c_k^rat + sqrt2 c_k^sq) B_k as real numbers. Certificate:
`W4pp_nice_certificate.json` (7,202,243 bytes; 5097 rows: atom 29, frame 28 [w2 10, w3 9, w4 9],
T0 5040 [std 4810, lift1 223, lift2old 7]).

## Evidence

| step | check | residual | digits |
|---|---|---|---|
| step02 | 7/7 NICE Li4 arguments convert via cmzv.find_cross; 29 atoms admissible weight-4 | — | exact |
| step03 | 29 atom vectors vs direct mpmath (log/polylog/zeta), worst | 3.18e-114 | >=45 required |
| step03 | eval(W4''_rat)+sqrt2 eval(W4''_sq) vs M4a G-value | 3.01e-90 | 90 |
| step04 | stage 1: rank 57 of 57 cols ([28 frames|29 atoms]); consistent at all 8 primes; CRT x7 + hold-out agrees, both components | exact mod p | heights <= 6 digits |
| step04 | c_k (solved) vs M4a w4pp_closed_form.json, term by term | 58/58 components EXACTLY equal | exact |
| step05 | R''_nice (direct mpmath) - W4'' (quad_90.json) | 3.01e-90 | 90 |
| step06 | stage 2 exact fmpz solve (5042x5042; y heights <= 382 digits, 5040 nonzero) + Fraction verification on ALL 5184 columns, both components | **0 (exact)** | exact |
| step08 | independent A3f4 evaluator on full G vectors (supports 4039/4007) | 8.97e-207 / 2.57e-205 | ~200 |
| verify_W4pp_nice.py | fresh rebuild of every row from first principles; residual on all words, both components | **0 (exact)** | exact |
| step09 | R''_nice (mpmath 1150 dps) - W4'' (M4a W4pp_1060.txt) | 1.07e-1060 | **1059.97** |
| step09 | direct eval of M4a w4pp_closed_form.json - certificate R''_nice | 1.74e-1150 | ~1150 |

## Verifier output (fresh run; verifier_nice.log)

```
target: gpl_w4pp.tsv == omega4.json re-derivation (432 terms); weight check True
EXACT rebuilt W4pp component=0 support=3524 residual_nonzero=0
EXACT rebuilt W4pp component=1 support=2638 residual_nonzero=0
EXACT certificate source families={'atom': 29, 'frame_w2': 10, 'frame_w3': 9, 'frame_w4': 9, 'std': 4810, 'lift1': 223, 'lift2old': 7} total_used_rows=5097
VERDICT: W4'' = R''_nice PROVED
```

## Inputs relied on (as in M4b / A3h5)

(i) the G-expression of W4'' (M4a, from the paper's exact IBP data omega4.json; verified to 90 digits
against quadrature there, and re-derived == TSV here); (ii) theorems: regularised double shuffle and
distribution for level-8 multiple polylog values, shuffle products of convergent iterated integrals,
function-level double shuffle in real Moebius frames with absolute-convergence domains;
(iii) correctness of the exact conversion code (cmzv path composition / Cayley pullback; frame_ds) —
audited and validated numerically (every atom to ~1e-114, whole identity to ~1e-205), not formally
verified. The M4a coefficient list is used ONLY as a comparison target: stage 1 solves for c_k
independently (free 57-column solve) and finds them equal.

## Files (work dir research/08_theory/M4b2_w4pp_nice/)

m4b2_common.py, step01_build_vectors.py, step02_build_atoms.py, step03_numeric_checks.py,
step04_stage1_modp.py, step05_numeric_Rpp.py, step06_stage2_exact.py, step07_make_certificate.py,
step08_numeric_sanity.py, step09_sanity_1060.py, step10_emit_display.py,
verify_W4pp_nice.py (standalone exact verifier), W4pp_nice_certificate.json (7.2 MB),
stage1_w4pp_vectors.pkl, stage2_atoms_nice.pkl, stage1_exact_cz_nice.json, stage2_G_nice.pkl,
stage2_y_nice.pkl, nice_display.md, *.out, progress.log, verifier_nice.log, result.md.

## Rerun commands

```
cd research/08_theory/M4b2_w4pp_nice
PY=<repo>/.venv/bin/python
$PY step01_build_vectors.py && $PY step02_build_atoms.py && $PY step03_numeric_checks.py
$PY step04_stage1_modp.py && $PY step05_numeric_Rpp.py
$PY step06_stage2_exact.py && $PY step07_make_certificate.py
$PY verify_W4pp_nice.py     # standalone exact verifier, prints PROVED
$PY step08_numeric_sanity.py && $PY step09_sanity_1060.py && $PY step10_emit_display.py
```
"""
open(here + '/result.md', 'w').write(res)
open(here + '/progress.log', 'a').write(
    'M4b2 COMPLETE: W4\'\' = R\'\'_nice PROVED. rank 57/57; c_k == M4a 58/58; stage2 exact residual 0 all 5184 cols both comps; verify_W4pp_nice.py fresh PROVED; 1060-digit sanity 1.07e-1060.\n')
import os
print('result.md bytes:', os.path.getsize(here + '/result.md'))
print(sorted(os.listdir(here)))
