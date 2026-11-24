import os
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
result = r'''# M4b: W4'' = R'' PROVED — exact closed form and exact linear-algebra certificate (no LLL)

## Verdict

**W4'' = R'' PROVED exactly**, in one linear-algebra computation (A3h5 template: stage-1 mod-p CRT with
hold-out + stage-2 exact solve + full 5184-column Fraction verification + standalone verifier).

W4'' = ∫_b^a [l1²l2 Ω21 + l1³ Ω30] dv = 1.4244536760088729215144224912572705356370998476755244717254677916845719691528061…
equals the following Q(√2)-combination of classical weight-4 constants
(λ = log 2, L = log(1+√2), b = √2−1; coefficient = c_rat + c_sq·√2):

| atom B_k | c_rat | c_sq |
|---|---|---|
| Li4((1+√2/2)/2) = Li4(1/2+√2/4) | −50/29 | −39/29 |
| Li4(−2+2√2) | −190/29 | −355/116 |
| Li4(−√2/2) | 474/29 | 350/29 |
| Li4(√2/2) | 6 | 6 |
| Li4(2−√2) | 40/29 | 8/29 |
| Li4(1−√2) | −289/29 | −17/58 |
| Li4(b) | −131/29 | 139/58 |
| λ⁴ | −923/2784 | −83/11136 |
| λ³L | 39/29 | 699/464 |
| λ²L² | −11/116 | −757/464 |
| λL³ | −21/58 | −59/174 |
| L⁴ | −237/232 | 55/696 |
| π²λ² | 971/1392 | −389/5568 |
| π²λL | −2561/1392 | 59/348 |
| π²L² | 1741/1392 | 31/464 |
| Li2(b)λ² | 0 | −7/4 |
| Li2(b)λL | −4 | 0 |
| Li2(b)L² | 5 | 1 |
| π⁴ | 50017/334080 | 4777/55680 |
| π²Li2(b) | −19/24 | −1/6 |
| Li2(b)² | 5 | 0 |
| ζ(3)λ | 1/16 | 77/64 |
| ζ(3)L | 63/32 | 3/8 |
| Li3(b)λ | 0 | 6 |
| Li3(b)L | 10 | 2 |
| Li3(−b)λ | 0 | −6 |
| Li3(−b)L | −10 | 2 |
| Li3(1/√2)λ | 0 | −11 |
| Li3(1/√2)L | −18 | 0 |

Numeric: R'' − W4'' = 3.01e-90 (all 90 digits of quad_90.json), direct mpmath at 110 dps.

## Proof structure (identical template to A3h5's proof of Conjecture 4)

Exact word-vector identity on the 5184 admissible weight-4 level-8 words, for X ∈ {rat, sq}:

  W4''_X = Σ_k c_k^X B_k + Σ_j z_j^X F_j + Σ_i y_i^X T0_{S_i}

- T0 = standard weight-4 level-8 relations (regularised double shuffle + Zhao distribution d∈{2,4,8} +
  branch identities + shuffle lifts): evaluate to 0.
- F_j = the 28 function-level double-shuffle frame lifts used by A3h5 (absolute-convergence domains,
  provenance re-derived and asserted): evaluate to 0.
- B_k = 29 classical weight-4 constants (22 products + 7 Li4 values), exact vectors via cmzv.atom/mul.

Hence W4'' = R'' = Σ_k (c_k^rat + √2 c_k^sq) B_k as real numbers.

## Atoms

- 22 weight-4 products of {log2, L}, {π², Li2(b)}, {ζ(3), Li3(b), Li3(−b), Li3(1/√2)} (m3_pslq.py list).
- Li4 candidates: 17 real S-units z<1 with 1−z an S-unit (j∈[−6,6], m∈[−4,4], sign ±).
  find_cross converts 14; does NOT convert 3: (−16+12√2), (1/2−3√2/8), (17−12√2)
  (not cross-ratios of the 8 real points with the convergence condition). The 7 W4 closed-form
  arguments (b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2) all convert. Stage 1 used 7 of the 14
  (columns 23,24,25,26,27,29,30 of the atom list; lexicographic independent set, frames first).

## Evidence

| step | check | residual | digits |
|---|---|---|---|
| step03 | 36 atom vectors vs mpmath (log/polylog/ζ), worst | 3.18e-114 | ≥45 required |
| step03 | eval(W4''_rat)+√2·eval(W4''_sq) vs M4a G-value (= quad at 90 dps) | 3.01e-90 | 90 |
| step04 | stage 1: rank 57 of 64 cols ([28 frames|36 atoms]); consistent at all 8 primes; CRT×7 + hold-out agrees, both components | exact mod p | heights ≤ 7 digits |
| step05 | R'' (direct mpmath) − W4'' (quad_90) | 3.01e-90 | 90 |
| step06 | stage 2 exact fmpz solve (5042×5042, y heights ≤ 383 digits, 5040 nonzero) + Fraction verification on ALL 5184 columns, both components | **0 (exact)** | exact |
| step08 | independent evaluator on full G vectors (supports 4042/4016) | 3.5e-206 / 2.2e-206 | ~200 |
| verify_W4pp.py | fresh rebuild of every row from first principles; residual on all words, both components | **0 (exact)** | exact |

Certificate source families: atom 29, frame 28 (w2 10, w3 9, w4 9), T0 5040 (std 4810, lift1 223, lift2old 7); 5097 rows; certificate 7.2 MB.

## Inputs relied on (as in A3h5)

(i) the G-expression of W4'' (M4a, from the paper's exact IBP data omega4.json; verified to 90 digits
against quadrature there, and re-derived == TSV here); (ii) theorems: regularised double shuffle and
distribution for level-8 multiple polylog values, shuffle products of convergent iterated integrals,
function-level double shuffle in real Möbius frames with absolute-convergence domains;
(iii) correctness of the exact conversion code (cmzv path composition / Cayley pullback; frame_ds) —
audited and validated numerically (every atom to ~1e-114, whole identity to ~1e-206), not formally verified.

## Files (work dir research/08_theory/M4b_w4pp_exact/)

m4b_common.py, step01_build_vectors.py, step02_build_atoms.py, step03_numeric_checks.py,
step04_stage1_modp.py, step05_numeric_Rpp.py, step06_stage2_exact.py, step07_make_certificate.py,
step08_numeric_sanity.py, verify_W4pp.py (standalone verifier), W4pp_certificate.json (7.2 MB),
stage1_w4pp_vectors.pkl, stage2_atoms.pkl, stage1_exact_cz.json, stage2_G.pkl, stage2_y.pkl,
progress.log, verifier.log, result.md.

## Rerun commands

```
cd research/08_theory/M4b_w4pp_exact
PY=<repo>/.venv/bin/python
$PY step01_build_vectors.py && $PY step02_build_atoms.py && $PY step03_numeric_checks.py
$PY step04_stage1_modp.py && $PY step05_numeric_Rpp.py
$PY step06_stage2_exact.py && $PY step07_make_certificate.py
$PY verify_W4pp.py        # standalone exact verifier, ~17 s
$PY step08_numeric_sanity.py
```
'''
open('result.md','w').write(result)
with open('progress.log','a') as f:
    f.write('2026-10-22 18:24:30 RESULT: W4'' = R'' PROVED. stage1 CRT+hold-out OK; stage2 exact 5184-col verification residual 0 both components; verify_W4pp.py fresh run PASS (16.7 s); R''-W4'' = 3.01e-90 (90 digits); atoms 22 products + 7 Li4; see result.md\n')
print('result.md written', len(result))
