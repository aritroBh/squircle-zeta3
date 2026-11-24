# M4b2: W4'' = R''_nice PROVED — NICE closed form with the seven Theorem-W4 Li4 arguments

## Verdict

**W4'' = R''_nice PROVED exactly**, reusing M4b's template (stage-1 mod-p CRT with hold-out + stage-2 exact solve + full 5184-column Fraction verification + standalone fresh verifier printing `VERDICT: W4'' = R''_nice PROVED`). The atom set is the same 22 weight-4 products plus Li4 at exactly **{b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2}** (b = √2−1) — all 7 convert via `cmzv.find_cross`/`cmzv.atom` (none fail). Stage 1 found the 57 columns [28 frame lifts | 29 atoms] **independent** modulo the proven relations, so c_k are unique, and they came out **exactly equal to the M4a coefficients, 58/58 components** — the fix-c/solve-z fallback (implemented in step04) was not needed.

## R''_nice (Theorem-W4-style display; λ = log 2, L = log(1+√2), b = √2−1)

$$
W_4'' = \bigl(1-104\sqrt2\bigr)\,\Li_4(b) - \bigl(1+104\sqrt2\bigr)\,\Li_4(-b) + \tfrac{41}{4}\sqrt2\,\Li_4(-b^2) + 2\sqrt2\,\Li_4\bigl(\tfrac1{\sqrt2}\bigr) - \bigl(\tfrac{37}{4}+\tfrac{1}{8}\sqrt2\bigr)\Li_4\bigl(\tfrac12\bigr) + \bigl(10+16\sqrt2\bigr)\Li_4\bigl(1-\tfrac1{\sqrt2}\bigr) - \bigl(10-16\sqrt2\bigr)\Li_4(2-\sqrt2)
$$
$$
-\bigl(\tfrac{37}{96}-\tfrac{13}{48}\sqrt2\bigr)\lambda^4 + \bigl(\tfrac{7}{12}+\tfrac{13}{12}\sqrt2\bigr)\lambda^3 L + \bigl(\tfrac{3}{2}+\tfrac{1}{8}\sqrt2\bigr)\lambda^2 L^2 - \bigl(\tfrac{7}{6}+\tfrac{7}{3}\sqrt2\bigr)\lambda L^3 - \bigl(\tfrac{19}{24}+\tfrac{1}{12}\sqrt2\bigr)L^4 + \bigl(\tfrac{19}{48}-\tfrac{157}{192}\sqrt2\bigr)\pi^2\lambda^2 - \bigl(\tfrac{23}{16}-\tfrac{7}{6}\sqrt2\bigr)\pi^2\lambda L + \bigl(\tfrac{49}{48}+\tfrac{11}{48}\sqrt2\bigr)\pi^2 L^2 - \tfrac{7}{4}\sqrt2\Li_2(b)\lambda^2 - 4\Li_2(b)\lambda L + \bigl(5+\sqrt2\bigr)\Li_2(b)L^2 + \bigl(\tfrac{1333}{11520}-\tfrac{391}{2880}\sqrt2\bigr)\pi^4 - \bigl(\tfrac{19}{24}+\tfrac{1}{6}\sqrt2\bigr)\pi^2\Li_2(b) + 5\Li_2(b)^2 + \bigl(\tfrac{1}{16}+\tfrac{77}{64}\sqrt2\bigr)\zeta(3)\lambda + \bigl(\tfrac{63}{32}+\tfrac{3}{8}\sqrt2\bigr)\zeta(3)L + 6\sqrt2\Li_3(b)\lambda + \bigl(10+2\sqrt2\bigr)\Li_3(b)L - 6\sqrt2\Li_3(-b)\lambda - \bigl(10-2\sqrt2\bigr)\Li_3(-b)L - 11\sqrt2\Li_3\bigl(\tfrac1{\sqrt2}\bigr)\lambda - 18\Li_3\bigl(\tfrac1{\sqrt2}\bigr)L
$$

## Coefficient comparison with M4a (solved c_k vs w4pp_closed_form.json): 58/58 components exactly equal

| atom | c_rat | c_sq | | atom | c_rat | c_sq |
|---|---|---|---|---|---|---|
| λ⁴ | −37/96 | 13/48 | | ζ(3)λ | 1/16 | 77/64 |
| λ³L | 7/12 | 13/12 | | ζ(3)L | 63/32 | 3/8 |
| λ²L² | 3/2 | 1/8 | | Li3(b)λ | 0 | 6 |
| λL³ | −7/6 | −7/3 | | Li3(b)L | 10 | 2 |
| L⁴ | −19/24 | −1/12 | | Li3(−b)λ | 0 | −6 |
| π²λ² | 19/48 | −157/192 | | Li3(−b)L | −10 | 2 |
| π²λL | −23/16 | 7/6 | | Li3(1/√2)λ | 0 | −11 |
| π²L² | 49/48 | 11/48 | | Li3(1/√2)L | −18 | 0 |
| Li2(b)λ² | 0 | −7/4 | | Li4(b) | 1 | −104 |
| Li2(b)λL | −4 | 0 | | Li4(−b) | −1 | −104 |
| Li2(b)L² | 5 | 1 | | Li4(−b²) | 0 | 41/4 |
| π⁴ | 1333/11520 | −391/2880 | | Li4(1/√2) | 0 | 2 |
| π²Li2(b) | −19/24 | −1/6 | | Li4(1/2) | −37/4 | −1/8 |
| Li2(b)² | 5 | 0 | | Li4(1−1/√2) | 10 | 16 |
| | | | | Li4(2−√2) | −10 | 16 |

(49 nonzero terms, max |coefficient| 104 — as M4a. Per-term record stored in `stage1_exact_cz_nice.json` → `comparison`.)

## Evidence table

| step | check | residual | digits |
|---|---|---|---|
| step02 | 7/7 NICE Li4 args convert via cmzv.find_cross; 29 atoms admissible weight-4 | — | exact |
| step03 | 29 atom vectors vs direct mpmath, worst | 3.18e-114 | ≥45 required |
| step03 | eval(W4'') vs M4a G-value | 3.01e-90 | 90 |
| step04 | stage 1: rank **57/57**; consistent at all 8 primes; CRT×7 + hold-out, both components | exact mod p | heights ≤ 6 digits |
| step04 | c_k vs M4a, term by term | **58/58 exactly equal** | exact |
| step05 | R''_nice − W4'' (quad_90) | 3.01e-90 | 90 |
| step06 | stage-2 exact fmpz solve (5042×5042; y heights ≤ 382 digits, 5040 nonzero) + Fraction check on ALL 5184 columns, both components | **0 (exact)** | exact |
| verify_W4pp_nice.py | fresh rebuild from first principles; residual on all words, both components | **0 (exact)** | exact |
| step08 | independent A3f4 evaluator on G_rat/G_sq | 8.97e-207 / 2.57e-205 | ~200 |
| step09 | R''_nice (1150 dps) − W4'' (M4a W4pp_1060.txt) | 1.07e-1060 | **1059.97** |
| step09 | direct eval of M4a closed form − certificate R''_nice | 1.74e-1150 | ~1150 |

Verifier output (fresh run): `target: gpl_w4pp.tsv == omega4.json re-derivation (432 terms)`; `EXACT rebuilt W4pp component=0/1 residual_nonzero=0`; families `{atom: 29, frame_w2: 10, frame_w3: 9, frame_w4: 9, std: 4810, lift1: 223, lift2old: 7}`, total 5097 rows; `VERDICT: W4'' = R''_nice PROVED`. Certificate `W4pp_nice_certificate.json` = 7,202,243 bytes.

## Files (new_results/M4b2_w4pp_nice/)

`m4b2_common.py`, `step01_build_vectors.py`, `step02_build_atoms.py`, `step03_numeric_checks.py`, `step04_stage1_modp.py`, `step05_numeric_Rpp.py`, `step06_stage2_exact.py`, `step07_make_certificate.py`, `step08_numeric_sanity.py`, `step09_sanity_1060.py`, `step10_emit_display.py`, `verify_W4pp_nice.py`, `W4pp_nice_certificate.json`, pickles/json (`stage1_w4pp_vectors.pkl`, `stage2_atoms_nice.pkl`, `stage1_exact_cz_nice.json`, `stage2_G_nice.pkl`, `stage2_y_nice.pkl`), `nice_display.md`, `result.md`, `progress.log`, `verifier_nice.log`, `*.out`.

## Rerun

```
cd new_results/M4b2_w4pp_nice; PY=<repo>/.venv/bin/python
$PY step01_build_vectors.py && $PY step02_build_atoms.py && $PY step03_numeric_checks.py
$PY step04_stage1_modp.py && $PY step05_numeric_Rpp.py && $PY step06_stage2_exact.py && $PY step07_make_certificate.py
$PY verify_W4pp_nice.py && $PY step08_numeric_sanity.py && $PY step09_sanity_1060.py && $PY step10_emit_display.py
```

Inputs relied on: identical to M4b (M4a G-expression of W4'' == omega4.json IBP data; regularised double shuffle + distribution at level 8; frame double shuffle; audited cmzv conversion, numerically validated to ~1e-114 per atom and ~1e-205 for the whole identity).
## Planner verification (W4 side session, 2026-10-22 18:4x)
- `diff ../M4b_w4pp_exact/verify_W4pp.py verify_W4pp_nice.py`: only the docstring, paths, log file and verdict text
  differ. The atoms are rebuilt from the certificate's encodings with cmzv, so it still rebuilds everything from
  first principles.
- Fresh planner re-run (`planner_rerun_verify.out`): "VERDICT: W4'' = R''_nice PROVED", exit 0.
- The coefficients equal M4a's LLL closed form (`../M4a_decomposition/w4pp_closed_form.json`). The planner had
  independently evaluated that JSON with mpmath and got |R'' − W4''| = 1.07e-1060. nice_display.md was spot-checked
  against the JSON.
- STATUS: W4'' (the weight-4 part of M4) = R''_nice is PROVED. It uses the same 7 Li4 arguments as Theorem W4, and
  carries the same proof caveat: correctness of the cmzv and frame_ds conversion code (audited, numerically validated).
