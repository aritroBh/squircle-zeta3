# ZQ5d — a cleaner weight-5 formula for L(χ₁₃,5) from the capped pool (all NUMERICAL/LLL-based, not proofs)

Workdir: `new_results/ZQ5d_zagier_d13_clean`. Libs copied (never edited): `w5blib.py`, `w5pipe.py`, `w5verify.py` (from ZQ5c), `planner_li5.py` (from ZQ5 planner). Data copied to `out/d13/`. Precision traps respected (`mp.mp.dps`, `gram='exact'`, Python ints). Convention throughout: **`T_coef·L(χ13,5)/√13 + z5_coef·ζ(5) + Σ li5 + Σ products = 0`** (same side as `d17_final.json`), with letters `lE13=(3+√13)/2, ls13=√13, l2=2, lp3=(5+√13)/2, lp3c=(5−√13)/2`.

## 1. Step 1 (gcd relation) — COMPLETE and verified

- The 9 exact relations of ZQ5c's `genuine_lattice_cap.json` (74 atoms `[13⁵·T, ζ(5), 32 combos K_j−P_j, 35 π²·log³, 5 π⁴·log]`) have T-coefficients (rows 4–8) `4155177511991, 4302263441619, 5184779019387, 19452114193303, 35594794969976`; **gcd = 36771482407 = 11·61·941·58237** (verified by exact integer gcd + `sympy.factorint`).
- Extended gcd: **`29·row4 − 28·row5`** gives T-coefficient exactly `g = 36771482407` (asserted in code). The T-free sublattice of the genuine lattice has **rank 8** (the 4 T-free LLL rows + combinations of the 5 T-rows); it contains the short exact-zero relations with atom max-coefs `3` (= `3K1+K3`), `16` (= `16K0−K2`), a ζ(5)-relation with coef `8952191640`, and one with `959163390`.
- **Atom-space minimisation** (LLL on the 8×74 T-free lattice + Babai nearest-plane + 3⁸=6561 enumeration; no further gain from a 5⁴×3⁴=50625 box + coordinate descent): vector **v_atom**, atom max coefficient **1550285430850157040** (1.55e18).
- **Expanded max-numerator minimisation** (LLL of the T-free lattice in the objective metric `O = 120·(275 Li5 + 166 products + ζ(5))`, transform recovered exactly via `U = R_C·C⁻¹`, det(U)=−1, `U·Omat == R` checked; Babai + enumeration + descent): vector **v_num**, expanded **max numerator 9070223133044173680** (9.07e18), lcm of all coefficient denominators **1**.

Both saved in the `d17_final.json` format with the convention line (product keys canonicalised as sorted multisets; no collisions occurred, 80 keys raw = 80 canonical):

| file | vector | T_coef | z5_coef | Li5 terms | products | max numerator | lcm | atom max coef |
|---|---|---|---|---|---|---|---|---|
| `d13_gcd.json` | v_atom (min atom coef) | **−13652994017342251 = −36771482407·13⁵** | 34510543973886225 | 143 | 80 | 13253867932032456368 | 1 | 1550285430850157040 |
| `d13_gcd_minnum.json` | v_num (min max numerator) | **−13652994017342251 = −36771482407·13⁵** | −7136597506206732780 | 143 | 80 | **9070223133044173680** | 1 | 7136597506206732780 |

So the cleanest exact formula in the capped (cap-275) pool is e.g.
`13652994017342251·L(χ13,5)/√13 = −(−7136597506206732780·ζ(5)…)` per the convention line — i.e. T-coefficient **36771482407·13⁵** (proved minimal for this lattice: the gcd is achieved). This is exact but **not** the clean `13⁵·T = …` form.

### Evidence table (everything below was printed by my code)

| check | result | precision/method | status |
|---|---|---|---|
| T-coef of v_atom and v_num | exactly ±36771482407 (asserted; = −13652994017342251 in the JSON convention) | exact integer arithmetic | ✔ |
| gcd factorisation | 36771482407 = 11·61·941·58237 | exact gcd + sympy.factorint | ✔ |
| EXACT symbol (weight-5, 5 letters) | **0 nonzero entries** for v_atom and for v_num | exact integer symbol from pool coords (vx, v1m) + product permutations | ✔ ZERO |
| S-unit factorisations (x and 1−x) | **143/143 checked, 0 bad** for both | exact QD arithmetic + `factor_sunit` | ✔ |
| residual, both vectors | **−6.31e-311** | 300 digits, `explicit_resid` (planner-algorithm Li5) + Hurwitz L, 13⁵-scaled | ✔ |
| residual of `d13_gcd_minnum.json` | **−1.1667e-1304** | **1300 digits**, independent JSON-driven verifier: `planner_li5` + own Hurwitz L (`Σ χ13(a)ζ(5,a/13)/13⁵`) + own letter logs | ✔ |
| residuals at 2200 digits (both files) and 1300 for `d13_gcd.json` | **launched but NOT captured** (background pid 99126, log `verify_gcd.log`, still running at turn end) | 2200 digits, same verifier | ⚠ pending — not claimed |
| max-numerator optimality | 3⁸=6561 + 5⁴·3⁴=50625 enumerations + coordinate descent found nothing below 9.07e18 | exhaustive local search on the rank-8 coset | ✔ (local optimum) |

## 2. Steps 2–4 — NOT completed in the turn budget

- **Step 2 (cap-350 pool)**: pool built and saved (`out/d13/candidates_cap350.json`, 350 smallest-height candidates: heights 74/142/105/29, 58 with x<0; first 275 verified index-by-index identical to ZQ5c's cap-275 pool). ZQ5c measured kernel dim 126, rank 224 for cap-350. **The symbol kernel, 2000-digit values, and the LLL runs (1000→1500→1800) were not run.** No statement can be made about the T-coefficient gcd in the cap-350 pool; in particular **no 13⁵-normalised relation was found** (the search was not executed).
- **Step 3 (cap-400)**: not run (moot without step 2).
- **Step 4 (best-relation JSON beyond the gcd relation)**: not applicable; the best verified relations remain the two gcd relations above. (`progress.log` has the step_01 line; `summary.md` rewrite was not done at turn end.)

## 3. Files and rerun commands (cwd = workdir; PY = `<project>/.venv/bin/python`)

- Producers (harness-saved step files): `step_28.py` (Bezout gcd + atom-space LLL/Babai/enum → `out/d13/gcd_vec_atom.json`), `step_31.py` (O-metric LLL + exact transform recovery + Babai), `step_32.py` (enumeration + descent → `out/d13/gcd_vec_num.json`), `step_34.py` (exact symbol 0, factorisations, resid@300), `step_35.py` (saved `d13_gcd.json`, `d13_gcd_minnum.json`), `step_37.py` (wrote `verify_relation.py`, launched verification), `step_38.py` (`out/d13/candidates_cap350.json`).
- Rerun verification: `nice -n 10 $PY verify_relation.py d13_gcd_minnum.json` (and `d13_gcd.json`); partial output already in `verify_gcd.log` (residual @1300 = −1.1667e-1304 for the minnum file).
- Next steps (not run): kernel for cap-350 (`w5pipe.kernel_and_comps` on `candidates_cap350.json`), values @2000d for indices 275–349 (reuse `vals_cap_2000.pkl` for 0–274), then LLL at 1000 digits with time extrapolation before 1500/1800 (each < 45 min; ZQ5c reference timings on a quieter machine: dim 74 @1800 = 79s, dim 106 @1800 = 346s, dim 152 @1500 = 757s; expected cap-350 dim ≈ 2 + effective + 40 ≈ 110–135).