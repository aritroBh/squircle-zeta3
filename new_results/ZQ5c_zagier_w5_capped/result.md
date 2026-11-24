# ZQ5c: weight‑5 formulas for L(χ_D,5), capped pools (D=29, D=13) — final report

All results are **NUMERICAL** (LLL-based, planner‑validated Li5), not proofs. Workdir: `new_results/ZQ5c_zagier_w5_capped`. Precision traps respected (`mp.mp.dps`, `gram='exact'`, Python ints).

## 1. Weight‑3 formulas (step 0) — DONE, exact JSONs saved

Re‑ran ZQ5b's weight‑3 pipeline, expanded kernel combinations to explicit terms, at **500 and 600 digits**, each form fresh‑verified at 700 digits; both precisions give valid forms differing only by exact‑zero atom relations (confirmed: difference residuals 1.03e‑732 (D=5), ‑6.52e‑731 (D=13) @700d). Exact symbol of each expanded form is **ZERO**.

| D | T_coef | ζ(3)_coef | Li3 terms | products | max numerator | lcm | fresh resid @700d | symbol |
|---|--------|-----------|-----------|----------|---------------|-----|-------------------|--------|
| 5 | **125 = 5³** | 0 | 38 | 11 | 1703 | 2 | 4.07e‑726 | ZERO |
| 13 | **2197 = 13³** | ‑4 | 257 | 39 | 5633 | 2 | ‑6.66e‑724 | ZERO |

- Pools: D=5: 50 exceptional S‑units (B=8, S={5,2}); D=13: 321 (h≤3, B=6, S={13,2,3}). Kernels 42 (D=5), 281 (D=13).
- Files: **`w3_D5.json`, `w3_D13.json`** (arguments `a+b*sqrtD` → rational coef, T_coef, z3_coef, products with `pi2*` keys, stats, residuals, symbol/agreement flags). Produced by `step_04_w3v2.py` (lib `w3lib_c.py`).

## 2. D=29, weight 5 — NONE FOUND (structural obstruction identified)

L(χ29,5) = 0.966085032430360658693721418661, T = L/√29 = 0.179397486804943938927849105668 (values @2000d; Li5 2nd‑way worst 6.68e‑2033; 12 duplication pairs worst 1.65e‑2039).

| S | pool | kernel | effective combos | LLL dim | precisions | result | reach |
|---|------|--------|------------------|---------|-----------|--------|-------|
| {29,2,5} | 77 (B=8; **provably exhausted**: B=10,12 give same 77) | 7 | **1** | 75 (π²·log³,π⁴·log) | 1500 (76s), 1800 (95s) | none common | 20 / 24 digits |
| {29,2,5,7} | 1088 (B=5) | h≤2:44, h≤4:78, full:106 | **1 in all cases** (K0 ≈ ‑0.507067) | 94 (D17 atoms, h≤2) | 1500 (161s), 1800 (220s) | none common | 16 / 19 |
| {29,2,5,7} | h≤2 | 44 | 1 | 152 (FULL atoms incl ζ(3)/T3) | 1500 (757s), 1800 (901s) | none common | 10 / 12 |

**Root cause:** for D=29 the zero‑symbol Li5 span collapses to a **single** effective combination K0 (robust across both S choices and all height caps), so T would have to equal c₀·ζ(5)+c₁·K0+products, and it does not.
**Wording:** *none found for D=29 with S={29,2,5} at 1800 digits (coefficients ≲ 1e19); none found for D=29 with S={29,2,5,7} at 1800 digits (coefficients ≲ 1e12).* (Not a statement that no formula exists.)

## 3. D=13, weight 5 — genuine exact relations found; clean 13^5 form not in capped pool

L(χ13,5) = 0.973328977797491433269523388413. Capped pool per instructions: 150→275 smallest‑height candidates to reach kernel ≥ 50 → **kernel dim 51, 32 effective combos**. Values @2000d verified (2nd‑way worst 6.64e‑2033; 42 duplication pairs worst 3.17e‑2039).

- LLL (D17 atom types, dim 74) at 1500 (59s) and 1800 (79s): **9 genuine relations** (residual at atom floor ~1e‑2008; the ~61 others at ~1e‑1746 are spurious): 4 T‑free + 5 with T. The two shortest **fresh‑verified EXACT at 2100 digits** (residuals 2.01e‑2107, ‑2.93e‑2108 — independent recomputation).
- **However** every genuine T‑relation has T‑coef a multiple of **36771482407 = 11·61·941·58237** (T‑coefs: 4155177511991, 4302263441619, 5184779019387, 19452114193303, 35594794969976). Minimal T‑coef in this span is therefore **36771482407·13⁵**, i.e. the cleanest exact formula here is
  `36771482407·13^5·L(χ13,5)/√13 = (55‑term combination, maxcoef ≈ 3.3e17)` — exact, but **not** the clean `13^5·T = …` form (analogue of D=5,13 weight‑3 and D=17 weight‑5).
- FULL atom list (dim 106, incl ζ(3)/T3 products) at 1500/1800 gave the **same** genuine relations — ζ(3)/T3 product atoms all get coefficient 0 in the genuine relations.
- Larger caps measured (kernel dims): **cap300→76, cap350→126, cap400→176, cap450→226, cap487→263** (rank 224 throughout). LLL on these bigger pools was **not run** (turn budget). The clean `13^5·T` relation, if it exists for S={13,2,3}, likely lives in one of these richer spans (dim would be ≤ ~160 for cap≤350, reach ~11–16 digits at 1500–1800).

## 4. Verification summary

| item | status |
|------|--------|
| w3 D=5: symbol ZERO; two‑prec forms differ by exact‑zero; fresh @700d 4.07e‑726 | ✔ |
| w3 D=13: symbol ZERO; two‑prec forms differ by exact‑zero; fresh @700d ‑6.66e‑724 | ✔ |
| d29 candidate factorizations (x and 1−x S‑units), 20/20 sign‑aware | ✔ |
| All Li5 values 2nd‑way cross‑checked (worst 6.7e‑2033); duplication relations in pools (worst 3.2e‑2039) | ✔ |
| D=13 capped genuine relations: two‑precision (1500/1800) + **fresh independent @2100d** (1e‑2107) | ✔ |
| D=13 clean 13^5 formula | ✘ not found in capped pool (min T‑coef 36771482407·13⁵) |

Note (step 4e): in D=13's genuine relations, all ζ(3)/T3 product atoms carry coefficient **0**.

## 5. Files and rerun commands (cwd = ZQ5c_zagier_w5_capped; PY = the venv python)

- Weight‑3: `w3_D5.json`, `w3_D13.json` — rerun `nice -n 10 PY step_04_w3v2.py` (libs `w5blib.py`, `w3lib_c.py`).
- D=29: data `out/d29/candidates_S7_B5.json`, kernels `out/d29/kernel5_S7h2.json`/`_S7h4.json`, values `out/d29/vals_S7_2000.pkl`, LLL out `out/d29/lll_S7h2*_P*.json` — rerun `step_05_d29_kernel.py`, `step_12_d29S7_vals.py`, `step_16_d29S7_lll.py`, `step_18_d29S7full_lll.py`.
- D=13: `out/d13/candidates_cap.json`, `kernel5_cap.json`, `vals_cap_2000.pkl`, **`genuine_lattice_cap.json`** (the 9 exact relations), `lll_capd17_P1500/P1800.json` — rerun `step_13_d13_cap.py`, `step_14_d13cap_kernel.py`, `step_19_d13cap_vals.py`, `step_20_d13cap_lll.py` (FULL=0/1), `step_24_d13_lattice.py`, fresh check `step_23_d13_fresh.py`. Expansion/verification helper: `w5verify.py`.
- Next step (not run): LLL on cap‑350 (kernel 126, dim ≈ 2+126+40 = 168) at 1800–2000 digits to look for the clean `13^5·T` relation; D=21 (S={3,7,2,5}) was not attempted since no clean D=29/D=13 relation was found.

`summary.md` and `progress.log` in the workdir carry this same state machine‑readably.