# M5b — W5″ closed (verified); W6−N6: symbol stage complete, kernel LLL at 2300 digits in progress

**Verdict W5″: NOT classical — it needs a depth-two N5″, and W5″ = R5″ + N5″ is fully closed and verified.**
**Verdict W6: symbol stage complete (rank 60, saturated kernel 89); the value-level kernel LLL did not close within the budget — largest LLL tried and finished: 1900 digits (spurious rows only, relation coefficients exceed ~5·10¹¹); a 2300-digit LLL on 2400-digit values is running now (log `w6_lll_2300.out`).**

## W5″ (planner-verified independently)

- **Precision**: Arb quadrature with log-endpoint substitutions + certified tail bounds (`m5b_precision.py`; the a-pole terms `c/(v−a)` are cancelled analytically on the right piece, `e·Ω(v) = −c_a + e·Ω̃`). Runs 4600/4800 bits agree to **1.8e-1382**; saved `W5pp_4600/4800.json` (radius ≤ 3.2e-1442, tail ≤ 3.2e-1462) and `W5pp_1350.json`.
- **G-expression** `gpl_w5pp.tsv` (2132 terms, weight check len+k=5): |G − quadrature| = **4.0e-100** (flint radius 1.1e-237) — `g_eval5pp.py`.
- **(3,2) pairings** (`w5pp_symbol.py`): **(15, −13√2)** — nonzero, so a depth-two part is required. CONTROL: same code gives W5 = (84, −38√2) ✓.
- **N5″** (`w5pp_cobracket.py`, M4c machinery; exact, residual asserted zero):
  N5″ = (3929/28 − 15√2/4)·Li₃,₂(b,1) − 120·Li₃,₂(1/√2,1) + (−68+26√2)·Li₃,₂(−b,1) + (−2025/28 + 119√2/4)·Li₃,₂(b,−1).
- **Classical part R5″**: symbol stage (71 atoms = products + Li5 at the 16 M4f-pool S-units + ζ5): rank 30, target in span, saturated kernel dim 41 (HNF, SNF all 1). Kernel LLL at **1100 and 1300 digits: identical relation** (residuals −3.5e-1211, 1.3e-1410; 12 exact dependencies removed first). Result: **69 classical terms, 14 Li5 arguments** {−b, 1−√2/2, b, 1/2, 2−√2, √2/2, −√2/2, b², −b/2, 2b, (2∓√2)/4, −b², 1−b⁴}, **max numerator 35,514,863, lcm denominator 37,739,520, L(χ8,5) coefficient 0**.
- **Verification** (`w5pp_verify.py`): 1. exact symbol of R5″ = symbol(W5″)−symbol(N5″) PASS; 2a. every Li5 → 0 in the (3,2) cobracket PASS; 2b. product control PASS; 2c. exact reduce(N5″−target)=0 PASS; 3. value R5″+N5″−W5″(4800 bits) = **−5.2e-1418**.
- **Independent direct-series check** (`m5b_indep_check.py`, no shared evaluator, 390 working digits): R5″+N5″ − quadrature = **5.0e-388**.

## W6 − N6 — how far it got

- **Precision**: W6 at 5600/5800 bits agree to **1.3e-1683** (radii ≤ 2.2e-1743, tails ≤ 2.6e-1760), matching `W6_400.txt` to 3.7e-401; extended to **6800 bits** (radius 1.9e-2044) and **8000 bits** (radius 1.4e-2405, 2378 digits). `W6_5400.json` (1650 digits) saved.
- **Exact symbol stage** (`w6_solve.py`): 149 atoms (weight-6 products incl. weight-5 factors = Li5 at the 10 R5_min3 args, ζ5, L(χ8,5), and the 4 Li₃,₂ atoms × λ, L; plus Li6 at W4's seven args + −b/2 + (2±√2)/4). Target = symbol(W6) − symbol(N6) (N6 words (0^{m−1},1/x,0^{n−1},1/(xy))). **Rank 60, target in span**; rref particular max|num| = 666,193, lcm 51,840; **saturated kernel dim 89** (HNF 11.3 s, SNF invariants all 1, LLL-reduced basis max entry 576). Saved `w6_symbol_solve.json`.
- **Values** (`w6_values.py`): all leaves at 2400 digits (cached `w6_leaves_2400.pkl`, `w6_atomvals_2400.pkl`). Controls: Li6 duplication (1/√2 pair vs 1/2, b pair vs b²) ≤ 3.5e-2511; mpl_general Li₃,₂(b,1) vs M4c's mpl = 3.7e-2441; N6 atoms Li₃,₃(b,1), Li₃,₃(1/√2,1), Li₃,₃(b,2), Li₄,₂(b,1/√2) vs M5a's values ≤ 4.7e-77 (their precision); Li₆ acb vs mpmath polylog 1.9e-82. L(χ8,5): the project χ8 is the **Q(√2) character (+1 on n≡±1 mod 8)**; my Hurwitz-ζ formula cross-checked against the planner's file to **4.7e-1811**, extended to 2500 digits (`Lchi8_5_2500.json`).
- **Kernel LLL** (`w6_lll_final.py`, 179 columns = target + 89 kv + 89 √2·kv, gram='exact'; planted-relation control PASS; dependency hunt at 500 digits found none):

| lattice | time | best rows (spurious) | genuine relation |
|---|---|---|---|
| 1300 digits | 546–1152 s | maxcoef ~2·10⁸, residual ~1e-1294 | none |
| 1600 digits | 850 s | maxcoef ~1·10¹⁰, residual ~1e-1592 | none |
| 1900 digits | 1274 s | maxcoef ~5·10¹¹, residual ~1e-1891 | none |
| **2300 digits** | **running (pid `w6_lll_2300.pid`)** | det^{1/n} ≈ 2.7e13, values at 2400 digits | — |

The relation was not isolated at ≤1900 digits because n=179 forces det^{1/n} = 10^{(D+103)/179}; genuine coefficients evidently exceed ~5·10¹¹. The 2300-digit run (≈40 min) should isolate coefficients up to ~10¹³. **No extra Li6 pool arguments have been needed so far** (extras −1/√2, b², 2b, b⁴, (4−3√2)/8, 1−b⁴ remain in reserve; `w6_solve.py <pool indices...>` adds them).

## Remaining recipe (all scripts written and syntax-checked)
```
cd new_results/M5b_w6_closed_form
PY='<project>/.venv/bin/python'
tail w6_lll_2300.out                          # wait for the relation (residual ~1e-2400 vs spurious ~1e-2290)
nice $PY -B w6_lll_final.py 1750              # second precision; then edit w6_assemble.py to the two jsons
nice $PY -B w6_assemble.py                    # R6 + stats (term count, Li6 args, max num, lcm)
nice $PY -B w6_verify.py                      # exact symbol, exact (3,3)/(4,2) pairings = 0, value >=1600 digits
nice $PY -B m5b_indep_check.py                # independent direct-series 300+ digit check (W5'' part already 5.0e-388)
```

## Key files (work dir `new_results/M5b_w6_closed_form/`)
- Scripts: `m5b_precision.py`, `w5pp_symbol.py`, `g_eval_lib.py`, `g_eval5pp.py`, `w5pp_cobracket.py`, `m5b_common.py`, `w5pp_solve.py`, `w5pp_lll.py`, `w5pp_assemble.py`, `w5pp_verify.py`, `w6_solve.py`, `w6_values.py`, `w6_lll.py`, `w6_lll_final.py`, `w6_assemble.py`, `w6_verify.py`, `m5b_indep_check.py` (+ copies of M4c/M4e templates: `w5_symbol.py`, `w5_cobracket.py`, `w5_search.py`, `w5_precision.py`, `nice_search.py`, `w4_symbol.py`, `w4_right_coproduct.py`).
- Data: `W6_5400.json`, `w6_5600/5800/6800/8000.json`, `W5pp_1350.json`, `w5pp_4600/4800.json`, `gpl_w5pp.tsv`, `symbol5pp.json`, `cobracket5pp.json`, `w5pp_symbol_solve.json`, `w5pp_lll_1100/1300.json`, `w5pp_closed_form.json`, `w6_symbol_solve.json`, `w6_leaves_{1650,2000,2400}.pkl`, `w6_atomvals_*.pkl`, `Lchi8_5_2100.json`, `Lchi8_5_2500.json`, `w6_lll_final_{1300,1600,1900}.json`, `w6_lll_2300.out`.
- progress.log documents every step (one line per finished step).