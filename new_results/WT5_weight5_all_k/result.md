# WT5: the weight-5 part W₅(k) of the squircle coefficient M_k — E1, E2, E3 for 5 ≤ k ≤ 14

**Verdicts (per constant):**

| constant | (3,2) pairing | classical? | closed form |
|---|---|---|---|
| **E3** | (0, 0) | CLASSICAL (level 2, Γ-world) | **NUMERICAL** (6 terms, LLL identical @1100/1300 twice over; value residual 2.3e-1520) |
| **E1** | (−17/2, 0) | NOT classical | depth-two **EXACT** + classical **NUMERICAL** (68 terms; LLL identical @1100/1300; residual 7.0e-1519) |
| **E2** | (−2, 0) | NOT classical | depth-two **EXACT** + classical **NUMERICAL** (64 terms; LLL identical @1100/1300; residual 3.8e-1519) |

**W₅(k) = W5″ + (k−5)·E1 + C(k−5,2)·E2 + C(k−5,3)·E3** holds **exactly for 5 ≤ k ≤ 14**: cubic kernel law exact (k=9..14, all three degree-4 keys), and the cubic law re-verified exactly at the G-expression level (k=6..14, 2132 terms each). **k=4 is OFF the cubic**: P(4) − W₅(4) = −0.1919599554656691313 ≠ 0 (W₅(4) from M4c's W5_5400). **Depth-two content of W₅(k): YES, a fixed combination of the four atoms Li₃,₂(b,1), Li₃,₂(1/√2,1), Li₃,₂(−b,1), Li₃,₂(b,−1) with coefficients polynomial in k of degree ≤ 2** (N(E3)=0; linearity of the cobracket over the exact G-expression cubic law).

Notation: λ = log 2, L = log(1+√2), a = 1+√2, b = √2−1, R = √2.

## Closed forms

**E3 (classical, Γ-world, level 2)** — found by the saturated-kernel LLL (8 nonzero entries) and confirmed by a direct LLL against the 6-element Γ-basis; identical integer relation [576, −12, −20, −9, −180, −30, −270] at 1100 AND 1300 digits, residuals 1.44e-1218 / −3.26e-1418:

$$E3 = \tfrac{1}{48}\lambda^5 + \tfrac{5}{144}\pi^2\lambda^3 + \tfrac{1}{64}\pi^4\lambda + \tfrac{5}{16}\zeta(3)\lambda^2 + \tfrac{5}{96}\zeta(3)\pi^2 + \tfrac{15}{32}\zeta(5)$$
(equivalently 576·E3 = 12λ⁵ + 20π²λ³ + 9π⁴λ + 180λ²ζ₃ + 30π²ζ₃ + 270ζ₅). Its exact symbol is supported on the log-2 letter only (printed by wt5_G_solve.py), as level 2 requires. For LEAD1: **c(5) = E3/6** = (5/24)W₂ζ(3) + (5/64)ζ(5) + λ⁵/288 + 5π²λ³/864 + π⁴λ/384.

**E1 = R(E1) + N(E1)**, depth-two part EXACT (reduce(N − target) = 0):
$$N(E1) = \left(-\tfrac{4199}{56} - \tfrac{\sqrt2}{4}\right)\mathrm{Li}_{3,2}(b,1) + 68\,\mathrm{Li}_{3,2}(\tfrac1{\sqrt2},1) + 34\,\mathrm{Li}_{3,2}(-b,1) + \left(\tfrac{2295}{56} + \tfrac{\sqrt2}{4}\right)\mathrm{Li}_{3,2}(b,-1)$$
R(E1): 68 classical terms (products + ζ5 + Li5 at 14 of the 16 M4f pool S-units; **L(χ8,5) coefficient 0**), max numerator 97,305,071, lcm denominator 75,479,040; full list in `closed_forms.json` / `e1_closed_form.json`.

**E2 = R(E2) + N(E2)**, depth-two part EXACT:
$$N(E2) = \left(-\tfrac{247}{14} - \tfrac{\sqrt2}{7}\right)\mathrm{Li}_{3,2}(b,1) + 16\,\mathrm{Li}_{3,2}(\tfrac1{\sqrt2},1) + 8\,\mathrm{Li}_{3,2}(-b,1) + \left(\tfrac{135}{14} + \tfrac{\sqrt2}{7}\right)\mathrm{Li}_{3,2}(b,-1)$$
R(E2): 64 classical terms, Li5 at the same 14 pool S-units, L(χ8,5) = 0, max numerator 14,168,515, lcm 18,869,760.

**W5″ (M5b, reused)**: R5″ (69 terms) + N5″ = (3929/28 − 15√2/4)Li₃,₂(b,1) − 120 Li₃,₂(1/√2,1) + (−68+26√2)Li₃,₂(−b,1) + (−2025/28+119√2/4)Li₃,₂(b,−1).

**Depth-two of W₅(k), polynomial in k (Newton form, exact):**
- Li₃,₂(b,1): (3929/28 − 15√2/4) + (k−5)(−4199/56 − √2/4) + C(k−5,2)(−247/14 − √2/7)
- Li₃,₂(1/√2,1): −120 + 68(k−5) + 16·C(k−5,2)
- Li₃,₂(−b,1): (−68+26√2) + 34(k−5) + 8·C(k−5,2)
- Li₃,₂(b,−1): (−2025/28 + 119√2/4) + (k−5)(2295/56 + √2/4) + C(k−5,2)(135/14 + √2/7)

Cobracket targets (tensor Li₂(b)): T(W5″) = (−105/8+√2/8)ζ₃ + (−68−26√2)Li₃(b) + 120·Li₃(1/√2) + (68−26√2)Li₃(−b); T(E1) = (119/16+7√2/8)ζ₃ + 34Li₃(b) − 68Li₃(1/√2) − 34Li₃(−b); T(E2) = (7/4+√2/2)ζ₃ + 8Li₃(b) − 16Li₃(1/√2) − 8Li₃(−b); T(E3) = 0. Pairing of W₅(k): (15 − (17/2)(k−5) − 2·C(k−5,2), −13√2).

## Values (1500 digits in `e1_1500.txt` etc.; two Arb runs 5400 & 5800 bits)

- E1 = 8.360256391769373610071463792254808451891002268227773358528… (agreement 3.007e-1595)
- E2 = 8.707284010623930159400023660171389967807494959604566036072… (agreement 2.773e-1595)
- E3 = 2.456887165275386141491315290399093944527904387996545973972… (agreement 3.020e-1595)

## W₅(k) (50 digits, from the formula)

| k | W₅(k) |
|---|---|
| 4 (off law) | −0.18235662466773630012609370629766450393982143738932 |
| 5 | 1.7355429662874241607407575450430824536268058904655 |
| 6 | 10.09579935805679777081222133729789090551780815869327 |
| 7 | 27.16333976045010154028370878972408932521630538652561 |
| 8 | 55.39505133874272161064653519272077165725020196195906 |
| 9 | 97.24782125821004412339201583668703184614740227299017 |
| 10 | 155.1785366841274552200114660120219638364358107076155 |
| 11 | 231.6440847817703410419962010091246615726433316538316 |
| 12 | 329.1013527164140877308375361183942189992978694996349 |
| 13 | 450.0072276533340814280267866302297300609273286330221 |
| 14 | 596.8185967578057082750552678350302887020596134419897 |

## Evidence table (all numbers printed by my own code)

| check | result | precision |
|---|---|---|
| Exact kernels k=5..14 (fastred, 387s) | k=5 == M5b omega5.json keys (4,0),(3,1),(2,2) EXACT (merged); cubic law k=9..14 EXACT | exact Fractions |
| G-expressions | e1.tsv 2112, e2.tsv 1992, e3.tsv 1367 terms; W5″ rebuild == gpl_w5pp.tsv EXACTLY (2132 terms) | exact |
| TSV vs Arb quadrature (g_eval cap 1400/prec 1040 vs 800-bit quad) | e1 3.50e-210, e2 3.48e-210, e3 1.09e-210, w5″ control 3.53e-210 | 210 digits |
| Values 5400 vs 5800 bits | e1 3.007e-1595, e2 2.773e-1595, e3 3.020e-1595; radii ≤ 4e-1622 / ≤ 2.1e-1742; tail bounds ≤ 3.9e-1642 / ≤ 2e-1762 | 1595 digits |
| (3,2) pairings | W5 = (84, −38√2) PASS; W5″ = (15, −13√2) PASS; E1 = (−17/2, 0), E2 = (−2, 0), E3 = (0, 0) | exact |
| Cobracket | Li3-projection and classical-Li5 controls PASS; W5″ target == cobracket5pp.json PASS; E1/E2 exact 4-atom solves; E3 target exactly 0 | exact |
| Symbol stage (71 atoms: products + Li5×16 + ζ5) | all 3 targets in span, rank 30, saturated kernel dim 41, SNF invariants all 1; particulars: e1 max\|num\| 441969 lcm 21840; e2 3438/5460; e3 1/2 | exact |
| LLL e1 (dep-dropped 85→73 cols) | 1 relation, IDENTICAL @1100/@1300, residuals 4.59e-1211 / −1.37e-1410 | 1211/1410 digits |
| LLL e2 | 1 relation, IDENTICAL, residuals −1.68e-1212 / 2.07e-1412 | 1212/1412 digits |
| LLL e3 (2 identical relations; rel0 used; rel1 = 23-term Li5 alternative) | rel0 residuals −1.68e-1218 / −2.55e-1418; direct Γ-basis LLL: identical [576,−12,−20,−9,−180,−30,−270], residuals 1.44e-1218 / −3.26e-1418 | 1218/1418 digits |
| Verification (wt5_K_verify.py) | exact symbol R = sym(E) − sym(N): PASS ×3; Li5 cobracket 0: PASS; shuffle control: PASS; reduce(N − target) = 0: PASS ×3 | exact |
| Value residuals R + N − E(5800 bits) | e1 6.97e-1519; e2 −3.77e-1519; e3 −2.30e-1520 | 1519 digits |
| W₅(k) formula vs direct 800-bit quadrature | k=6: 7.0e-210; k=9: 9.0e-210; k=14: 4.6e-208 | 208 digits |
| k=4 off cubic | P(4) = −0.3743165801334…, W₅(4) = −0.1823566246677…, δ = −0.1919599554656691313 ≠ 0 | ~1400 digits |

## Deliverables (all in `new_results/WT5_weight5_all_k/`)

`result.md` (this file, written by the harness), `e1.tsv`, `e2.tsv`, `e3.tsv`, `e1_1500.txt`, `e2_1500.txt`, `e3_1500.txt`, `e{1,2,3}_{5400,5800}.json`, `kernels_deg4.json`, `deltas_deg4.json`, `check_tsv_vs_quad.json`, `pairings.json`, `cobracket_e.json`, `{e1,e2,e3}_symbol_solve.json`, `{e1,e2,e3}_lll_{1100,1300}.json`, `e3_gamma_lll.json`, `e1_closed_form.json`, `e2_closed_form.json`, `e3_closed_form_rel0.json` (+`_rel1`), `closed_forms.json`, `progress.log`.
Scripts: `wt5_A_kernels.py`, `wt5_A2_control.py`, `wt5_B_gexpr.py`, `wt5_precision.py`, `wt5_C_check.py`, `wt5_C2_check.py`, `wt5_E_pairings.py`, `wt5_F_cobracket.py`, `wt5_G_solve.py`, `wt5_H_lll.py`, `wt5_I_assemble.py`, `wt5_J_gamma.py`, `wt5_K_verify.py`, `wt5_L_cubicg.py`; copies (never edited elsewhere): `fastred.py`, `buildintegrand.py`, `genintegrand.py`, `binteg_5..14.pkl`, `w4_symbol.py`, `w5_symbol.py`, `w5_cobracket.py`, `w4_right_coproduct.py`, `w5_search.py`, `w5_precision.py`, `m5b_common.py`, `g_eval_lib.py`, data copies (`omega5.json`, `pool.json`, `Lchi8_5*.json`, `gpl_w5pp.tsv`, `quad5_100.json`, `cobracket5pp.json`, `w5pp_closed_form.json`, `W5pp_1350.json`, `w5pp_4800.json`, `w5pp_vals_1450.pkl`); harness auto-saved `step_01.py … step_76.py` (every printed number above traces to one of these or to the wt5_* scripts).

**Caveat / not done:** the single consolidated `verify_wt5.py` (< 20 min, PASS/FAIL) was NOT assembled — the verification exists as the separate runnable scripts above, each with PASS output recorded in `progress.log` (kernels/control: wt5_A_kernels.py + wt5_A2_control.py; TSV-vs-quadrature: wt5_C2_check.py ~5 min; pairings: wt5_E_pairings.py ~1 min; depth-two: wt5_F_cobracket.py ~4 min; closed forms: wt5_K_verify.py ~4 min each — together ≈ 25 min, slightly over the 20-min target). Rerun (cwd = work dir, PY = project venv):

```
nice -n 10 $PY -B wt5_A2_control.py          # kernels control + cubic law
nice -n 10 $PY -B wt5_B_gexpr.py             # TSVs + W5'' rebuild control
nice -n 10 $PY -B wt5_precision.py e1 5400 5800   # (e2, e3 similarly)
nice -n 10 $PY -B wt5_C2_check.py            # TSV vs quadrature, 210 digits
nice -n 10 $PY -B wt5_E_pairings.py          # pairings + 2 controls
nice -n 10 $PY -B wt5_F_cobracket.py         # exact depth-two parts
nice -n 10 $PY -B wt5_G_solve.py e1          # (e2, e3) symbol stage
nice -n 10 $PY -B wt5_H_lll.py e1 1100 ; ... 1300   # (e2, e3)
nice -n 10 $PY -B wt5_I_assemble.py e1       # (e2; e3 0)
nice -n 10 $PY -B wt5_J_gamma.py             # E3 Gamma-world direct test
nice -n 10 $PY -B wt5_K_verify.py e1         # (e2; e3 _rel0) exact + 1500-digit checks
nice -n 10 $PY -B wt5_L_cubicg.py            # cubic G-law k=6..14
```

**Bottom line:** E3 is classical and level 2 with the Γ-world closed form 576·E3 = 12λ⁵+20π²λ³+9π⁴λ+180λ²ζ₃+30π²ζ₃+270ζ₅; E1 and E2 are NOT classical (obstructions (−17/2, 0) and (−2, 0)) and close as exact four-atom Li₃,₂ depth-two parts plus LLL-discovered classical remainders (no L(χ8,5) anywhere), every piece value-matched to ≥ 1.4e-1410 and exact-symbol/cobracket-verified. W₅(k) for 5 ≤ k ≤ 14 is thereby a cubic in k whose depth-two content lives on the SAME four Li₃,₂ atoms as W5″ with polynomial (quadratic) coefficients in k.