import os, json
work='<project>/research/08_theory/WT4_weight4_all_k'
md = r'''# WT4: the weight-4 part W_4(k) of the squircle coefficient M_k — closed form for every k (Kimi K3)

## Verdict

| constant | status | closed form |
|---|---|---|
| D2 (second difference) | **PROVED** (exact certificate, standalone verifier prints `VERDICT: D2 = R PROVED`) | level 2: −λ⁴/24 − π²λ²/24 − π⁴/160 − ζ(3)λ/4 |
| D1 (slope) | **PROVED** (same method, `VERDICT: D1 = R PROVED`) | level 8, 16 terms, **no √2, no Li4/Li2** (below) |
| W4″ = W_4(4) | PROVED earlier (M4b2_w4pp_nice) | see M4b2 |
| W_4(k) = W4″ + (k−4)D1 + C(k−4,2)D2 for 4 ≤ k ≤ 14 | **PROVED** (kernel polynomial law exact for k ≤ 14, KSPAN + re-verified here) | assembled below |
| k = 3 | **does NOT lie on the quadratic**: P(3) − W_4(3) = 0.15446753106687400637871015148911002803587706268747… ≠ 0 (both sides exact/proved, difference 800+ digits) | — |

The all-k statement (every k ≥ 4) still needs: the kernel polynomial law K^(k) = K4 + (k−4)Δ1 + C(k−4,2)Δ2
proved for *all* k (verified exactly here only for k ≤ 14) — NOT this job's task.

## The closed forms (λ = log 2, L = log(1+√2), b = √2−1)

**D2** = −λ⁴/24 − π²λ²/24 − π⁴/160 − ζ(3)λ/4.  Level-2 (alternating MZV of weight 4); Li₄(1/2) coefficient 0.
   Remark (exact rational identity, `wt4_G_assemble.py`): with W₂ = (π²+6λ²)/24 and S = ζ(3)/8 + π²λ/24 + λ³/12,
   D2/2 = c₂ = W₂² − λS − (7/16)ζ(4).

**D1** = −(5/48)λ⁴ − (1/12)λ³L − (1/4)λ²L² + (1/6)λL³ + (1/3)L⁴ − (7/96)π²λ² + (5/24)π²λL − (5/24)π²L² − (1/192)π⁴
        − (3/8)ζ(3)λ − (7/16)ζ(3)L − λLi₃(b) − 2L·Li₃(b) − λLi₃(−b) + 2L·Li₃(−b) + 4L·Li₃(1/√2).
   Level 8 but purely rational (no √2, no Li₄, no Li₂(b)); note λ+2L = log(2a²), λ−2L = log(2b²).

**W_4(k)** = W4″ + (k−4)·D1 + C(k−4,2)·D2, with W4″ as proved in M4b2 (combination of Li₄ at
{b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2} and 22 products). Equivalently W_4(k) = c₀ + c₁k + c₂k² with
c₂ = D2/2, c₁ = D1 − 9D2/2, c₀ = W4″ − 4D1 + 10D2 (wait: c₀ = W4″ +4D1·(−1)... see values JSON;
the anchored form above is the proved one).

## How it was done (controls in parentheses)

1. **Exact kernels** (`wt4_A_kernels.py`): fastred.reduce_k_fast on binteg_3..14; degree-3 keys (2,1),(3,0).
   (CONTROL: k=4 kernels == M4a omega4.json **exact**.) Deltas Δ1 = K⁽⁵⁾−K⁽⁴⁾, Δ2 = K⁽⁶⁾−2K⁽⁵⁾+K⁽⁴⁾
   (`deltas_deg3.json`): all coefficients rational (q=0), poles {0, a, −a, −b}.
   **K^(D2)_(2,1) = Rw2**, the weight-2 kernel — the structural reason D2 is level 2.
   (CONTROL: quadratic law K⁽ᵏ⁾ = K4+(k−4)Δ1+C(k−4,2)Δ2 **exact** for k = 7..14, both keys.)
2. **G-expressions** (`wt4_B_gexpr.py`): d1.tsv (420 terms), d2.tsv (388 terms), alphabet {0,±a,±b}.
   (CONTROL: rebuilt W4″ G-expression == M4a gpl_w4pp.tsv **exactly**, 432 terms.)
3. **TSV vs quadrature** (`wt4_C_check.py`, flint ctx.prec=700 vs mpmath dps 130):
   D1 |G−quad| = 2.474e-94, D2 = 1.651e-94, W4″ control = 6.01e-95 (flint radii ~1e-208);
   control quad(W4″) − M4a quad_90 = 3.005e-90 (the 90-digit reference limit).
4. **860-digit values** (`wt4_D_hp.py`, tanh-sinh, two working dps 900/850, err_est 1e-908):
   D1: |@900−@850| = 7.742e-812; D2: 1.498e-810. Files d1_860.txt, d1_810.txt, d2_860.txt, d2_810.txt.
5. **LLL closed forms** (`wt4_E_lll.py`, flint lll gram='exact', identical relations at 700 and 780 digits):
   D2 in the level-2 basis (D2 coef 480, max|coef| 480, residuals 3.0e-730 / 1.2e-809, check −3.18e-861);
   D1 in the full M4a 60-number basis (D1 coef −192, max|coef| 768, residuals 1.0e-730 / 1.9e-809, check 1.65e-861;
   16 terms, no √2, no Li4/Li2).
6. **Exact certificates** (M4b2 pipeline adapted: `wt4_cert_step01/02/04/06/07.py`, common `wt4_cert_common.py`):
   stage-1 mod-p solve at 8 primes, CRT×7 + hold-out (hold-out agrees, heights ≤ 3 digits) — the solved c_k
   **equal the LLL closed forms 58/58 components** for both targets; stage-2 exact fmpz solve + Fraction
   verification on ALL 5184 level-8 weight-4 words, residual 0. Both word vectors have **√2-part identically 0**,
   and **no frame lifts are needed** (all z_j = 0): D2 = 4 atoms + T0 rows (4906 used), D1 = 16 atoms + T0 (5030 used).
   Standalone fresh verifiers: `verify_wt4_cert.py D2` → `VERDICT: D2 = R PROVED`; `verify_wt4_cert.py D1` →
   `VERDICT: D1 = R PROVED` (each rebuilds T0 relations, frame rows, atom vectors, and the target from the TSV
   cross-checked against deltas_deg3.json, then checks the exact identity with Fractions).
7. **Sanity** (`wt4_F_sanity.py`, mp.mp.dps=950): R_D1 − D1(860d) = 1.646e-861; R_D2 − D2(860d) = 3.176e-861;
   W_4(3) quadrature (130 dps) − thm:W4 closed form = 1.468e-130; formula vs direct quadrature for k = 5..8:
   max |diff| 2.4e-108 (quadrature limit).

## Simpler proof for D2?

D2's G-expression uses letters {0, ±a, ±b} (level 8) — it does **not** reduce to letters {0, ±1}. But its proof
is materially simpler than W4″/D1's: (i) the (2,1)-kernel of D2 is **exactly Rw2** (the weight-2 kernel);
(ii) its level-8 word vector is purely rational (√2-part zero); (iii) its certificate needs only the 4 classical
level-2 atoms and the standard T0 relations (regularized double shuffle + distribution at level 8) — **no frame
double-shuffle lifts**; (iv) its value is a level-2 (alternating MZV) constant. So: the proof lives entirely in
the level-8 MZV quotient with classical atoms.

## k = 3

W_4(3) = 0.79578505983224081132035981962972676254652741387102… (quadrature of k=3 kernels == thm:W4 closed form
to 1.5e-130). Quadratic prediction P(3) = W4″ − D1 + D2 = 0.95025259089911481769906997111883679058240447655849….
δ = P(3) − W_4(3) = **0.15446753106687400637871015148911002803587706268747…** ≠ 0: **k = 3 is NOT on the
k ≥ 4 quadratic** (both sides exact constants — W_4(3) proved in A3h5/thm:W4, P(3) from proved pieces — so the
inequality is rigorous).

## W_4(k) values (50 digits, from the proved formula; k=5..8 confirmed by direct quadrature)

| k | W_4(k) |
|---|---|
| 3 (off law) | 0.79578505983224081132035981962972676254652741387102 |
| 4 | 1.4244536760088729215144224912572705356370998476755 |
| 5 | 0.87435084219557566488537389418208741787647316846362 |
| 6 | −0.70005591054077695218807582010671256269947556107721 |
| 7 | −3.298766582200184929705926651609129406090746340947 |
| 8 | −6.9217811727826482676681786003251631122973391711457 |
| 9 | −11.569099682288166966074831666254813681319254051673 |
| 10 | −17.24072211071674102492588584939808111315649098253 |
| 11 | −23.936648458068370444221341149754965407809049963715 |
| 12 | −31.65687872434305522396119756732546656527693099523 |
| 13 | −40.401412909540795364145455102109584585560134077073 |
| 14 | −50.170251013661590864774113754107319468658659209245 |

## verify_wt4.py (standalone, 14/14 PASS, ~2 min)

[1] k=4 kernels == omega4.json (exact); [2] quadratic law k=7..10 exact (11..14 in wt4_A); [3] TSV rebuilds exact;
[4] TSV vs quadrature dps 60 (D1 6.2e-61, D2 1.7e-60, W4″ 2.3e-60); [5] closed forms vs quadrature 200 dps
(2.9e-198, 5.7e-198); [6] certificates D2/D1 → PROVED; [7] W_4(3) quad == thm:W4 (1.7e-119); [8] δ(k=3) ≠ 0.

## Files

- Kernels/deltas/G-expressions: `kernels_deg3.json`, `deltas_deg3.json`, `d1.tsv`, `d2.tsv`
- Values: `d1_860.txt`, `d1_810.txt`, `d2_860.txt`, `d2_810.txt`, `hp_results.json`
- Closed forms: `closed_forms.json`, `lll_results.json`
- Certificates: `D1_certificate.json`, `D2_certificate.json`, `verify_wt4_cert.py`,
  pipeline `wt4_cert_common.py`, `wt4_cert_step01.py`, `wt4_cert_step02.py`, `wt4_cert_step04.py`,
  `wt4_cert_step06.py`, `wt4_cert_step07.py`; pickles `stage1_{D1,D2}_vectors.pkl`, `stage2_atoms_nice.pkl`,
  `stage2_{G,y}_{D1,D2}.pkl`, `stage1_exact_cz_{D1,D2}.json`
- Checks: `check_tsv_vs_quad.json`, `sanity_k3_table.json`, `verify_wt4.py`
- Steps: `wt4_A_kernels.py`, `wt4_B_gexpr.py`, `wt4_C_check.py`, `wt4_D_hp.py`, `wt4_E_lll.py`,
  `wt4_F_sanity.py`, `wt4_G_assemble.py`, `verify_wt4.py` (auto-saved harness steps step_NN.py alongside)

Proof caveats: same as M4b2/A3h5 — correctness of the cmzv and frame_ds conversion code (audited, numerically
validated to ~1e-114 per atom and ~1e-205 per identity there); T0 = regularized double shuffle + distribution
at level 8. Shared inputs: A3h5 `rref_T0_*.pkl`, `stage2_rows_S.pkl`, `stage1_exact_z.json`; A3h4 frame pickles;
A3c relations; M4a `gpl_w4pp.tsv`/`omega4.json` (controls only).
'''
open(os.path.join(work,'result.md'),'w').write(md)
open(os.path.join(work,'progress.log'),'a').write('result.md final written; verify_wt4.py 14/14 PASS\n')
print('result.md written,', len(md), 'chars')
