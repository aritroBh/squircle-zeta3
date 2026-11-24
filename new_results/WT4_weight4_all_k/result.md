# WT4: the weight-4 part W_4(k) of the squircle coefficient M_k — exact closed form for every k

**Verdict: D1 PROVED, D2 PROVED ⇒ W_4(k) = W4″ + (k−4)·D1 + C(k−4,2)·D2 is a PROVED closed form for every integer k with 4 ≤ k ≤ 14 (W4″ was proved in M4b2; the kernel polynomial law is exact for k ≤ 14). k = 3 does NOT lie on the quadratic (rigorous: δ ≠ 0 with both sides exact constants).**

All work is in `new_results/WT4_weight4_all_k/`. Notation: λ = log 2, L = log(1+√2), b = √2−1, a = 1+√2.

## The proved closed forms

**D2** (level 2; Li₄(1/2) coefficient 0):
$$D2 = -\tfrac{1}{24}\lambda^4 - \tfrac{1}{24}\pi^2\lambda^2 - \tfrac{1}{160}\pi^4 - \tfrac{1}{4}\zeta(3)\lambda$$
Remark (exact rational identity, checked in `wt4_G_assemble.py`): with W₂ = (π²+6λ²)/24 and S = ζ(3)/8 + π²λ/24 + λ³/12, D2/2 = c₂ = W₂² − λS − (7/16)ζ(4).

**D1** (level 8, 16 terms, but **no √2, no Li₄, no Li₂(b)**; both level-8 word-vector components: rat only, √2-part identically 0):
$$D1 = -\tfrac{5}{48}\lambda^4 - \tfrac{1}{12}\lambda^3 L - \tfrac{1}{4}\lambda^2 L^2 + \tfrac{1}{6}\lambda L^3 + \tfrac{1}{3}L^4 - \tfrac{7}{96}\pi^2\lambda^2 + \tfrac{5}{24}\pi^2\lambda L - \tfrac{5}{24}\pi^2 L^2 - \tfrac{1}{192}\pi^4 - \tfrac{3}{8}\zeta(3)\lambda - \tfrac{7}{16}\zeta(3)L - \lambda\,\mathrm{Li}_3(b) - 2L\,\mathrm{Li}_3(b) - \lambda\,\mathrm{Li}_3(-b) + 2L\,\mathrm{Li}_3(-b) + 4L\,\mathrm{Li}_3(1/\sqrt2)$$

**W_4(k)** = W4″ + (k−4)D1 + C(k−4,2)D2 (W4″ as proved in M4b2_w4pp_nice). Equivalently c₀+c₁k+c₂k² with c₂ = D2/2, c₁ = D1 − 9c₂, c₀ = W4″ − 4D1 + 10D2.

## How it was proved (pipeline, all controls passed)

1. **Exact kernels** (`wt4_A_kernels.py`): fastred.reduce_k_fast on binteg₃..₁₄; degree-3 keys (2,1),(3,0). CONTROL: k=4 kernels == M4a `omega4.json` **exact**. Deltas (`deltas_deg3.json`): D1: K₂₁ = 1/v + (3/4)/(v−a) − 2/(v+a) − (3/4)/(v+b), K₃₀ = (−1/12)/v − (1/2)/(v−a) + (5/6)/(v+a) − (1/6)/(v+b); D2: K₂₁ = ½/v + ½/(v−a) − 1/(v+a) − ½/(v+b) = **exactly Rw2 (the weight-2 kernel)** — the structural reason D2 is level 2 — and K₃₀ = (−1/3)/(v−a) + (1/3)/(v+a). Quadratic kernel law **exact for k = 7..14**, both keys.
2. **G-expressions** (`wt4_B_gexpr.py`): `d1.tsv` (420 terms), `d2.tsv` (388 terms), alphabet {0,±a,±b}. CONTROL: rebuilt W4″ == M4a `gpl_w4pp.tsv` **exactly** (432 terms).
3. **TSV vs quadrature** (`wt4_C_check.py`, flint ctx.prec=700, mpmath dps 130): D1 |G−quad| = 2.474e-94, D2 = 1.651e-94, W4″ control = 6.01e-95 (flint radii ~1e-208; mp.mp.dps printed = 130).
4. **860-digit values** (`wt4_D_hp.py`, tanh-sinh, err_est 1e-908): D1@900 vs D1@850 agree to 7.742e-812; D2 to 1.498e-810. Files `d1_860.txt`, `d1_810.txt`, `d2_860.txt`, `d2_810.txt`.
5. **LLL discovery** (`wt4_E_lll.py`, flint lll `gram='exact'`, identical relations at 700 and 780 digits): D2 in the level-2 basis (coef 480, max 480, residuals 3.0e-730/1.2e-809, closed form − value = −3.18e-861); D1 in the M4a 60-number basis (coef −192, max 768, residuals 1.0e-730/1.9e-809, check 1.65e-861; 16 terms, no √2/Li4/Li2).
6. **Exact certificates** (M4b2 pipeline adapted: `wt4_cert_common.py`, `wt4_cert_step01/02/04/06/07.py`): stage-1 mod-p solve at 8 primes, CRT×7 + hold-out prime (agrees; heights ≤ 3 digits); the solved coefficients equal the LLL forms **58/58 components**. Stage-2 exact fmpz solve + Fraction check on **all 5184 level-8 weight-4 words: residual 0**, both components, both targets. **No frame lifts needed** (z = 0): D2 = 4 atoms + 4906 T0 rows; D1 = 16 atoms + 5030 T0 rows. Standalone fresh verifiers (rebuild T0 relations, frame rows, atoms, and the TSV target cross-checked against `deltas_deg3.json`): `verify_wt4_cert.py D2` → `VERDICT: D2 = R PROVED`; `verify_wt4_cert.py D1` → `VERDICT: D1 = R PROVED`.
7. **Sanity** (`wt4_F_sanity.py`, dps 950): R_D1 − D1(860d) = 1.646e-861; R_D2 − D2(860d) = 3.176e-861; formula vs direct quadrature k = 5..8: max 2.4e-108 (quadrature limit); W_4(3) quadrature − thm:W4 closed form = 1.468e-130.

## Evidence table

| check | residual | digits / precision |
|---|---|---|
| k=4 kernels vs omega4.json | 0 | exact |
| quadratic law k=7..14 (both keys) | 0 | exact |
| G-expr rebuilds vs TSVs / gpl_w4pp.tsv | 0 | exact |
| TSV-vs-quadrature D1 / D2 / W4″ | 2.474e-94 / 1.651e-94 / 6.01e-95 | dps 130, flint rad ~1e-208 |
| D1, D2 two-precision quadrature | 7.742e-812 / 1.498e-810 | 860/810 digits |
| LLL D2 lvl2 (700/780) | 3.0e-730 / 1.2e-809; form−value −3.18e-861 | identical relations |
| LLL D1 60-basis (700/780) | 1.0e-730 / 1.9e-809; form−value 1.65e-861 | identical relations |
| Stage-1 CRT hold-out (both) | exact mod p; c_k == LLL 58/58 | 8 primes |
| Stage-2 exact verify (both, all 5184 cols) | 0 | exact Fractions |
| Standalone certificates | residual_nonzero=0 both components | exact; `VERDICT: D1/D2 = R PROVED` |
| R_D1, R_D2 (950 dps) vs 860-digit values | 1.646e-861 / 3.176e-861 | 861 digits |
| W_4(3) quadrature vs thm:W4 | 1.468e-130 | dps 130 / 950 |
| W_4(k) formula vs quadrature, k=5..8 | ≤ 2.4e-108 | dps 110 |
| δ = P(3)−W_4(3) ≠ 0 | δ = 0.15446753106687400637871015148911002803587706268747… | ~810 digits, rigorous inequality |
| verify_wt4.py standalone | **14/14 PASS** | ~2.5 min |

## k = 3 and the all-k statement

W_4(3) = 0.79578505983224081132035981962972676254652741387102… (quadrature == thm:W4 to 1.5e-130). Quadratic prediction P(3) = W4″ − D1 + D2 = 0.95025259089911481769906997111883679058240447655849…. **δ = 0.15446753106687400637871015148911002803587706268747… ≠ 0: k = 3 is NOT on the k ≥ 4 quadratic** (W_4(3) proved in A3h5/thm:W4, P(3) from proved pieces, so the inequality is rigorous). Bonus (NUMERICAL only, `delta_k3_lll.json`): δ is a 46-term relation in the same 60-number basis (LLL identical at 700/780 digits, candidate − δ = 2.33e-849, unnormalized coefficients up to 11520·469).

What an all-k statement still needs (NOT this job's task): the kernel polynomial law K⁽ᵏ⁾ = K4 + (k−4)Δ1 + C(k−4,2)Δ2 proved for **all** k ≥ 4 (here verified exactly for k ≤ 14; KSPAN verified the rank/degree structure for k ≤ 14).

## W_4(k) values (50 digits, from the proved formula)

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

## Simpler proof for D2?

D2's G-expression uses letters {0,±a,±b} (level 8) — it does **not** reduce to letters {0,±1}. But its proof is materially simpler: (i) K^(D2)_(2,1) = Rw2 exactly (weight-2 kernel); (ii) its word vector is purely rational (√2-part zero); (iii) its certificate needs only the 4 classical level-2 atoms and the standard T0 relations (regularized double shuffle + distribution at level 8) — **no frame double-shuffle lifts**; (iv) the value is an alternating-MZV (level-2) constant.

## Deliverables (all in `new_results/WT4_weight4_all_k/`)

`result.md`, `d1.tsv`, `d2.tsv`, `kernels_deg3.json`, `deltas_deg3.json`, `d1_860.txt`, `d1_810.txt`, `d2_860.txt`, `d2_810.txt`, `hp_results.json`, `closed_forms.json`, `lll_results.json`, `delta_k3_lll.json`, `check_tsv_vs_quad.json`, `sanity_k3_table.json`, `D1_certificate.json`, `D2_certificate.json`, `verify_wt4_cert.py`, `verify_wt4.py` (+`verify_wt4.out`, 14/14 PASS), certificate pipeline `wt4_cert_common.py`/`wt4_cert_step01.py`/`wt4_cert_step02.py`/`wt4_cert_step04.py`/`wt4_cert_step06.py`/`wt4_cert_step07.py`, pickles `stage1_{D1,D2}_vectors.pkl`, `stage2_atoms_nice.pkl`, `stage2_{G,y}_{D1,D2}.pkl`, `stage1_exact_cz_{D1,D2}.json`, step files `wt4_A_kernels.py`, `wt4_B_gexpr.py`, `wt4_C_check.py`, `wt4_D_hp.py`, `wt4_E_lll.py`, `wt4_F_sanity.py`, `wt4_G_assemble.py`, `wt4_H_delta_lll.py` (plus harness-auto-saved `step_NN.py`), `progress.log`.

Proof caveats: same as M4b2/A3h5 — correctness of the audited cmzv/frame_ds conversion code; T0 = regularized double shuffle + distribution at level 8. Shared inputs: A3h5 `rref_T0_*.pkl`/`stage2_rows_S.pkl`/`stage1_exact_z.json`, A3h4 frame pickles, A3c relations; M4a `gpl_w4pp.tsv`/`omega4.json` (controls only).