# new_results

The scripts, outputs and write-ups behind the later results of *Investigating ζ(3)* (Sections X to XIV of the paper). The main `README.md` of the repository maps these results to the paper; this file says what is in each folder and how far each result is established.

Numbers such as "Theorem 6" are those of the paper. Statuses follow the research ledger: **proved** (a written proof, checked), **numerical** (agreement to the stated number of digits, no proof) or **conjecture**. *Proved†* means that the full written proof is in the research notes, not in this repository.

## What is in a folder

| File | What it is |
|---|---|
| `result.md` | The write-up of the job that produced the result, as it was written. Later corrections are listed below and in the table, and they win over `result.md`. |
| `verify.md`, `verify_num.md`, `verify_proof.md`, `verify/` | A second reader re-checked the result with its own code. `D10c_corner_general/verify_proof.md` has the corrected proof that the paper uses. |
| `planner_check*.py` | A separate re-computation of key numbers, done with its own code and not with the job's scripts. |
| `step_NN.py` | Numbered working scripts that `result.md` cites. Short throwaway steps that no write-up cites were left out. |
| `*.out`, `*.json`, `*.tsv`, `*.txt` | Outputs of the runs, each under 2 MB. |

Left out on purpose: job transcripts, Python pickles, and the short scratch steps. Absolute paths in the copied text were replaced by `research/…` (the project tree) or `<scratch>/…`. The scripts are kept as they were run, so some read inputs that sit in the research tree and not here (a script that loads a pickle needs the step that wrote it to be re-run first). Each folder is meant to be run from inside itself, with the packages in `../requirements.txt`. The precision and the running time are stated in `result.md`.

## The folders

| Folder | Paper item | Status | Main files |
|---|---|---|---|
| `D16_pure_weight/` | Theorem 6: pure weight in τ = 1 − 1/p, for the rhombus and for the surface of Σ\|xᵢ\|^p = 1 in ℝⁿ; the 2D series is exact through τ¹² (three routes) | proved†; checked against the direct perimeter at p = 1.01–1.04 through τ⁸ | `d16_engine.py`, `planner_check_tau.py` |
| `D11c_zeta_first_order/` | Proposition 2: ω_n = (4(−1)ⁿ/n)(1 − 2^(1−n)V_n), the coefficient of ζ(n) in aₙ, so ω₃ = 1 and ω₅ = −3/4; the rhombus version; ζ(3) = π²/3 + P‴(1)/(6√2) | proved; ω_n ≠ 0 is checked for odd n ≤ 21, not proved for all n | `result.md`, `step_01.py` … `step_08.py` |
| `F1_fourier_squircle/` | Theorem 7: Fourier coefficients of the support function h_p at the diamond end, every k, through τ⁴ (and P = π·a₀(h_p)) | proved†; τ⁵ term numerical | `f1_01_duality_check.py` … `f1_06_comb_and_figure.py`, `verify/v1_diamond_richardson.py` (k ≤ 7), `planner_check_fourier.py` |
| `K3_fourier_all_orders/` | Proposition 3: the coefficient of τⁿ is a polynomial Θₙ(k²) of degree exactly ⌊n/2⌋ − 1 (called Qₙ in `result.md`); closed forms Θ₅–Θ₇ | degree proved†; closed forms numerical (PSLQ, residual 1e−94) | `planner_check_Q7.py`, the cited `step_NN.py` |
| `D13b_Mk_constant/` | Observation 1: M_k ≈ −(√π/k) e^(−√(πk)) cos √(πk), amplitude √π, phase π, and the first correction κ₁ in closed form | numerical: the table `mk_deg9.tsv` (k ≤ 400) and independent Cauchy-integral values up to k = 2400. The saddle step has no written error bound, so this is not a theorem | `mk_deg8.tsv`, `mk_deg9.tsv`, `result.md` |
| `D13c_two_ends/` | Section XI-D: the p = ∞ coefficients M_k follow from the Borel data at p = 1; the Stokes jump goes like ε^(1/2) (an earlier statement had ε^(3/2)) | numerical, mechanism derived | `result.md`, `d11b_util.py`, the cited `step_NN.py` |
| `D12c_ruler_family/` | Proposition 4: L₁(p) = 8 and L_∞(p) = 8·2^(−1/p) for every p. Section XI-E: K₀(q) for q ≥ 3 involves periods of the Fermat curve; radius ≥ 1 at every fixed ruler q, radius 0 for q ∼ λp | proposition proved†; K₀(3), K₀(4): numerical (PSLQ at 110 digits finds no classical relation); ζ(3) coefficient as a function of q: conjecture | `planner_check_angle.py`, `result.md` |
| `D10c_corner_general/` | Theorem 13: M₁(α) at every interior angle, with real dilogarithms | proved. The ζ(2) coefficient is not canonical (+1 in the real form, −1 in the first form of `result.md`) | `verify_proof.md` (corrected proof), `verify_num.md`, `planner_check_paper_thm.py`, `planner_check_quad.py`, `corner_series.py` |
| `G1_zeta3_any_angle/` | Theorem 14: M_k(α) = Σ μ_{k,n} c^(2n) about the flat angle, with the ζ(3) part of μ_{3,n} equal to −ζ(3)/(4n²−1); Theorem 15: M₂(α) = W₂ + W₃; Proposition 5; M₂(α) with Li₃ at eight arguments; three bases for M₂(2π/3) | theorems proved†; the Li₃ closed forms numerical (300 digits, 9 angles). Splitting M₃ into rational, ζ(2) and ζ(3) parts is conditional, because the ℚ-independence of 1, π², ζ(3) is open | `flat_symbolic_proof.py`, `m2_w2_identity.py`, `m2_classical_assemble.py`, `verify/v_cauchy.py`, `verify/v_taylor_exact.py`, `planner_check_W2.py` |
| `D15c_edge_zeta3/` | Theorem 17: C_m = 48m sin⁴h cos h (3 sin²h − 1), h = π/2m | proved† for the closed form of the local integral. C_m is the whole ζ(3) coefficient of P‴(1) only for m = 2, 3 (see the correction) | `step_09_verify.py`, `sym_G.py`, `poly_coeffs.py`, `z3coeff_m2_8.json` |
| `E1_edge_law/` | Proposition 6: Z_e = 6ξ₀ξ₁(h²/L)(3(ξ₀L₀ + ξ₁L₁)/L − 2), the parallelogram and affine-regular hexagon laws; the hexagon {\|x\|^p + \|y\|^p + 3\|x−y\|^p ≤ 4^p} | proved† for the laws. The general edge law is refuted (see the correction); the hexagon is numerical (PSLQ at 110 and 150 digits, finite differences at 60) | `s04_exact_local.py`, `planner_check_parallelogram.py`, `planner_check_counterexample.py`, `verify/t3_hexagon.py`, `verify/t7b_total150.py` |
| `D17_zeta5_polygons/` | Section XII-C: C_m⁽⁵⁾ = 2880m sin⁴h cos h (65 sin⁶h − 63 sin⁴h + 17 sin²h − 1) and C_m^(2,3) = −960m cos h sin⁴h (45 sin⁶h − 27 sin⁴h + 6 sin²h − 1) | numerical: matches the perturbative extractor for m = 2…8 to 4e−58; exact at m = 2 (−90√2 ζ(5) − 210√2 ζ(2)ζ(3)). It is the local part only | `closed_form_J5.json`, `pert_extract.py`, `pert_extract_m2_8.json`, `universal_weights.json` |
| `P4c_area_xnyn/` | Theorem 19: A(n) = 1 + √2/n + R(n)/n² + R₃(n)/n³ + O(log⁴n/n⁴) with closed forms for β, γ and d₃; ζ(3) enters at n⁻³ | proved† (adversarial proof check), and the constants agree with free fits to n ≤ 10²⁰ | `planner_check_constants.py`, `verify_num/fit_all.py`, `verify_num/pslq.py`, `verify_proof.md`, `verify_num.md` |
| `H1_octahedron_end/` | Theorem 21: S(1+d) = 4√3[1 + d + (1 − π²/9)d² + (1 + π²/9 − (16/9)ζ(3))d³ + …], exact through d¹², single zeta values in every dimension | proved†; Borel singularities at (−log(n−1) ± iπ)/2 numerical; zero radius for n ≥ 3 is a conjecture | `01_series_routeA.py`, `02_series_routeB.py`, `04_compare_series_vs_direct.py`, `planner_check_area.py` |

## Corrections that come after some `result.md` files

- **The edge law holds only for fully local polygons.** `D15c_edge_zeta3/result.md` and `E1_edge_law/result.md` (its first part) say that the ζ(3) coefficient of P‴(1) is Σ_e Z_e for every centrally symmetric polygon. That is proved only when every term is local: every parallelogram and every affine-regular hexagon, so C₂ = 6√2 and C₃ = −9√3/8 are complete. For K_p = {\|x\|^p + \|y\|^p + 3\|x−y\|^p ≤ 4^p} the ζ(3) coefficient is (231/125)√5 − (15/64)√2, while Σ_e Z_e = (72/125)√5 − (15/64)√2 (`E1_edge_law/planner_check_counterexample.py` computes P‴(1) = 19.90026636314805, which agrees with the verifier's closed form to 2e−14). The terms that mix a corner with a far kink line carry no ζ(3) one by one, but their sum can. In general Σ_e Z_e is only the corner–corner part.
- **Regular 2m-gons with m ≥ 4.** C_m is the local part. Whether it is the whole ζ(3) coefficient of P‴(1) is open, and the "cross = total" statement of `D15c_edge_zeta3/result.md` is checked only for m = 2, 3. The size is C_m = −3π⁴/m³ + 25π⁶/(8m⁵) − …, not 9π⁵/(16m⁴).
- **D17.** `D17_zeta5_polygons/result.md` is a stub: the job ended without a final answer. Its findings are in `closed_form_J5.json` and `pert_extract_m2_8.json`. The same caveat applies: whether C_m⁽⁵⁾ is the whole ζ(5) coefficient of P⁽⁵⁾(1) is open for m ≥ 3, and exact only for the diamond.
- **D13c.** The Stokes jump goes like ε^(1/2); an earlier statement of it had ε^(3/2).
- **D10c.** The paper's version of the proof of Theorem 13 is the "Corrected proof" in `verify_proof.md`. In `result.md`, the claim that only ±½, ±1, ±2, ±c/2, ±c, ±2c occur as dilogarithm coefficients is false (−4c, 1+c and 2−2c also occur).
- **G1.** In a fixed angle the number M₂(α) can show ζ(3) in one basis and not in another (at 2π/3 it is 3/16, 123/16 or −29/4 in three valid bases). Always say "as a power series in c² about the flat angle".

## Not here yet

The paper's verification map also names `W1_corner_family/` (it holds `b3_allN_proof.py`, the all-n nested harmonic sums behind Theorem 14 and Proposition 5) and `D11b2_stokes_jump/`. They are not in this commit. Verifier scripts for `D10c_corner_general` and `P4c_area_xnyn` that exist only in scratch space (apart from `fit_all.py` and `pslq.py`) are not here either.
