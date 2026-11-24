The current map, for draft 50 of the paper, is [VERIFICATION_MAP_draft50.md](VERIFICATION_MAP_draft50.md). This file describes the state at draft 10.

# Verification map for "Investigating ζ(3)"

This is the material that was moved out of Table II of the paper (the old Table III and its continuation): the *Method*,
*Re-check* and *Script* columns, and the list of folders. It is linked from the README of this repository and maps every row of Table II to its scripts.
Numbers in the *Result* column are those of the theorems, propositions, sections and tables in the IEEE copy of the
paper (`paper/ieee/main.tex`). The repository is public.

## Table A: results of Sections IX to XIV (the old Table III)

| Result | Status | Method | Digits (as in the text) | Re-check | Script |
|---|---|---|---|---|---|
| $M_1$; Thm. 1 | Proved | exact reduction over $\mathbb Q(\sqrt2)$ (App. A); two quadratures | exact; 1200 | 69; 1200 | `m1_symbolic_proof`, `mk_numerics` |
| Thm. 2 ($M_2$) | Proved† | exact reduction (App. B); two quadratures | exact; 1200 | 68; 1200 | `m2_weight3_proof` |
| Thm. 3 ($M_3=W_4+W_3'+T$) | Proved† | exact reduction; quadrature | exact | 68 | `m3_tails_exact` |
| Thm. 4 ($W_4-\mathcal R\in\mathbb Q\pi^4+\mathbb Q\sqrt2\pi^4$) | Proved† | symbol, coproducts (exact); Arb balls | $\|c_0\|<10^{-322}$, $\|c_1\|<10^{-321}$; zero to 1200 | certificate re-run; 1200 | `w4_symbol`, `w4_coproduct`, `w4_kernel`, `w4_certified` |
| Prop. 1 (weight two, $k\le5$); Step 3 for $k\le5$ | Proved | exact reducer | exact | exact re-run ($k=4,5$) | `ibp_reduce` |
| Lemma 1; values of $M_0,\dots,M_3$ | Proved; check | written proof; Cauchy integral on $\|\varepsilon\|=0.9$ | — | 41–43 | `verify_all mk` |
| Thm. 5 ($P'''(1)$, any rhombus) | Proved | log-Gamma expansion; finite differences, 5 aspect ratios | exact | 42–43 | `r2_analysis`; `verify_all rhombus` |
| Prop. 2 ($a_0,\dots,a_7$; $\omega_n$) | Proved | exact series; own generating function; formal expansion | exact; $\omega_n\ne0$ for odd $n\le21$ | 89–90; extra | `rhombus_exact_coefficients`; `verify_all borel` |
| Shape of the coefficients for any rhombus (Sec. X-C) | Proved | exact series | exact | orders $\le3$ (any $a,b$), $\le8$ (square) | `r2_analysis`; `verify_all borel` |
| Thm. 8 (radius 0) | Proved† | written proof; growth of coefficients | — | Borel radius $\pi/2$ to 5–6 | `verify_all borel` |
| Growth law of $f_n$ (Sec. XI-A) | Numerical | Richardson extrapolation, $n\le64$ | constant $1$ within $6\cdot10^{-6}$ | same | `verify_all borel` |
| Padé–Borel sum (Sec. XI-A) | Numerical | $[90/90]$ Padé of the Borel transform, $f_0,\dots,f_{180}$ at 450 digits | 44, 31, 25, 30 | re-run: 44.1, 30.8, 24.8, 29.6 | `pade_borel_check` |
| Thm. 9 (radius 1) | Proved | from Lemma 1 and Thm. 8 | — | — | — |
| Obs. 1 ($M_k$ decay; $\varpi_1$) | Numerical | $M_k$, $k\le400$, two methods; Cauchy integrals to $k=2400$; fits; two derivations of $\varpi_1$ | $\ge31$; $A_0=1\pm1.3\cdot10^{-6}$, $\phi_0=\pi(1\pm2\cdot10^{-7})$ | ratio $0.997$–$1.010$ at $k=150,200,250$ | `verify_all mk` (table to $k=400$) |
| Thm. 10 ($\pi_p$: $N_0,\dots,N_3$; $N_4,N_5$) | Proved†; $N_4,N_5$ modulo 5 Euler sums | series arithmetic in $1/p$; Euler sums to $\ge100$ digits | exact | 59–60 | `verify_all pip` |
| Thm. 11 (exact $L_\theta(n)$) | Proved† | vs. arc length, 7 angles $\times$ 8 values of $n$ | $\|\Delta\|\le1.6\cdot10^{-42}$ | — | `v2_exact_length` |
| Thm. 12 ($a_1,\dots,a_4$; remainder); Ahmed's Table IV | Proved† | Cauchy coefficients, 8 angles; exact formula at $n=10^6$ | $\|\Delta\|\le5.4\cdot10^{-71}$; all printed digits | extra | `v3_taylor_coeffs`, `v4_asymptotics`, `v8_ahmed_table` |
| Thm. 13 ($M_1(\alpha)$, all $\alpha$) | Proved | vs. quadrature, 6 angles (App. C) | $\|\Delta\|\le5\cdot10^{-38}$ | 35–40 | `corner_M1_theorem`; `verify_all cornerthm` |
| (octagon form), (pi/6 form), $M_1(\frac{2\pi}3)$ | Numerical | PSLQ/LLL; own quadrature | 500; about 189 for $\frac{2\pi}{3}$ | 68–70 | `verify_all corner` |
| Thm. 16 ($P'(1)$) | Proved† | integration by parts | exact | 18 | `hexagon_p1_check` |
| Thm. 18, its expansion, the slope-$2/e$ law and contact of order $m$ | Proved† | Richardson fit, 6 functions | 94–102 | — | `analyze`, `degenerate`, `contact_all_m` |
| Bounds for $A(n)$ (Sec. XIV-A) ($x^ny^n$) | Proved | bounds vs. own quadrature, $n\le10^6$ | bounds hold | 6 for $D(10^5)$; bounds hold | `area_xnyn`; `verify_all area` |
| Thm. 20 ($S(p)$, $Q_0$, $s_1$, $s_2$) | Proved† | exact reduction; fit of $S(p)$ | exact; 60, 50 | surface-area integral at $p=1,2$ | `13_proof_Q0`, `16_richardson_check`; `verify_all ellipsoid` |
| $s_3$ (value), $Q_1$ | Numerical (formula proved) | three routes for $Q_1$; fit of $S(p)$ | 40; $Q_1$: 50 (400 by one route) | — | `17_Q1_subpieces`, `16_richardson_check` |
| $S_n(p)$ to order $p^{-2}$ | Proved† | same argument; direct quadrature at $n=4$; independent $p^{-2}$ coefficient | $5\cdot10^{-32}$; $10^{-65}$ | — | `20_ndim_check`, `21_ndim_s2_check` |

† Proof sketched here; the complete written proof is in our research notes (Section XV of the paper). "Digits" is the precision stated in the text. "Re-check" is the agreement found independently by `verify_all.py` (the section named in the last column; for the first four rows, the sections `squircle`, `vform`, `hi` and `cert`) or, where it says "extra", by `extra_checks.py`. For Prop. 1 and the Padé–Borel row it is instead a re-run of the named script (for Prop. 1 through `k5_weight2_rerun` and `ipoles_k5_check`, which call `ibp_reduce`), not an independent check. A dash means the row is checked only by its own script. Scripts are named without the `.py` and live in the folders `squircle_square_end/scripts`, `weight4_M3`, `rhombus_end`, `corners_ahmed_P2`, `corners_p_smoothed`, `polygons_p1`, `graphs_fn_ahmed_P3`, `surfaces` (for `area_xnyn`) and `surfaces/superellipsoid` (for the numbered scripts) of the repository, or at its top level (`verify_all`, `extra_checks`). The new results are checked in Table B.

## Table B: results of Sections X to XIV added in the second round (the old "Table III (Continued)")

| Result | Status | Method | Digits (as in the text) | Re-check | Script |
|---|---|---|---|---|---|
| Thm. 6 (pure weight in $\tau=1-1/p$: rhombus; $\{\sum_i\|x_i\|^p=1\}$ in $\mathbb R^n$); 2D series to $\tau^{12}$ | Proved† | $\log\Gamma$ grading; three exact routes | exact | direct perimeter, $p=1.01$–$1.04$: through $\tau^8$ | `d16_engine`; `planner_check_tau` |
| Thm. 7 (Fourier coefficients of $h_p$ to $\tau^4$, every $k$; $P=\pi a_0(h_p)$); Prop. 3 (degree of $\Theta_n$) | Proved† | expansion of an exact integral; master identity and Watson coefficients | exact; forms: $10^{-50}$; $n\le12$ | Richardson quadrature, $k\le7$: $8.5\cdot10^{-51}$; D11's exact coefficients: $10^{-97}$ | `f1_01_duality_check`, `f1_03_diamond_fourier`, `step_12`, `step_16`; `v1_diamond_richardson` |
| $\tau^5$ term; closed forms $\Theta_5$–$\Theta_7$ | Numerical | PSLQ; $\tau^5$ for $k\le3$ | residual $10^{-94}$ ($\Theta_n$); $2.5\cdot10^{-60}$ ($\tau^5$) | quadrature, $k\le2$; $2.7\cdot10^{-50}$ ($\tau^5$, $k\le7$) | `step_14`, `step_29`, `planner_check_Q7`; `f1_03_diamond_fourier`, `v1_diamond_richardson` |
| Two-ends dictionary; $\varpi_1$ (Sec. XI-D); Prop. 4; $K_0(q)$ (Sec. XI-E) | Numerical (mechanism derived); Proved† (ruler); $K_0$ at $q=3,4$ numerical | Domb–Sykes fit of $f_n$; dictionary for $\varpi_1$; exact representation; PSLQ | radius $\pi/2$: $10^{-9}$; $\varpi_1$: 81; exact; 110 | $M_k=\frac{\sqrt2}{2}[\varepsilon^k]F(\varepsilon-1)$: $10^{-57}$–$10^{-22}$; own quadrature: 16–28 | step scripts of `D13c_two_ends`, `D11b2_stokes_jump`; `planner_check_angle` |
| Thm. 14; Thm. 15; Prop. 5 | Proved† | symbolic in $n$; nested harmonic sums for all $n$ | exact | own Cauchy contours, 13 angles: $2\cdot10^{-46}$; exact Taylor check, $n\le300$ | `flat_symbolic_proof`, `m2_w2_identity`, `b3_allN_proof`; `v_cauchy`, `v_taylor_exact` |
| $M_2(\alpha)$ with classical $\Li_3$ at eight arguments; three bases for $M_2(\frac{2\pi}3)$ | Numerical | identities checked numerically; PSLQ | 300 (9 angles); 300 | formula typed by hand, 13 angles: $6\cdot10^{-47}$; 320 | `m2_classical_assemble`, `report_formula_check`, `m2_2pi3_pslq`; `v_closed_forms`, `v_basis_change` |
| Thm. 17 ($C_m$); Prop. 6 (edge law; parallelogram and hexagon laws; hexagon $P''$, $P'''$) | Proved† | exact symbolic integration and Beta moments; independent jet code; finite differences | exact; $C_m$ for $m\le8$: $10^{-120}$; hexagon values first found by PSLQ at 170 digits | jet code, 3 polygons: $4.4\cdot10^{-55}$; new parallelogram: 13; hexagon: 17, 13 | `step_09_verify`, `s04_exact_local`, `t1_sanity`, `planner_check_parallelogram`; `vjet`, `verify_all hexagon` |
| More $\zeta(3)$ than the edge law for $\{\|x\|^p+\|y\|^p+3\|x-y\|^p\le4^p\}$ | Numerical | PSLQ per arc and on the total | residual $8\cdot10^{-141}$ | finite differences, 170 digits: all printed | `t3_hexagon`, `t7b_total150`; `t4_fd` |
| Corner–corner $\zeta(5)$ law $C_m^{(5)}$ (Sec. XII-C) | Numerical (exact at $m=2$) | exact integrals vs. extractor, $m\le8$ | $4\cdot10^{-58}$ | $-90\sqrt2$ from Prop. 2: exact | `step_28`, `step_33`; `pert_extract` |
| Thm. 19 (Ahmed's Problem 4: $A(n)$ through $n^{-3}$) | Proved† | decomposition, Watson's lemma, Mellin moments; free fits to own values, $n\le10^{20}$ | 47–63 | own 25-digit values, $n\le10^5$: through $n^{-3}$ | `fit_all`, `pslq`; `planner_check_constants` |
| Thm. 21 (octahedron end: $S(1+d)$ to $d^{12}$; single zeta values, every dimension) | Proved† | exact Gamma moments, two routes; direct quadrature | exact; $10^{-33}$ | 2-D quadrature: error $\propto d^5$ | `01_series_routeA`, `02_series_routeB`, `04_compare_series_vs_direct`; `planner_check_area` |

† and "Digits" as in Table A; "Re-check" is an independent computation, by another method or parametrization, done by the scripts after the semicolon. Scripts are named without the `.py`; they are in subfolders of `new_results` in the repository, named after the computation (`F1_fourier_squircle`, `K3_fourier_all_orders`, `D16_pure_weight`, `G1_zeta3_any_angle`, `W1_corner_family`, `D15c_edge_zeta3`, `E1_edge_law`, `D17_zeta5_polygons`, `P4c_area_xnyn`, `H1_octahedron_end`, `D13c_two_ends`, `D11b2_stokes_jump`, `D12c_ruler_family`).

## Notes for whoever assembles the repository

* The scripts of Tables B and C live in subfolders of `new_results/`.
* The `P4c` scripts `fit_all.py` and `pslq.py` are in `new_results/P4c_area_xnyn/verify_num/`.
* `python verify_all.py` re-derives, from its own integral representations, the quantities behind most entries of the *Re-check* column
  (about four minutes on a laptop); `extra_checks.py` covers the entries marked "extra".

## Table C: rows added in draft 6 (last six rows of Table II)

| Result | Status | Method | Digits | Re-check | Folder (in `new_results/`) |
|---|---|---|---|---|---|
| Level-8 relation at weight two; count $9=3^2$ | Numerical; exact count | LLL hint (325 digits), then Rogers five-term instance at $(\sqrt2/2,\sqrt2-1)$; exact rank over $\mathbb Q$ | residual $10^{-329}$ | planner re-evaluation, $2.8\cdot10^{-41}$ at 40 digits (`planner_check_w2.py`) | `A3f_level8_relations` (engine in `Q/`) |
| Bound for $M_k^{\rm main}$, $100\le N\le10^{10}$ | Computer-assisted | arb interval arithmetic, box assembly (`step_37.py`, `step_38.py`) | $C=914$ on $[100,400]$, $617118$ on $[400,10^{10}]$ | planner rerun of `step_38.py`; deterministic 43-cell sweep `planner_sweep.py` | `D13f_sweep`, `D13g_all_N` |
| Decay of $M_k(\alpha)$: rate $\sqrt{2\alpha N}$; amplitude, phase | Numerical; Conjecture | Cauchy/FFT tables $k\le300$ at 5 angles | 21–59 (tables) | least-squares fits `planner_check_fit.py`, `planner_check_fixed_lead.py` | `K2_corner_decay` |
| Prop. 7 (non-local $\zeta(3)$ of an edge) | Derived; not yet reviewed | exact symbolic reduction (J-table, fixed $E_3$) | 100 (two hexagon families) | `planner_check_law.py` (sympy, exact) | `K1d_engine_fixed` |
| Regular octagon: $P(1)$, $P'(1)$, $C_4$; no relation for $P'''(1)$ | Numerical | LLL over the S-unit basis of $\mathbb Q(\sqrt2)$ | 200; none up to $8.7\cdot10^{12}$ at 520 | `planner_check_closed_forms.py` | `OCT3_octagon_basis` |
| Cuboctahedron: $S(1+d)$ to $d^3$ | Numerical | tanh–sinh surface integral; integer relations after removing $4^{d/(1+d)}$ | 48–51 | `planner_check_pslq.py` | `K4b_cubo_finish` |

## Table D: rows added or changed in draft 7

| Result | Status | Method | Digits | Re-check | Folder (in `new_results/`) |
|---|---|---|---|---|---|
| Bound for $M_k^{\rm main}$, all $N\ge100$, $C=254.23$ | Computer-assisted; reviewed twice | arb box assembly (29 + 85,181 cells) + analytic final cell; bugs B1/B2 fixed | — | planner reruns of `step_18.py`, `step_10_J192.py`, 10 sweep cells; reviews D13r, D13t | `D13h_final_cell`, `D13s_fix_recertify` (fixed code in `D13s_work/fixed/`), `D13r_lemma_review`, `D13t_review_fixes` |
| Corner law with $\kappa_1(\alpha)$ | Derived; checked at 5 angles | Hankel + saddle analysis | tables $k\le300$ | `planner_check_kappa1.py` (residual decays as $N^{-1}$; conjugate fails) | `K2b_corner_phase`, `K2_corner_decay` |
| 3D series to $\tau^6$ | Proved; independent numerical check to $\tau^8$ | exact Gamma-product expansion (`d16_engine.py`, n = 3) | exact | `planner_run_n3.py`; GLM-5 $S(p)$ data, residual $\propto\tau^9$ (`planner_check_fit.py`) | `D16_pure_weight`, `GLM2_pureweight3d` |
| Tetrakis hexahedron, $d^2$ term | Numerical | kink-split tanh–sinh (`k4c_driver.py`), PSLQ | 58 (three grids) | `planner_check_T2.py` | `K4d_identify`, `K4c_3d_edge_law` |
| Octagon double polylog $\mathcal D$ (superseded in draft 8: $\mathcal D$ reduces, see Table E) | Numerical; one identity proved | letters of $E_3$; integral reduction | 40+ | planner quadrature (diff $4\cdot10^{-42}$) | `OCT4_octagon_letters` |
| Level-8 weight-3 relations; count 27 | Numerical (not proved) | Glanois basis $B^8$, exact LLL (sympy) at 800 digits | residuals $\le6\cdot10^{-781}$ | `planner_check_relations.py` | `A3f4_glanois_basis` |

## Table E: rows added or changed in draft 8

| Result | Status | Method | Digits | Re-check | Folder (in `new_results/`) |
|---|---|---|---|---|---|
| Prop. 8 (order-five edge law $Z_e^{(5)}$, $Z_e^{(2,3)}$; $C_m^{(5)}$, $C_m^{(2,3)}$; parallelograms, regular hexagon) | Proved$^\dagger$ | exact symbolic integration of the weight-5 integrand (`steps/step_01`–`step_05`); law minus raw integrals $=0$ | exact; 57 (sheared parallelogram, `steps/step_07.py`) | `planner_check_beta.py`: independent exact log-Beta expansion of the arc length (diamond, sheared parallelogram: every coordinate of $P^{(5)}(1)$ matches) | `E5_edge_law` |
| Regular octagon: $P'''(1)=48uS$ in closed form; $\mathcal D$ and its swap reduce | Numerical; $A$-atom reductions proved | exact assembly of the $d^3$ integrand, IBP for mixed integrals, LLL at 745 digits | 412 (independent quadrature); $\mathcal D$: 745 | `planner_check_closed_form.py` (reads `table_Stot.pkl`; vs the 200-digit value in `D15b_polygon_law/vals_m2_8_200.json`: $1.6\cdot10^{-196}$) | `OCT6_octagon_closed`, `D15b_polygon_law` |
| Tetrakis hexahedron, $d^2$ term to 85 digits | Numerical | kink-split tanh–sinh, $h=1/256$, dps 135 (`jobA_hp.py`) | 85 (finest of four grids) | `planner_check_T2_85.py` ($d^2$ residual $1.2\cdot10^{-85}$) | `K4e_tetrakis_d3` |

## Table F: rows added or changed in drafts 9 and 10

| Result | Status | Method | Digits | Re-check | Folder (in `new_results/`) |
|---|---|---|---|---|---|
| **Theorem 5: W₄ = 𝓡, so M₃ has a closed form** (Conjecture 4 proved) | Proved (computer-assisted) | exact certificate: 5068 relations (4810 standard, 230 lifts, 28 lifts of 13 frame double-shuffle relations), exact check on 5184 level-8 words, π⁴ coefficient 0 | exact | `A3h5_exact_certificate/verify_D.py` (runs from this repo, ~4 s, prints "Conjecture 4 PROVED"; rebuilds every row via `A3h4_frame_ds/frame_ds.py`, `A3c_claude_elementary/Q`); domain audit `A3h4_frame_ds/PLANNER_AUDIT.md` | `A3h5_exact_certificate`, `A3h4_frame_ds`, `A3c_claude_elementary`, `A3_W4_proof`, `A3f3_level8_compute` |
| ω_n ≠ 0 for every odd n ≥ 3 (Prop. 2) | Proved | parity of V_n (Hurwitz series) | exact (V_n odd checked to n = 61) | `CX2_review_omega_k1e/review.py` | `CX2_review_omega_k1e` |
| Theorem 19: P‴(1) of every centrally symmetric hexagon | Proved† | exact symbolic reduction + Nielsen–Kölbig identity | 3 independent direct checks (21–94 digits) | `K1e_hexagon_closed/planner_check_direct.py`, `planner_check_irrational.py`; `CX2_review_omega_k1e/review.py` | `K1e_hexagon_closed` |
| Octagon depth-2 reductions (𝒟 and swap) | Proved (two routes) | derivative families G(t), H(t) (exact); exact level-8 certificate with frame rows | exact | `PRV1_depth2_proofs/verify_full.py`; `OCT_R34_level8/certificate_frames_claude.py`; `CX1_octagon_depth2_proof/verify.py` | `PRV1_depth2_proofs`, `OCT_R34_level8`, `CX1_octagon_depth2_proof` |
| Regular decagon P‴(1), 15 terms | Numerical (depth-2 reductions and A-atoms proved) | series-engine integrand + exact reductions | 676 (quadrature), 197 (independent value) | `PRV1_depth2_proofs/verify_full.py` (reductions) | `DEC1_decagon`, `PRV1_depth2_proofs` |
| χ₃ bridge: 48χ₃(√2−1) = 32√2 L(χ₈,3) − 3π²L + 4L³ | Proved | exact level-8 certificate (470 standard rows) | exact | `SQ_chi3_bridge/verify_bridge_claude.py` | `SQ_chi3_bridge` |
| Weight-two part of M_k for k ≤ 15 | Exact computation | second exact reducer (fastred.py) | exact | `WT2_weight2_all_k/planner_check_k6_original.py` (original reducer at k = 6) | `WT2_weight2_all_k` |
| Level-8 weight three: R1–R5 new and independent | Proved (exact ranks over ℚ(√2)) | regular representation, flint rank | exact (544, 545, 549) | `A3h3_exact_w3/planner/planner_check_exact.py` | `A3h3_exact_w3` |
| Non-local ζ(5) on a hexagon edge | Numerical | integer relation over 20 alternating MZVs | 336 | `E5e_pslq_340` steps 08–12 | `E5e_pslq_340` |
| Generic octagon: depth 2 stays on one edge | Numerical (bound ~1e5) | LLL/PSLQ over 90 single-polylog constants | 450 search, 700 check | `OCT7_generic_octagon` steps | `OCT7_generic_octagon` |
