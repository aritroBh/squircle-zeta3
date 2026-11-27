# Investigating ζ(3): the code

This repository holds the code for *Investigating ζ(3)* (Aritro Bhattacharjee; mentor Dr. Bob Hingtgen, UCSC COSMOS Cluster 1, 2024). It is one project, not several.
The original version of *Investigating ζ(3)* ended on a non-linear homotopy that was still in progress: the squircle |x|^p + |y|^p = 1, whose perimeter P(p) runs from the rhombus (p = 1) to the square (p → ∞).
Everything here continues that thread, and each folder is one part of it.
At the square end we write P(p)/8 = 1 + Σ M_k p^(−k−1): the scripts find M₁ and M₂ in closed form, reduce M₃ to an explicit weight-4 integral, and pin that integral down up to a term (c₀ + √2 c₁)π⁴ with c₀, c₁ rational.
At the rhombus end, ζ(3) shows up in P‴(1), and the scripts compute that coefficient for general semi-axes a, b.
The same squircle tools then go into three of Ahmed's open problems (power-law corners at any angle, the length of the graph of fⁿ, and the area of z = xⁿyⁿ, now through the n⁻³ term) and into the first terms for the superellipsoid |x|^p + |y|^p + |z|^p = 1.
Further parts: the exact coefficients and the divergence of the expansion at the rhombus end, the decay of the M_k, π_p (the squircle measured with its own ℓᵖ ruler, where ζ(3)/4 appears at order p⁻⁴), smoothed corners of any angle, and smoothed polygons at p = 1.
The later results (the pure-weight theorem, Fourier series of the squircle, the polygon and ζ(5) laws, the area of xⁿyⁿ to order n⁻³, the octahedron end, and more) are in `new_results/`; numbers such as "Theorem 13" are those of the paper.
Every result in the table below comes with the script that checks it. The *Check* column says whether that check is exact (a proof by exact arithmetic) or numerical (agreement to the stated number of digits).

## State of the repository

The code here is as of draft 58 of the paper (45 pages).

- `new_results/` now has 238 folders, 2 more than before. They cover the later results of the paper.
- `proofs_supplement_scripts/` holds the check scripts that the written proofs supplement of the paper runs. The supplement itself is not here.
- `VERIFICATION_MAP_draft58.md` is the current map from every theorem and every row of Table II to its folders, reports and scripts. `VERIFICATION_MAP.md` is the older map, from draft 10.
- `OMITTED_LARGE_FILES.md` lists the files that were left out because they are too large.

Most scripts were written to run inside the research tree, so some need their paths adjusted. `verify_all.py` and the certificates named in the map run from this layout.

## The parts

All of these sit under *Investigating ζ(3)*:

| Folder | Part of the paper |
|---|---|
| `squircle_square_end/` | The squircle at the square end (p → ∞): M₀ … M₃, the exact ℚ(√2) polylogarithm engine, and the stored high-precision values |
| `weight4_M3/` | The weight-4 piece W₄ of M₃: symbol and coproduct checks |
| `rhombus_end/` | The squircle at the rhombus end (p = 1), where ζ(3) appears in P‴(1) |
| `corners_ahmed_P2/` | Ahmed's Problem 2: power-law corners at any angle, and the link back to the squircle corner constant |
| `graphs_fn_ahmed_P3/` | Ahmed's Problem 3: the length of the graph of fⁿ |
| `surfaces/` | Going up a dimension: the area of z = xⁿyⁿ (Ahmed's Problem 4) and the superellipsoid |
| `corners_p_smoothed/` | Smoothed corners of any interior angle α: the closed form of M₁(α) (Theorem 13) |
| `polygons_p1/` | Smoothed regular 2m-gons at p = 1: the P′(1) law and ζ(3) at the hexagon corners |
| `new_results/` | Scripts, outputs and write-ups for the later results, one folder per computation (see *New results* below) |
| `verify_all.py` | One command that re-checks the main results independently (see below) |
| `extra_checks.py` | Checks added in review round 1: Faulhaber's formula, the coefficient of ζ(n) in aₙ, Ahmed's Table IV from the exact corner formula, the corner-asymptotics remainder, the oblique-corner length (`python extra_checks.py`, seconds) |

## Setup

```sh
python3 -m venv .venv && . .venv/bin/activate      # Python >= 3.11
pip install -r requirements.txt
```

Each script uses paths relative to its own folder, so `cd` into that folder before running it.
Reference output from the runs cited in the paper is in `squircle_square_end/outputs/`, in the `*.out` files next to the scripts in `corners_ahmed_P2/` and `graphs_fn_ahmed_P3/`, and in the `out_*.txt`, `outputs_13_*.txt` and `outputs_25_proof_Zb.txt` files in `surfaces/superellipsoid/`.

## One command

```sh
python verify_all.py              # about 2 to 4 minutes on a laptop; prints one PASS/FAIL line per check
python verify_all.py rhombus pip  # or only some sections
```

`verify_all.py` re-derives the main quantities from its own integral representations (it imports nothing from the other folders) and compares them with the closed forms of the paper: the square-end coefficients M₀…M₄ and the 1200-digit stored values, the rhombus end (all aspect ratios, the coefficients a₀…a₈, the growth law), the decay of M_k, π_p, the corner and polygon results, the area of z = xⁿyⁿ and the superellipsoid representation. Its section `cert` re-runs the certified ball-arithmetic bound `weight4_M3/w4_certified.py`. Sections: `squircle vform lemmas rhombus borel mk area corner hi cert pip ellipsoid cornerthm hexagon`. The long exact computations (the 1200-digit values, the exact ℚ(√2) proofs) have their own scripts, listed below.

## Results and the scripts that check them

| Result | Folder: script | Command | Check |
|---|---|---|---|
| **Square end** | | | |
| M₀…M₃ to 1200 digits, computed by two independent parametrisations (A, B) plus a black-box large-p fit (C) | `squircle_square_end/scripts/mk_numerics.py` | `python mk_numerics.py A 1200` (likewise `B 1200` and `C 60`; about 37 min each at 1200 digits; `A 50` takes seconds) | numerical: A and B agree to all 1200 digits; C agrees with them to 3.7e−42 for M₃ |
| Closed form of M₁ (Theorem 1; κ = 8M₁ answers Ahmed's Open Problem 1) | `squircle_square_end/scripts/m1_symbolic_proof.py` | `python m1_symbolic_proof.py` | exact |
| The weight-2 part of M₂ is (π² + 6 log²2)/24 | `squircle_square_end/scripts/m2_decomposition_and_tail.py` | `python m2_decomposition_and_tail.py` | exact |
| Closed form of M₂ (Theorem 2: ζ(3), Li₃(√2−1), Li₃(1−√2), Li₃(1/√2), plus products) | `squircle_square_end/scripts/m2_weight3_proof.py` (uses `polylog_exact.py`, `kummer.py`, `dilog_exact.py`) | `python m2_weight3_proof.py` | exact |
| Which Li₃ identity reduces each Li₃ value (Lemma trilog) | `squircle_square_end/scripts/li3_table.py` | `python li3_table.py` | exact |
| After the Euler substitution, no simple poles at ±i (k ≤ 4); only the letters 0, ±(1+√2), ±(√2−1) remain | `squircle_square_end/scripts/ibp_reduce.py` | `python ibp_reduce.py 2 3 4` | exact |
| The same for k = 5, and the weight-2 part of M₄ and M₅ is (π² + 6 log²2)/24 (Proposition 1, k ≤ 5; re-runs of `ibp_reduce`, not independent) | `squircle_square_end/scripts/ipoles_k5_check.py`, `k5_weight2_rerun.py` | `python ipoles_k5_check.py 2 3 4 5`, `python k5_weight2_rerun.py` (about 2 min each; output in `../outputs/`) | exact |
| M₃ = W₄ + W₃′ + (π² + 6 log²2)/24, with W₃′ explicit (Theorem 3) | `squircle_square_end/scripts/m3_tails_exact.py` | `python m3_tails_exact.py 3` (and `4` for M₄) | exact |
| W₄ by direct quadrature of its integral equals M₃ − W₃′ − T | `squircle_square_end/scripts/w4_direct_check.py` | `python w4_direct_check.py 120` | numerical, 120 digits |
| The Li₄ closed form 𝓡 of W₄, found by LLL (Theorem 4 says W₄ − 𝓡 = (c₀ + √2 c₁)π⁴; c₀ = c₁ = 0 is numerical) | `squircle_square_end/scripts/m3_pslq.py`, `m3_closed_form.py`, `m3_assemble.py` | `python m3_pslq.py ../outputs/M3_methodA_1200.txt 1150`, then `python m3_closed_form.py` and `python m3_assemble.py` | numerical: the same relation from methods A and B at 700 and 1150 digits; agreement to 1202 digits |
| L(χ₈,3) and Cl₃(π/4) forms of M₂ (these rest on an identity that is still unproved) | `squircle_square_end/scripts/bridge_check.py` | `python bridge_check.py 1200` | numerical only, 1200 digits |
| **Weight 4 (M₃)** | | | |
| Full symbol of W₄ equals the symbol of 𝓡, the right-hand side of Theorem 4 (16 coefficients) | `weight4_M3/w4_symbol.py` | `python w4_symbol.py` (`--export-gpl` rewrites `gpl_expression.tsv`) | exact |
| Δ₃,₁(W₄ − 𝓡) = 0 and Δ₁,₃(W₄ − 𝓡) = 0 (Theorem 4) | `weight4_M3/w4_coproduct.py`, `w4_right_coproduct.py` | `python w4_coproduct.py`, `python w4_right_coproduct.py` | exact |
| The right weight-2 symbol is injective, so Δ₂,₂ = 0 by coassociativity | `weight4_M3/w4_kernel.py` | `python w4_kernel.py` | exact |
| Together with the motivic argument (Brown, Prop. 5.2; Borel) this gives W₄ − 𝓡 = (c₀ + √2 c₁)π⁴ with c₀, c₁ rational (Theorem 4) | (written argument, not code) | | c₀ = c₁ = 0 is still open |
| Fresh quadrature against the RHS, and the stored 1200-digit values | `weight4_M3/w4_numeric.py` | `python w4_numeric.py 160` | numerical: 160 digits, and 1.8e−1202 against the stored values |
| Certified bound: \|c₀\| < 10⁻³²², \|c₁\| < 10⁻³²¹ (ball arithmetic, every log called with the holomorphy flag) | `weight4_M3/w4_certified.py` | `python w4_certified.py 1100 1200` (a few seconds) | certified (python-flint/Arb balls) |
| **Rhombus end** | | | |
| P′(1) = 0 and P″(1) = (√2/3)π² (both known, see below); P‴(1) = 6√2 ζ(3) − 2√2 π²; Gamma-function operator formula for P | `rhombus_end/rhombus_formula.py` | `python rhombus_formula.py` | numerical: formula − direct ≈ 2e−21 at p = 1.01 and 0.99; third-order residual 1.27e−12 at d = 10⁻³ (≈ d⁴) |
| Semi-axes a, b (Theorem 5): P‴(1) = 24a²b²/r³ · (ζ(3) − π²/3) with r = √(a²+b²); the Taylor series at p = 1 has radius zero | `rhombus_end/r2_analysis.py` | `python r2_analysis.py` (writes `result.md`) | symbolic identities + finite differences at 4 aspect ratios; the zero-radius proof is a written argument that the script checks numerically |
| Exact coefficients a₀…a₁₀ of P/√2 = Σ aₙ (p−1)ⁿ (Proposition 2 states a₀…a₇); ζ(n) first appears in aₙ; γ cancels | `rhombus_end/rhombus_exact_coefficients.py` | `python rhombus_exact_coefficients.py` | exact (sympy); checked independently by `verify_all.py borel` to 89–90 digits |
| Growth law fₙ ~ √2 cos(π(3−2n)/4) n! (2/π)ⁿ n^(−3/2) and the Borel radius π/2 (Section XI-A) | `verify_all.py` | `python verify_all.py borel` | numerical: Richardson constant 1 ± 6e−6 |
| The divergent series is Borel summable (Section XI-A): a [90/90] Padé–Borel sum from f₀…f₁₈₀ (450 digits) reproduces direct quadrature of I(η) at η = 0.1, 0.2, 0.3, −0.2 | `rhombus_end/pade_borel_check.py` | `python pade_borel_check.py` (about 12 min; output in `pade_borel_check.out`) | numerical: 44, 31, 25, 30 digits |
| Radius exactly 1 at the square end (Theorem 9); decay M_k ≈ −(√π/k) e^(−√(πk)) cos √(πk) (Observation 1) | `verify_all.py` (data: `squircle_square_end/outputs/mk_deg9.tsv`) | `python verify_all.py mk` | the radius is a proof (no code); the decay law is numerical |
| **π_p (the squircle measured with its own ℓᵖ ruler)** | | | |
| π_p = 4(1 + N₀/p + N₁/p² + N₂/p³ + N₃/p⁴ + …), N₃ = log⁴2/24 + (π²/24)(1 − log 2) + ζ(3)/4 (Theorem 10) | `verify_all.py` | `python verify_all.py pip` | numerical check of the proved closed forms, 59–60 digits; N₄, N₅ and the divergent growth |
| **Ahmed's Problem 2: corners at any angle** | | | |
| Mellin transform of Ahmed's density h_c as a Legendre integral | `corners_ahmed_P2/v1_mellin.py`, `v5_mellin_barnes.py` | `python v1_mellin.py` | numerical, ≤ 1.6e−56 |
| Exact arc length L_θ(n) for every angle (Theorem 11; θ = π/2 gives Ahmed's Theorem 2) | `corners_ahmed_P2/v2_exact_length.py` | `python v2_exact_length.py` | numerical: against direct arc-length quadrature at 7 angles × 8 values of n, worst case 1.6e−42 |
| Coefficients a₁…a₄ in closed form (Theorem 12), and the all-order multiple-polylog series | `corners_ahmed_P2/v3_taylor_coeffs.py`, `v4_asymptotics.py` | `python v3_taylor_coeffs.py`, `python v4_asymptotics.py` | numerical, ≤ 5.4e−71; remainder ≤ 6.3e−67 at n = 10¹⁴ |
| Link to the squircle corner constant K(α) | `corners_ahmed_P2/v6_K_bridge.py` | `python v6_K_bridge.py` | numerical, at θ = π/3, π/2, 2π/3, 3π/4: moment and finite-M identities to ≤ 1.4e−47; the full residue sum Σ Λ′_m = 2K(π−θ) by extrapolation, to ≤ 8.5e−42 |
| Angular forms and the constants table; Ahmed's Table IV reproduced | `corners_ahmed_P2/v7_angular_and_table.py`, `v8_ahmed_table.py` | `python v7_angular_and_table.py`, `python v8_ahmed_table.py` | numerical: all of Ahmed's printed digits |
| Mellin–Barnes formula for any monotone corner curve | `corners_ahmed_P2/v9_general_curve.py` | `python v9_general_curve.py` | numerical, 2.3e−41 |
| **Smoothed corners of any angle and polygons** | | | |
| M₁(α) for every interior angle: ζ(2) + 27 real dilogarithms + 5 log products (Theorem 13) | `corners_p_smoothed/corner_M1_theorem.py` | `python corner_M1_theorem.py` | numerical check of the proved formula at 6 angles, ≤ 5e−38 |
| Reduced closed forms of M₁ at 3π/4 (octagon) and π/6, and at 2π/3 | `verify_all.py` | `python verify_all.py corner` | numerical: PSLQ/LLL forms, 68–69 digits here (500 digits when found) |
| Smoothed 2m-gons (Theorem 16): P(1) = 4m(sec(π/2m) − cos(π/2m)), P′(1) = P(1)·R₁; hexagon P″(1), P‴(1) with −(9√3/8)ζ(3) | `polygons_p1/hexagon_p1_check.py` | `python hexagon_p1_check.py` | P′(1) law proved; P″, P‴ numerical (PSLQ at 170 digits), finite differences here to 19/17/13 digits |
| **Ahmed's Problem 3: the length of the graph of fⁿ** | | | |
| 2 − L_n = (ln n + ln f′(1) + 1 − ln 2)/(f′(1)(n−1)) + O(ln²n/n²), and the expansion to all orders (Theorem 18) | `graphs_fn_ahmed_P3/run_defects.py`, then `analyze.py` | `python run_defects.py` (about 2 min on 12 cores; `data_defects.json` is already included), then `python analyze.py` | numerical: D_n to ≥ 60 digits for 6 functions, n = 10³…10²⁰ |
| Layer-cake representation, Mellin transform, moments | `graphs_fn_ahmed_P3/verify_lemmas.py` | `python verify_lemmas.py` | numerical at 50 digits: worst relative residual 1.2e−40 (one complex Mellin point), the rest ≤ 1e−50 |
| Q₁…Q₄ in closed form; Ahmed's Prop. A reproduced; the slope-2/e principle | `graphs_fn_ahmed_P3/symbolic_Qk.py`, `compact_forms.py`, `slope_principle_sym.py` | `python symbolic_Qk.py` (about 3 min) | exact (sympy) |
| Contact of order m ≠ 1 | `graphs_fn_ahmed_P3/degenerate.py`, `contact_all_m.py` | `python degenerate.py` | numerical, up to n = 10⁴⁰ |
| Corners at any angle, with general f; direct cross-check; locality | `graphs_fn_ahmed_P3/angle.py`, `crosscheck.py` | `python angle.py`, `python crosscheck.py` | numerical (crosscheck: 9.4e−82) |
| Growth of Q_k(0) (a remark on divergence) and the figure | `graphs_fn_ahmed_P3/growth_Qk.py`, `make_figure.py` | `python growth_Qk.py` | numerical observation |
| **Surfaces** | | | |
| Ahmed's Problem 4 (Theorem 19): the area of z = xⁿyⁿ is A(n) = 1 + √2/n + O(log²n/n²) | `surfaces/area_xnyn.py` | `python area_xnyn.py` (about 3 min) | numerical, 25 digits, n up to 10⁵ (this script covers the leading term; the terms through n⁻³ are in `new_results/P4c_area_xnyn/`) |
| Elementary bounds for R(n) = n²(A(n) − 1 − √2/n), giving R(n) = −½log²n + O(log n log log n) (Section XIV-A) | `verify_all.py` | `python verify_all.py area` | bounds checked against an independent quadrature for n = 10…10⁶ |
| Superellipsoid (Theorem 20): an exact one-dimensional formula for S(p), also in n dimensions | `surfaces/superellipsoid/07_area_checks.py`, `20_ndim_check.py` | `python 07_area_checks.py` | numerical: S(1) = 4√3 and S(2) = 4π to 1e−35; direct surface integrals to 20 digits |
| No log p terms; s₁ = 48 M₀; s₁, s₂, s₃ read off directly from S(p) | `surfaces/superellipsoid/09_expansion_check.py`, `16_richardson_check.py` | `python 16_richardson_check.py` | numerical: 60, 50 and 40 digits |
| Closed form of the vertex constant Q₀ (dilogarithms over ℚ(√2)) | `surfaces/superellipsoid/13_proof_Q0.py` (uses `dilog_engine.py`, `qs2.py`); found first with `04`, `06`, `10`, `11` (PSLQ at 300 digits); Q₀ computed two ways with `01`, `02` | `python 13_proof_Q0.py` | exact; every step is also checked to 60 digits |
| s₂ in closed form | `surfaces/superellipsoid/19_final_numbers.py` | `python 19_final_numbers.py` | numerical, residual 6e−410 |
| Q₁ (for s₃) by three independent routes; three of its five pieces identified | `surfaces/superellipsoid/08_Q1.py`, `14_Q1_pieces.py`, `15_Q1_highprec.py`, `17_Q1_subpieces.py`, `18b_lll_Q1_squircle.py` | `python 17_Q1_subpieces.py 400`, `python 18b_lll_Q1_squircle.py` | numerical, 400 digits (LLL identification) |
| The piece Z_b of Q₁ is 2M₁ − 2M₀ + π²/12 − (log²2)/2 | `surfaces/superellipsoid/25_proof_Zb.py` (uses `qs2.py`, `dilog_engine.py`); its numerical inputs checked by `26_check_2d_theta.py` | `python 25_proof_Zb.py`, `python 26_check_2d_theta.py` | exact (difference is exactly 0); the steps also agree to 40 digits |
| The value of Q₂ (for s₄); the second coefficient in n dimensions | `surfaces/superellipsoid/12_Q2.py`, `21_ndim_s2_check.py` | `python 12_Q2.py 45` | numerical: 43 digits; 10 digits for n = 4, 5 |

Prior work: P′(1) = 0 is in Montilla, Rodríguez and Linares (arXiv:2607.04592, Cor. 3.5). P″(1) = (√2/3)π² is in Rodríguez and Montilla (arXiv:2607.09048, Cor. 5.3). Problems 1 to 4 are Ahmed's (Zenodo, 2026).

**Exploratory, no claim depends on them:** `surfaces/superellipsoid/18_pslq_Q1_squircle.py` (the PSLQ run that came before `18b`), `22_lll_vertex_try.py` and `23_V13_hp.py` (an attempt to identify the remaining pieces of Q₁; `python 22_lll_vertex_try.py 1250` reads `V1_1250.json` and `V3p_1250.json` and finds no relation, see `out_22_lll_1250.txt`), `27_lll_Za_M2.py` (an LLL test of whether Z_a and V₄ are M₂ plus simpler constants; no relation turned up), and `05_sunit_search.py` and `24_sunit_zeta8.py` (S-unit searches for basis candidates, over ℚ(√2) and ℚ(ζ₈); `24` writes `sunit_zeta8.json`).
`corners_ahmed_P2/v6_funcs.py` is a helper (Λ′_m through the degree derivative of P_ν, and K(α)) for an extrapolation experiment that is not in this repository.
`squircle_square_end/scripts/check_tex_formulas.py` re-checks every displayed formula in the note's TeX source. That source is not in this repository, so the script only runs with `SQUIRCLE_TEX=/path/to/main.tex`. When `SQUIRCLE_TEX` is set, `weight4_M3/w4_symbol.py` uses the same file to check that the coded right-hand side matches the typeset formula.

## New results (`new_results/`)

The folders in `new_results/` hold the later computations behind Sections X to XIV of the paper, one folder per computation, named after the job that produced it. Each has the key scripts, the small outputs (`*.out`, `*.json`), `result.md` (the write-up of the job), `verify*.md` or `verify/` where a second reader re-checked it, and `planner_check*.py` (a separate re-computation).
Transcripts and pickles are not included, nor are short scratch steps that no write-up cites, and absolute paths in the copied text were replaced by `research/…`. `new_results/README.md` has the details.
Numbers are those of the paper. *Proved†* means that the full written proof is in the research notes, not in this repository.

| Paper item | Status | Folder |
|---|---|---|
| Theorem 6: pure weight in τ = 1 − 1/p (rhombus; the surface of Σ\|xᵢ\|^p = 1 in ℝⁿ); the 2D series is exact through τ¹² | proved† | `D16_pure_weight/` |
| Proposition 2: ω_n, the coefficient of ζ(n) in aₙ (ω₃ = 1, ω₅ = −3/4), and its rhombus version | proved; ω_n ≠ 0 checked for odd n ≤ 21 | `D11c_zeta_first_order/` |
| Theorem 7: the Fourier coefficients of the support function h_p at the diamond end, every k, through τ⁴ | proved† | `F1_fourier_squircle/` |
| Proposition 3: the degree ⌊n/2⌋ − 1 of Θₙ(k²); the closed forms Θ₅–Θ₇ and the τ⁵ term | degree proved†; closed forms and τ⁵ numerical (PSLQ) | `K3_fourier_all_orders/` |
| Observation 1: M_k ≈ −(√π/k) e^(−√(πk)) cos √(πk), with the first correction in closed form | numerical (independent Cauchy integrals up to k = 2400) | `D13b_Mk_constant/` |
| Section XI-D: the p = ∞ coefficients M_k from the Borel data at p = 1 (Stokes jump ∝ ε^(1/2)) | numerical, mechanism derived | `D13c_two_ends/` |
| Proposition 4: L₁(p) = 8 and L_∞(p) = 8·2^(−1/p); Section XI-E: K₀(q) at q = 3, 4 | proposition proved†; K₀(3), K₀(4) numerical (no classical relation at 110 digits); the ζ(3) coefficient as a function of q is a conjecture | `D12c_ruler_family/` |
| Theorem 13: M₁(α) for every angle | proved | `D10c_corner_general/` |
| Theorems 14, 15 and Proposition 5: the flat-angle expansion, M₂(α) = W₂ + W₃; M₂(α) with Li₃ at eight arguments | proved†; the Li₃ closed form numerical (300 digits, 9 angles) | `G1_zeta3_any_angle/` |
| Theorem 17: C_m, the mixed-term ζ(3) coefficient of P‴(1) for the smoothed regular 2m-gon | proved† (see the correction below) | `D15c_edge_zeta3/` |
| Proposition 6: the local terms Z_e of an edge; the parallelogram and affine-regular hexagon laws; a hexagon with more ζ(3) than ΣZ_e | proved†; the hexagon counterexample numerical (residual 8e−141) | `E1_edge_law/` |
| Section XII-C: the corner–corner ζ(5) law C_m⁽⁵⁾ | numerical (exact at m = 2: −90√2) | `D17_zeta5_polygons/` |
| Theorem 19: Ahmed's Open Problem 4, A(n) through n⁻³ | proved† | `P4c_area_xnyn/` |
| Theorem 21: the octahedron end, S(1+d) to d¹², single zeta values in every dimension | proved†; the Borel singularities numerical; zero radius for n ≥ 3 a conjecture | `H1_octahedron_end/` |

**Read this before the polygon folders.** Some write-ups in `D15c_edge_zeta3/` and `E1_edge_law/` say that all of the ζ(3) in P‴(1) comes from the corner–corner terms. That is true, and proved, only for *fully local* polygons (every parallelogram, every affine-regular hexagon), so C₂ and C₃ are the whole ζ(3) coefficient for the diamond and the regular hexagon.
For a general polygon ΣZ_e is only the corner–corner part. For K_p = {\|x\|^p + \|y\|^p + 3\|x−y\|^p ≤ 4^p} the ζ(3) coefficient of P‴(1) is (231/125)√5 − (15/64)√2, while ΣZ_e = (72/125)√5 − (15/64)√2.
For regular 2m-gons with m ≥ 4, C_m is the local part, and whether it is the whole ζ(3) coefficient is open. Its size is C_m = −3π⁴/m³ + O(m⁻⁵), not 9π⁵/(16m⁴) as `D15c_edge_zeta3/result.md` first said.

## How verification works

**Exact checks (proofs).** Every algebraic number in the square-end computations lies in ℚ(√2). The code stores it exactly as p + q√2 with Python `Fraction`s (`qsqrt2.py`; `qs2.py` for the superellipsoid).
Every logarithm that appears is the logarithm of an S-unit ±2^(j/2)(1+√2)^m. It is therefore stored exactly as a vector over log 2 and L = log(1+√2).
Integrals are evaluated with exact antiderivatives. Each antiderivative is checked by symbolic differentiation in `sympy`.
Polylogarithm values are reduced only with classical functional equations:
- for Li₂: reflection, inversion, Landen, duplication and Abel's five-term relation;
- for Li₃: inversion, the three-term relation, duplication and the Spence–Kummer equation on its real domain.

A result counts as proved when the final difference is exactly zero in this formal ring.
Floating point appears only in separate sanity checks. The weight-4 scripts go further: they make any conversion of a ℚ(√2) number to a float raise an error (`ban_float`).

**Numerical checks (evidence, not proof).** Values are computed with `mpmath` at the stated working precision (`mp.dps`), using tanh–sinh quadrature with breakpoints at the singular points.
Each key value is computed by two independent methods, or at two precisions, and the table above reports how far they agree.
Integer relations are found with PSLQ (`mpmath.pslq`) or LLL (`python-flint`, `fmpz_mat.lll`). A relation is accepted only if:
- it comes back unchanged at several precisions and from independent data;
- its coefficients are far smaller than a random relation would need at that precision.

`numpy` and `matplotlib` are needed only for the figure scripts.

## Lean proofs (`lean/`)

`lean/` is a Lean 4 project (Mathlib) that machine-checks 580 statements taken from the paper (577 distinct), in 71 files under `lean/Zeta3/`. Each one is proved from the axioms of Lean with no `sorry`.

- `Wave1/` (66 statements): exact finite claims. Examples: the two expansions at the diamond end agree through order five; C₂ = 6√2, C₃ = −9√3/8 and C_m < 0 for m ≥ 3; the edge sums over the regular 2m-gon equal C_m, C_m⁽⁵⁾ and C_m⁽²,³⁾; L₁(p) = 8 and L_∞(p) = 8·2^(−1/p).
- `Wave2/` (57 statements): a dilogarithm library; Lemma `lem:dilog`; the closed form of M₁ from its integral; M₁(α) at every angle (the table form and the four-dilogarithm form); ω_n ≠ 0 for every odd n ≥ 3; R₁ > 0 for m ≥ 3.
- `Wave3/` (30 statements): the integral for P(p) is the true length of |x|^p + |y|^p = 1; P(p)/8 has a convergent expansion in 1/p for every p > 1; its first two coefficients are M₀ = log 2 − √2 log(1+√2) and the closed form of M₁ (Theorem `thm:M1`), with the printed decimals to 10⁻²⁰; the length of the graph of fⁿ (Theorem `thm:fn`).
- `Wave4/` (74 statements): for every rhombus, P is three times differentiable at p = 1 with P′(1) = 0, P″(1) = 4π²A²B²/(3r³) and P‴(1) = (24A²B²/r³)(ζ(3) − π²/3) (Theorem `thm:rhombus`); the diamond coefficients a₀…a₅ with ζ(3) and ζ(5); the expansion of π_p through N₃; convergence of the ℓ^q rulers; P′(1) = P(1)R₁ for the smoothed regular 2m-gon; the area between xⁿ and yⁿ through order n⁻³ with the explicit bound 45(1+log n)⁴/n⁴ (Theorem `thm:xnyn`).
- `Wave5/` (63 statements): the closed form of M₂ with ζ(3), Li₃(±b), Li₃(1/√2), and its printed decimal to 10⁻²⁰ (Theorem `thm:M2`); the exact length of the graph of xⁿ at every angle and its expansion to O(log n/n³) (Theorems `thm:P2`, `thm:P2asym`); the graph of fⁿ at every angle (Theorem `thm:fnangle`); the Fourier coefficients of the support function through τ⁴ (Theorem `thm:fourier`); the octahedron end of the surface area through τ⁶ (Theorem `thm:oct`, for the face integral); the Taylor series of every rhombus perimeter at p = 1 has radius zero, and P is not analytic there (Theorem `thm:zero`).
- `Wave6/` (78 statements): M₃ = W₄ + W₃′ + T from the curve (Theorem `thm:M3`); the series ΣM_kεᵏ has radius exactly one and ε = 1 is its only singular point on the unit circle (Theorem `thm:one`); the Fourier coefficients through τ⁵ with Θ₅; the four special values of Proposition `prop:nonlocal` and the Nielsen–Kölbig identity; the area formula of the superellipsoid with its first two terms at the cube end (Theorem `thm:S`); the next two orders of Theorem `thm:P2asym` with explicit constants; M₂(α) at every angle (Theorem `thm:M2alpha`); the printed decimals of M₀, κ and M₂ are truncations.
- `Wave7/` (36 statements): Cauchy's formula P(p) = π·a₀(h_p) for the true perimeter, and with it the series of Theorem `thm:pw` through τ⁵; Θ₆; S(1) = 4√3 and S(2) = 4π; Theorem `thm:flat` for k = 0, 1, 2, 3, with the coefficient −ζ(3)/(4n²−1) in M₃(α).
- `Wave8/` (66 statements): Theorem `thm:pw` to every order, for the squircle and for every rhombus: log(P/(4√2)) = Σ w_k τᵏ with each w_k a rational polynomial in ζ(2),…,ζ(k) of weight k; Theorem `thm:flat` for every k, and the first rung of the zeta ladder; Θ₇; Theorem `thm:pip` in full: the expansion of π_p in 1/p to every order, the closed forms of N₄ and N₅, the growth law of N_k for odd and even k, and the divergence of the series for every p; the decay law, Theorem `thm:decay`, in full: |M_k + (√π/N)e^(−√(πN))cos√(πN)| ≤ 1.21·10⁹ N^(−3/2) e^(−√(πN)) for N = k+1 ≥ 100, with every lemma of its proof.
- `Wave9/` (56 statements): the first correction to the decay law, Theorem `thm:varpi`, in full, with constant 1,000,668; the zeta ladder, Theorem `thm:zladder`, in full (every zeta monomial, with the rungs for ζ(4), ζ(5), ζ(6), ζ(2)ζ(3), ζ(3)²); Θ₈; Proposition `prop:degree` for every n (degree exactly ⌊n/2⌋−1, the leading coefficients, the atoms as iterated integrals); the last sentence of Theorem `thm:pw` (pure weight for the surface area in every dimension) and both weight claims of Theorem `thm:oct`; the expansion of the superellipsoid area to every order with no log p terms, and s₂, s₃ (Theorem `thm:S`).
- `Wave10/` (54 statements): Θ₉, which completes Proposition `prop:theta` through τ⁹; Theorem `thm:M4alpha` with its ζ(3), ζ(4) and ζ(2)² parts in closed form (the parts R₄ and Ψ₂ as convergent series); Theorem `thm:hexagon`: P‴(1) as a sum over edges, the eleven-term formula for every hexagon and every ratio r > 0 (the five printed coefficient blocks as printed, the other six as integer polynomials that exist), the ζ(3) block, the link to the non-local values, the hexagon law at r = 1; the field rule, Theorem `conj:field`, in the form of its last sentence, for every k: M_k(α) is a finite combination of harmonic polylogarithms H_w(c) with letters 0, 1, −1, with zeta monomials and rational functions of c as coefficients and total weight at most k+1; the closed form of the vertex constant Q₀ (Theorem `thm:S`), with s₂ in closed form and its decimals −4.6892155184…

Lean does not check everything. The identity W₄ = 𝓡 (Theorem `thm:W4R`, a computer certificate), parts (a) and (b) of the field rule, the weight-two theorem `thm:wt2`, Theorem `thm:hexagon5`, the closed forms of R₄ and Ψ₂ and the results about symbols are not in Lean; they rest on the written proofs and the scripts in this repository. `lean/LEAN_CHECKLIST.md` lists, target by target, what is proved and what is not covered.
Each file states its claims between `BEGIN STATEMENTS` and `END STATEMENTS`; `lean/audit/` holds the statement blocks as accepted by a separate audit against the paper, and `lean/tools/gate.py` checks that a file still has those statements, compiles, has no `sorry`, and uses only the standard axioms.
To rebuild: install `elan`, then in `lean/` run `lake exe cache get` and `lake build`.

## Files in each part

```
squircle_square_end/  scripts/ (M_k at p -> infinity: exact Q(sqrt2) polylog engine, proofs, numerics)  outputs/ (stored values and logs of the runs)
rhombus_end/          rhombus_formula.py, r2_analysis.py, rhombus_exact_coefficients.py, pade_borel_check.py (+ .out)
weight4_M3/           w4_symbol.py, w4_coproduct.py, w4_right_coproduct.py, w4_kernel.py, w4_numeric.py, w4_certified.py, gpl_expression.tsv
corners_ahmed_P2/     common.py, closed_forms.py, v6_funcs.py, v1..v9 scripts with reference *.out
graphs_fn_ahmed_P3/   fnlib.py, symbolic_Qk.py, analyze.py, ... with reference *.out and data_defects.json
surfaces/             area_xnyn.py; superellipsoid/ (01..27 scripts, libraries, stored values)
corners_p_smoothed/   corner_M1_theorem.py
polygons_p1/          hexagon_p1_check.py
verify_all.py         one-command independent re-check of the main results
new_results/          one folder per later computation (README.md there lists them)
extra_checks.py       checks added in review round 1 (Faulhaber, coefficient of ζ(n) in aₙ, Ahmed's Table IV, corner remainder, oblique corner)
```

## Verification map

Every row of the paper's verification table (Table II) is mapped to its method, re-check and scripts in [VERIFICATION_MAP.md](VERIFICATION_MAP.md).
