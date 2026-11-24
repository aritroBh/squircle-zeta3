# Referee report: prop:nonlocal and thm:hexagon (main.tex lines 1630-1660)

All code here is the referee's own. Nothing from the project is imported. The project's data are only read: the text of `formula_table.txt`, and, for comparison, `E3_fixed.pkl`, `Jtab_t.pkl` and `jloc.json`.
Python: `.venv/bin/python` (mpmath, sympy).

## Verdicts

- **prop:nonlocal: PROVED (computer-assisted)**, with one wording repair (R4).
- **thm:hexagon: PROVED AFTER REPAIRS (computer-assisted)**. Repairs R1 to R6 are below. None of them changes a number.

## 1. Statement checks (numerical, independent)

**Direct method (`direct.py`).** The polar perimeter integral is
P(p) = ∫ F^{-1/p} sqrt(1+(F_θ/(pF))²) dθ, with F = Σ_j w_j^{1-p} |c_j cross(V_j,u)|^p.
This is exactly the weighted family w_j|k_j⟨v_j,x⟩|^p with w_j k_j = c_j.
- It differentiates under the integral sign. At each node the integrand is a truncated power series in d = p-1, to order 3.
- The d³ coefficient is integrated by tanh-sinh, split at every kink direction. Each sector is integrated from both of its corners, so there is no cancellation near a corner.
- Finite differences in p (5-point stencil with Richardson extrapolation, `fd.py`) check the interchange of derivative and integral.

**The closed form compared against.** main.tex does not print the coefficients of E_e. It prints only the structure, the ζ(3) coefficient Z_e, and the two Li₃ coefficients. So the full closed form was transcribed from the text of `formula_table.txt` (`formula.py`), with "ln(1-x)" read as ℓ = -log r, as main.tex defines ℓ. The edge data (B, r, q0, L², |V0|², A0, A1, the weights) are computed from scratch from the corners.

| hexagon (corners V0,V1,V2; weights) | ratios r | P'''(1) direct (dps 100) | table − direct | digits agreeing |
|---|---|---|---|---|
| H1 regular, w=1,1,1 | 1,1,1 | -5.14545501713561583717168113112397507111555 | 1.3e-98 | ~98 (also equals the printed regular-hexagon closed form) |
| H2 paper example (1,0),(2,2),(0,1), w=1,3,1 | 2, 1/2, 1 | 19.9002663631480669804245454749161564779543 | 6.2e-98 | ~98 |
| H3 far from regular: frame a=2, b=3/2, map [[3,1],[0,1/2]], w=2,7,3 | 2, 2/3, 3/4 | -46.1611422454261048715611608268283107936406 | 2.8e-97 | ~97 |
| H4 frame a=13/10, b=2, map [[1,-2],[1,1]], w=5,1,2 | 13/10, 1/2, 20/13 | -81.8681370772164540567168691618347239632655 | 1.1e-97 | ~97 |
| H5 irrational: a=√2, b=√3-1/2, map [[1,√5],[0,1]], w=3/2,4,1 | √2, ..., ... | 34.5009456800334132134845645691948928087294 | 1.4e-97 | ~97 |
| H6 a=4, b=16/5, w=1,2,3 | 4, 5/16, 4/5 | -7.31516326661736979166530279796929613429888 | 2.3e-97 | ~97 |
| H7 a=2, b=3/2, map [[2,1],[1,3]], w=1,5,3 | 2, 2/3, 3/4 | -11.6758700154996667469020291021529154694693 | 2.1e-97 | ~97 |

"Frame a,b" means the corners (1,0),(a,b),(0,1), then mapped by the 2×2 matrix shown.

- Each edge also agrees separately, to the same precision (`log_tests_dps60.txt`, `log_tests_dps100.txt`).
- At dps 60 every residual is about 1e-57. At dps 100 it is about 1e-97. So the residuals are working precision, not a real discrepancy.
- Finite differences, which do not differentiate under the integral: H3 agrees with the table to 5.5e-21 and H6 (r=4) to 2.2e-20 (`log_fd_*.txt`).
- The paper's own three test hexagons (`paper_cases.py`) were recomputed. They reproduce the printed ratios and agree with the table to about 1e-57:
  - (1,0),(3,3),(0,1), w=3,5,3: -44.57991977168253341002772110472172...
  - (2,0),(3,2),(0,1), w=2,5,7: 11.17600084483001768696522886994283...
  - golden hexagon, w=2,3,5: -30.85385931932535892431770994199399...
- **Read literally** (ln(1-x) = log r, as the labels in `formula_table.txt` say), the table is wrong by O(1) to O(60) on every edge with r ≠ 1. See R1.

**Test of prop:nonlocal, independent of the engine (`pslq_nonlocal.py`, dps 120).** For each edge, X = E_e·L comes from the direct jet integral. An integer relation is found on [X, ζ3, π²log2, log³2, π², log²2, log2, 1]. Its ζ(3) coefficient minus the printed Z_e is compared with the printed N_e:

| edge | r | ζ(3) coefficient (integer relation) | Z_e (printed) | N_e (printed) | difference | relation residual |
|---|---|---|---|---|---|---|
| H2 edge 0 | 2 | 231√5/500 | 18√5/125 | 159√5/500 | 0 | 8e-116 |
| H2 edge 1 | 1/2 | 231√5/500 | 18√5/125 | 159√5/500 | 0 | 8e-116 |
| H2 edge 2 | 1 | -15√2/128 | -15√2/128 | 0 | 0 | 4e-117 |
| H7 edge 0 (fresh) | 2 | 6507√170/157216 | 945√170/314432 | 12069√170/314432 | 0 | 3e-112 |
| H8 edge 1 (fresh: a=3/2, b=2, map [[1,1],[-1,2]], w=3,5,1) | 1/2 | 17253√26/21632 | 30375√26/281216 | 96957√26/140608 | 0 | 2e-113 |

The relations use coefficients up to about 2.4e7 on 8 numbers. That uses about 60 of the 113 available digits.

## 2. Derivation check: how each step is established

Pipeline in the notes: per-edge d³ integrand E3 (K1/K1d) → split by powers of W = log(1+κt) → local integrals (Beta function) and the J-table (`esum_engine`) → assembly with A eliminated (K1e step_23) → 11-term table.

| step | how it was established in the notes | what the referee did |
|---|---|---|
| Setup: arc of ∂K_p by radial projection, x(t) = y(t)·F_p(y(t))^{-1/p}; P''' = 2Σ over three edges by central symmetry | by hand; correct | The polar sectors equal the radial-projection arcs, because the kink lines pass through the origin. Confirmed by the per-edge agreement. |
| E3 is a polynomial in t, log t, log(1-t), log(1+κt) | exact computer algebra (K1 step_33; a 1/8 vs 1/4 bug was fixed in K1d) | Re-derived independently in sympy, using F' = (1+d)Σ a_j' e^{d u_j} and F = 1 on the edge. It matches `E3_fixed.pkl` exactly at 4 random rational parameter points (`integrand_check_fast.py`). The fully symbolic comparison was too slow and was stopped. |
| Local/non-local split | exact: a partition of the monomials by the power of W | trivial; no issue |
| Local integrals ∫tⁿ logᵃt logᵇ(1-t) | exact Beta-function derivatives | all 40 entries of `jloc.json` match quadrature to 8e-45 |
| J-table (c ≥ 1) | exact engine algebra: power series of logᶜ(1+κt) via Stirling/harmonic sums, partial fractions; Euler sums reduced to Li₁, Li₂, Li₃ and A | Not re-derived symbolically. The engine logic was read and the formulas it uses (M(m;a,b), the series of log², log³) are standard. All 40 entries match the referee's quadrature at κ = -0.9, -0.5, 0.4, 1, 1.7, 3.5, 9 (worst relative error 1.5e-42). The atoms Σ H⁽²⁾ₙxⁿ/n and Σ Hₙ²xⁿ/n match at 3 points. |
| Validity of the J-table for r ≥ 2 (κ ≥ 1) | **not argued.** The κ-series only converges for |κ| ≤ 1. | Repair R2: one sentence on analytic continuation (identity theorem). Numerically confirmed at κ = 1.7, 3.5, 9 and on hexagons with r = 3, 4. |
| Eliminating A with the Nielsen–Kölbig identity | classical; the notes checked it numerically; CX2 proved it by differentiation | Re-proved: the sympy derivative difference is 0 (with Li₂ Landen), and both sides vanish at 0. A = Li₃ + S₁,₂ continuation checked at x = -2, -5, 0.6. |
| Assembly into 11 terms | exact sympy (K1e step_23) | Not re-run symbolically. Checked numerically as an identity in free parameters: table vs the referee's own integrand integrated in t at 5 random parameter sets with r up to 5.65 (`random_params.py`, run stopped after 5 of 10 points for time), residuals about 1e-43. |
| Coefficient of ζ(3) = Z_e; coefficient of Li₃(x) = 2Φ_A+Φ_L; coefficient of Li₃(x/(x-1)) = Φ_A | exact | Re-done exactly from the table text against the **printed** Z_e, Φ_A, Φ_L: all three differences are 0 (`symbolic_checks.py`). |
| Special values A(-1), A(1/2), Li₃(-1), Li₃(1/2) | classical; numerical in the notes | Re-checked to 50 digits. They are classical results. |
| Differentiating under the integral in p | **not argued** | Repair R3. Finite differences agree to 2e-20. |
| Orientation independence of E_e | implicit | Table checked: E(forward) = E(reversed) to 1e-39. |

**No coefficient of the theorem was found by fitting.** Integer relations appear only in the checks. In K1b the ζ(3) coefficients of the paper's example were first found that way; K1d later derived them from the engine, and here they are confirmed independently on fresh edges.

**Residual trust.** Two steps were not re-derived symbolically by the referee: the J-table engine algebra and the final symbolic assembly. Both are exact computer algebra in the notes, and both are confirmed to 40 or more digits over a wide range of parameters.

## 3. Hypotheses

- The theorem says "every centrally symmetric convex hexagon". The derivation needs a non-degenerate (strictly convex) hexagon, so that 0 < B < 1, r > 0 and A0 = 1 - Br > 0. A "hexagon" with a straight angle is a parallelogram and is not covered (log A0 diverges). Within the strictly convex class, the formula covers every r > 0, provided the continuation of R2 is stated. Tested at r from 5/16 to 4, irrational r, and random r up to 5.65.
- The result is for the weighted family. P'''(1) depends on the weights w_j: the same affine-regular shape gives 4.3254... with w = 1,1,1 and -44.32... with w = 2,3,5. The paper says this in the paragraph before the theorem. The theorem statement itself should say "for the weighted family above".

## 4. Repairs and gaps

- **R1 (notes, required).** `formula_table.txt` (and the "paper paragraph" in K1e result.md) labels the atom ln(1-x). It is actually ℓ = -log r. Read literally, four blocks have the wrong sign. Only a header note says this. Since the paper points readers to "the table of coefficients in our research notes", relabel the blocks as ℓ (or negate those four blocks). Also state there that the denominator is L² to the power 5/2 (L⁵), and that lA_j = log(A_j/w_j) are constants. main.tex itself is consistent: it defines ℓ = -log r.
- **R2 (proof sketch).** Add: "The table is derived for |1-r| < 1 by expanding log^c(1+κt) in powers of κ. Both sides are analytic in r > 0 (for x = 1-r < 1), so they agree for every r > 0." This is needed for r ≥ 2, including the r = 2 case of prop:nonlocal (an Abel boundary case) and the paper's test with r = 3.
- **R3 (proof sketch).** Add one sentence justifying P'''(1) = ∫ ∂³_p(integrand). For |p-1| < 1/4 the p-derivatives are dominated by C·(1 + |log t|³ + |log(1-t)|³)·t^{-1/4}(1-t)^{-1/4}, which is integrable.
- **R4 (wording; both results).** "The ζ(3) of P'''(1)" and "In this basis the coefficient of ζ(3)" are unique only if the constants (respectively the 11 functions) are linearly independent, which is open (the paper says so elsewhere). Suggested wording:
  - prop:nonlocal: "Written in terms of ζ(3), π² times logarithms and products of logarithms, the coefficient of ζ(3) beyond Z_e is N_e = ..."
  - thm:hexagon: "In this presentation the coefficient of ζ(3) is exactly Z_e."
  - The 11 terms are not proved to be a basis. Use "combination" or "presentation" instead of "basis".
- **R5 (reduction at r = 1).** "The theorem reduces ... to the hexagon law of Proposition edgelaw at r=1" is true only when all three ratios are 1 **and all w_j = 1** (the unweighted convention of prop:edgelaw). With w = 2,3,5 on an affine-regular hexagon the law is off by 48.6.
- **R6 (hypothesis).** Say "strictly convex", or "non-degenerate", centrally symmetric hexagon.

## 5. Sentences in main.tex that claim more than is proved

1. Line 1650: "In this basis the coefficient of ζ(3) is exactly Z_e, for every r". The 11 terms are not shown to be a basis (R4).
2. Line 1650: "The theorem reduces ... to the hexagon law of Proposition edgelaw at r=1". This needs all w_j = 1 and r = 1 on every edge (R5).
3. Line 1645: "For every centrally symmetric convex hexagon". This should say strictly convex, and refer to the weighted family (R6, §3).
4. Line 1631: "the ζ(3) of P'''(1) beyond Z_e is N_e". This is unique only modulo the open independence of the constants (R4).
5. Line 1658: "The table of coefficients is in our research notes." The table there is mislabeled (R1).

## Files

- `direct.py`: polar jets and perimeter
- `formula.py`: table transcription and printed formulas
- `run_tests.py`, with `log_tests_dps60.txt` and `log_tests_dps100.txt`
- `fd.py`, with `log_fd_H3.txt` and `log_fd_H6.txt`
- `pslq_nonlocal.py`, with `log_pslq_nonlocal.txt`
- `symbolic_checks.py`, with `log_symbolic.txt`
- `integrand_check_fast.py`, with `log_integrand_fast.txt`
- `jtable_check.py`, with `log_jtable.txt`
- `random_params.py`, with `log_random_params.txt`
- `paper_cases.py`, with `log_paper_cases.txt`
- `log_S12.txt`
