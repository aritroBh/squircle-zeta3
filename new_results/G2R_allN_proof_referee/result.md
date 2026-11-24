# G2R referee verdict: the all-n proof of Theorem B4 is **PROVED** (no repairs needed)

Refereed object: `G2_M4alpha/step_11.py` (+ `step_10.py`, `g2_fit.py`, the five term lists, `mu4_formula.pkl`), following the Z-sum method of `b3_allN_proof.py`/`zsum.py`. All work in `new_results/G2R_allN_proof_referee` (`result.md`, `progress.log`, scripts + `.out` files below). Every number quoted was printed by a script I ran.

## 1. What step_11.py proves (my own words)

For each zeta sector P ∈ {1, ζ(2), ζ(3), ζ(4), ζ(2)²} separately, it proves for **all n ≥ 1**: `[c^{2n}] B4|_P = μ₄,ₙ|_P` and `[c^{2n+1}] B4|_P = 0`, plus the leftover coefficients c^{−3}…c³ checked directly (targets 0 except [c²] = μ₄,₁|_P). Coverage is complete: even N ≥ 2 / odd N ≥ 3 symbolically, N ≤ 3 directly, and no term has valuation < −3. With Theorem A (M₄ = Σ μ₄,ₙ c^{2n}, |c| < 1) and polynomial boundedness of μ, equal coefficients of analytic functions ⇒ B4 = M₄ on |c| < 1. The per-part split is stronger than needed (no independence assumption on the zeta monomials).

- **Z-sums/upper limits**: μ side uses H⁽ʳ⁾_{n−1}, H⁽ʳ⁾_{2n−1} (r ≤ 4); HPL coefficients give Z_v(N−1) at N = 2n+δ, δ = parity−e ∈ {−1,…,4}; everything canonicalised to Z_w(2n) (observed indices {−4..−1, 2..4}, depth ≤ 4; 29/5/2/1/1 words per part; odd forms cancel to empty).
- **Identities used**: quasi-shuffle (stuffle with a∘b = sgn(a)sgn(b)(|a|+|b|)); even/odd index splitting H⁽ʳ⁾_{n−1} = 2^{r−1}(Z_r(2n)+Z_{−r}(2n)) − n^{−r}; index shifts Z_{a,v}(2n+j) expanded recursively (j ∈ {−2,…,3}); the HPL Taylor recursion h_{a,w}(N) = coef·a^{N±1}/N·Z_{(aεk),v}(N−1). Partial fractions occur only upstream in step_07's F₁,F₂ derivation (re-verified: both π⁴ identities and the apart forms True).
- **Validity for every n ≥ 1, boundaries**: my denominator audit over all ten canonical forms found only factors {n, n+1, 2n−1, 2n+1, 2n+3} — no zeros at integers n ≥ 1; summation ranges stay ≥ 0. Boundary-verified exactly: 96 words × δ ∈ {−1,…,4} × n = 1,2,3 (1728 checks), shifts for 29 Z-words × j ∈ {−2,…,3} × n = 1,2 (348 checks), harmonic conversions — all True.
- **"Identical canonical forms ⇒ equality"**: yes — the check is `diff == {}` (the literally-zero ℚ(n)-combination) after exact rewrites only; relations among Z-sums are irrelevant in this direction, and the other direction is not needed.
- **c^e prefactors, e ∈ {−3,…,1}**: [c^N]c^eH_w = h_w(N−e) with correct ε^δ sign and He-partner sign (−1)^{e+nz(w)} (consistent with g2_fit and with the planner's independent reading). Cross-validated: my own evaluator == `hpl.coeffs_exact` (96 words, m < 84) == `hpl_coeff_at` at n = 1,2,3.
- **Low orders c^{−3}…c³**: checked exactly per part; independently re-checked (N = −3…1) against the definition.

## 2. The run, assertion points, non-vacuity

Pristine copy `run_base/`: rc = 0, **3.7 s**, final line `B4 PROVED for all n (all five parts): True`. Equality is asserted at `len(zsum.lin((1, ev), (-1, muE))) == 0`, `len(hpl_side(terms, 1)) == 0`, `low_order_check(...)` — with `ev` from the JSONs (B4 candidate) and `muE` from `mu4_formula.pkl` (step_05's independent derivation). **Not vacuous**: the odd identity of part '1' cancels 446 (word,coefficient) insertions to 0 surviving words; at n = 5 both sides equal the exact fraction −7596252093627541/219583673971200000 (= flat_exact's μ_{4,5} rational part — three-way agreement); and all four mutants below are caught. Regenerated `mu4_zsum.pkl` == stored artifact.

## 3. Mutation tests — all four FAILED as required (proof detects each)

| mutant | perturbation | even residual (word: ℚ(n) coef) | low-order mismatch | verdict |
|---|---|---|---|---|
| (a) R₄ weight-5 coeff | `c^0 He[0,0,0,0,1]`: 1→2 | Z_() : 1/(16n⁵) | c²: −1007/2160 vs −571/1080 | FAILED ✓ |
| (b) Ψ₂ coeff | `c^0 He[1]`: 5/4→19/12 | Z_() : 1/(3n) | c²: 247/432 vs 103/432 | FAILED ✓ |
| (c) F₁ coeff | artanh 7c/16→5c/16 | Z_() : −1/(16n−8) | c²: 17/120 vs 4/15 | FAILED ✓ |
| (d) prefactor exponent | `c^-2 He[1,0,1]`→`c^-3` (coef 2) | Z_() : (4n²+2n−1)/(8n⁵+36n⁴+62n³+51n²+20n+3); Z_(2) : −2/(2n²+5n+3) | c⁰: −1/3 vs 0; c²: −811/1080 vs −571/1080 | FAILED ✓ |

In each mutant only the affected part fails (per-part sensitivity); odd identities stay True (mutated objects remain even). No mutant passes ⇒ no bug found.

## 4. Independent spot check (task 4) — PASS

Own parser + own exact HPL evaluator vs μ₄,ₙ from the definition (G1's unmodified `flat_exact.py 40 4`, run by me, 92 s): **(4a) [c^{2n}](B4 part) == μ₄,ₙ part exactly in ℚ, per zeta monomial, all n = 1…40, all parts: True**; (4b) odd c¹…c⁷⁹ vanish: True; (4c) c^{−3}…c¹ vanish: True; (4d) `mu4_formula.pkl` (the proof's μ input) == definition for n ≤ 40: True.

Corroboration reproduced here: full `step_10.py` (147 s): **EXACT through c^1280 all parts True** + numerics 3.99e-92/2.55e-92/2.87e-92/1.75e-93/4.87e-93/1.56e-93 at six angles and |B4−M4_800| = **9.64e-121** at c=1/√2; planner's independent `planner_check_b4.py` (dps 130): **8.6732e-131** at the square, 4.8651e-93 (2π/3), 2.8715e-92 (π/3); `zsum.py` self-tests True.

## 5. Overall verdict: **PROVED**

Every method step: CORRECT (see verdict table in `result.md` §5). The only defect found is cosmetic: G2's `result.md` §2 states R₄'s weight census as "3/8/28/70" while the actual JSON census is **1/8/27/73** (total 109 either way) — documentation only, does not touch the proof. Scope notes (inputs, not gaps of step_11): Theorem A and step_05's symbolic μ formula — the latter cross-validated vs the independent `flat_exact.py` for n ≤ 640 (by step_05) and n ≤ 40 (by me).

## 6. Files that produced this

- **`verify_g2r.py`** (< 5 min; actual 7 s): pristine pipeline run (must say True) + mutant (a) (must say False) → `verify_g2r.out`: **PASS**, exit 0.
- `referee_checks.py`/`.out`: non-vacuity (446 insertions), boundary checks (1728+348), denominator audit, 96-word evaluator cross-check.
- `spot_check.py`/`.out` + `flat_run/flat_exact_N40_K4.json`: task-4 exact spot check.
- `run_base/` + `run_base.out` (step_11), `run_base_step10.out` (full step_10); `mut_a…d/` + `mut_a…d.out`; `progress.log` (one line per step); `result.md`. (The small `step_NN.py` files in the folder are the environment's auto-saved scratch cells; the named scripts above are the deliverables.)