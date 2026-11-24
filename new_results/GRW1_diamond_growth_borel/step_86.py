wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
result = r'''# GRW1 — the diamond end: growth law of Taylor coefficients, and Borel summability

**Verdict.** Both goals achieved. The growth law is **PROVED** (square with the two-term refinement and all-n explicit constants; every rhombus with explicit constant and phase — the paper's missing uniformity step is written out in full and gives an *exponentially* small correction). Borel summability on both real directions is **PROVED**, with the Borel sum shown equal to I(η) for all η > −1 (rhombus: |η| < ρ²/|L|, which is ≥ π). The paper's Padé–Borel digits are reproduced at η = 0.1/0.2/0.3 (44.6/31.3/25.4 vs paper's 44/31/25) and exceeded with the H-split (up to 57 digits); the paper's "30 digits at η = −0.2" is **not** reproduced by the plain [90/90] Padé (we measure 23.8 at every order 60–90, two quadratures; decision in §2) but is exceeded by the H-split (42.7). The Stokes constant across arg η = π/2 is 2√π with power ε^{1/2} and phase e^{−3iπ/4}, matching D13c's dictionary, and is also computed as an exact Laplace integral. All reference values of I(η) are **CERTIFIED** (arb) to 61 digits.

Labels: PROVED / CERTIFIED / NUMERICAL / EXACT (exact rational arithmetic). Full proofs in `proof.md` (lemma ledger at the end). Reproduce with `verify_grw1.py` (≈4 s) and the step files.

---

## 1. The growth law (Goal 1) — PROVED

Setup (Lemma 1, PROVED; verified to 30–44 digits at 10 (η,q) pairs, step_02):
  I(η) = ∫₀¹[t(1−t)]^{η/2} g(ηl)dt = ∫₀^∞ e^{−u} H(ηu) K(η,u) du,
  g(z) = √(cosh z + q sinh z), H(z) = e^{−z/2}[g(z)+g(−z)], K(η,u) = (1+e^{−u})^{−2−η}.

Exact coefficient formula (Lemma 2, PROVED; verified 3e-18…3e-39, step_04): f_n = Σ_{j=0}^n H_{n−j} M_{n,j} with M_{n,0} = n!η_D(n), |M_{n,j}| ≤ (n−j)!/(j!(j+1)^{n−j+1}).

**Delta-domain Darboux transfer with explicit remainder (Lemma 3.2, PROVED; this is the reviewer's item 1–2).** H analytic in the Delta domain {|z| < 2.9ρ} slit from z± = ±iπ/2 away from 0; local form H(z±+w) = 2w^{1/2}(1 − w/2 + w²E(w)) with |E| ≤ E* = 5/24 (CERTIFIED on a 4000-point grid on both slits, steps_39–41). Hankel expansion gives
  H_n = −√2(2/π)ⁿn^{−3/2}cos((2n+1)π/4)(1 + (3−3π(−1)ⁿ)/(8n)) + R_n,
  |R_n| ≤ 4.45·ρ^{−n}n^{−7/2} + 135·(2.9ρ)^{−n}   (n ≥ 10),
and, with the EXACT computation on n ≤ 190,
  |H_n| ≤ C_H(2/π)ⁿn^{−3/2}, C_H = 5.874 for all n ≥ 1;   |H_n| ≥ 0.8745·(2/π)ⁿn^{−3/2} for all n ≥ 10 (proven ≥ 0.968 for n ≥ 190).
(General rhombus: same form with the four slits from ±L ± iπ/2, c_j ∈ {A√2/r, B√2/r}, per-case E*_j.)

**Uniformity step (Lemma 4, PROVED, explicit constants):**
  |f_n − H_n n! η_D(n)| ≤ C_H e^{2ρ} ρ^{−n} n! 2^{−n−1}/n,
so the K-correction is O(n^{1/2}2^{−n}) relative to the main term — exponentially small, far below the required O(n^{−1/2}).

**Theorem 1 (PROVED).**
- Square:  f_n = √2 cos(π(3−2n)/4)·n!(2/π)ⁿ n^{−3/2}·(1 + (3−3π(−1)ⁿ)/(8n) + O(n^{−2})). Sign pattern (−,+,+,−) ✓. |cos| ≡ √2/2, no cancellation.
- Rhombus (A ≠ B): f_n ∼ −(√2/(r√π)) ρ^{1/2} n! ρ^{−n} n^{−3/2}·[B(−1)ⁿ cos((n−1/2)θ₀) − A sin((n−1/2)θ₀)], θ₀ = arg(L+iπ/2).
- Radius zero: limsup |I_n/n!|^{1/n} = 1/ρ, ρ = √(log²(B/A)+π²/4) (paper's thm:zero).

Numerical confirmation:
- square: ratios f_n/law → 1 in all four classes n ≤ 180 (n=180: 0.995415); Richardson constants 0.9999962 in every class; f_n/(n!H_n) = 1 − O(2^{−n}) down to 6.6e-55 (n=180); two-term residual ≤ 13.07·(2/π)ⁿn^{−5/2} (weaker than the proven 4.45·n^{−7/2}-form for n ≥ 1);
- rhombus, three independent cases (reviewer's item 4 — note f_n(A,B) = f_n(B,A) by the swap symmetry, verified to 40 digits, so A=2,B=1 ≡ A=1,B=2):
  q = −3/5 (A=1,B=2): two-term law worst |1−f_n/(n!·law2H)| = 4.8e-4 at safe n ≤ 240;
  q = 5/13 (A=3,B=2): worst 1.8e-3 at safe n ∈ [60,125], residual ×n² ≈ 0.81;
  q = −4/5 (A=1,B=3): worst 4.9e-3, ×n² ≈ 1.27;
  one-term law worst |ratio−1| = 0.02 / 0.039 / 0.083; root tests approach 1/ρ slowly (as predicted).

## 2. Borel summability (Goal 2) — PROVED

B(ζ) = Σ f_n ζⁿ/n! (radius ρ by Theorem 1). Write H^{[j]}(w) = Σ_m H_m w^{m+j}/(m+j)!.

- Lemma 5 (PROVED): H(z) = √2√{1+e^{−2z}} = √2 Σ c_k e^{−2kz}, c_k = binom(1/2,k), Σ|c_k| = 2; for x ≥ 0, H^{[0]}(x) = √2 Σ c_k J_0(2√{2kx}), |H^{[0]}| ≤ 2√2. (Confirmed 2.4e-10, step_03.)
- Lemma 6 (PROVED): for x ≤ 0, H^{[0]}(x) = √2 Σ c_k J_0(2√{(1−2k)x}) = √2[I_0(2√{|x|}) + Σ_{k≥1}c_k J_0(2√{(2k−1)|x|})], so |H^{[0]}(x)| ≤ √2(I_0(2√{|x|})+1) ≤ √2(e^{2√{|x|}}+1). (Confirmed ~2e-10…7e-10 at x = −1…−200, step_38.)
- Lemma 6b (PROVED, indicator): |H^{[0]}(w)| ≤ C_ε e^{(2/π)|Im w|+ε|w|} (rotated-rectangle proof written out).
- Lemma 7 (PROVED): 𝔅_j(ζ) = ∫₀^∞ e^{−u}k_j(u)u^{−j}H^{[j]}(ζu)du is analytic in the open strip |Im ζ| < (j+1)π/2; B = Σ_{j≥0}𝔅_j on |ζ| < ρ, so Σ𝔅_j is the analytic continuation of B along arg ζ = 0 and arg ζ = π; B − H analytic in |Im ζ| < π; |B(ζ)| ≤ Ce^{2√|ζ|} (ℝ⁺), ≤ C(1+√{π|ζ|})e^{|ζ|} (ℝ⁻) — at most exponential growth.
- Lemma 8 (PROVED): ∫₀^∞ e^{−t}H^{[j]}(wt)dt = w^j H(w) for all real w.
- **Theorem 2 (PROVED):** S(η) = ∫₀^∞ e^{−t}B(ηt)dt = I(η) for all η > −1 (η ≠ 0), S(0) = 1 = I(0); Gevrey-1 asymptotic Σ f_n ηⁿ with rate (2/π+ε).

**Singular directions (reviewer's item 3 — final statement with proof status).**
(i) **PROVED: in the open strip |Im ζ| < π, the only singularities of B are the two square-root branch points at ±iπ/2.**
(ii) PROVED decomposition: B(ζ) = H(ζ) − 1 + Σ_{m≥2}(−1)^{m−1}[H(ζ/m) − 2] + R₂(ζ), R₂ = Σ_{j≥1}𝔅_j analytic in |Im ζ| < π; every singular point of the pieces H(ζ/m) lies on the imaginary axis at {±imπ(2k+1)/2}; the point ±iπ is contributed only by H(ζ/2) (square-root branch, consistent with the PROVED asymptotic r_n = [ζⁿ](B−H) = −2^{−n}H_n(1+O(1/n)), radius exactly 2ρ = π).
(iii) The R₂-pieces may a priori have natural boundaries on the lines |Im ζ| = (j+1)π/2 (j ≥ 1); possible off-axis points there are not pinned down by our method.
**Summary:** all identified singular points of B lie on the imaginary axis; singular directions are exactly arg ζ = ±π/2 (square) — rhombus: arg ζ ∈ {±θ_k, ±(π−θ_k)}, θ_k = arctan(π(2k+1)/(2L)) — and **arg ζ = 0, π are not singular directions** (B analytic in strips around both rays with at most exponential growth; rhombus: strips become diamonds (π/2)|Im ζ| + |L||Re ζ| < (j+1)ρ² containing the whole real axis, and S(η) = I(η) holds for |η| < ρ²/|L| ≥ π).

Numerical confirmation:
- references I(η) **CERTIFIED** by arb at η = ±0.05, ±0.1, ±0.2, ±0.3, total radius ≤ 5.4e-62 (step_34);
- direct Padé [90/90] (build certified 2.25e-358): **44.6 / 31.3 / 25.4 / 23.8** digits at η = 0.1 / 0.2 / 0.3 / −0.2, and 39.5 / 15.6 at η = −0.1 / −0.3 (steps_18, 44, 45) — reproduces the paper's 44/31/25 at the positive points;
- **the paper's "30 digits at η = −0.2" (reviewer's item 5 — decision):** the plain Padé–Borel [M/M] gives 23.8 digits at M = 90 and 19.3–22.8 at M = 60–85 (two independent quadratures agree: DE 23.83, adaptive 20.60); so the paper's 30 is not matched by the plain Padé — it presumably reflects a different resummation detail (quadrature/reference/subtraction). With the H-split (B = H + R, R of radius 2ρ = π measured 0.305→0.318) we measure **42.7** at η = −0.2, comfortably exceeding the paper's 30; the claim "the divergent series resums to dozens of digits at η = −0.2" is thus validated, and our direct-Padé number 23.8 is the honest plain-Padé value;
- H-split Padé: **57.0 / 43.7 / 35.4 / 42.7** digits at 0.1 / 0.2 / 0.3 / −0.2; 55.6/54.3 at ±0.05/−0.1; 34.0 at −0.3 (steps_20–21);
- Borel sum converges to I(η) on the whole range η ∈ (−1,∞): ±0.5: 26.9/24.6 digits, ±0.7: 22.4/18.8, ±0.9: 19.4/5.4 (step_29); the task's points η = ±0.1, ±0.3 are covered with 44.6/39.5/25.4/15.6 (direct) and 57.0/54.3/35.4/34.0 (H-split) digits;
- B(ζ) on the real axis via the independent Σ𝔅_j integrals matches series/Padé to 5e-18 … 2e-9 at ζ = ±0.5, ±1, ±2, ±3 (step_23);
- Padé poles all |z| > π/2 (innermost 1.57175081), none within 0.05 of ℝ (step_28).

## 3. Stokes constant across arg η = π/2 (Goal 3) — PROVED formula + NUMERICAL

Theorem 3 (proof.md Part V): with the cut going up from iπ/2, ζ = i(π/2+s): Disc B = −2H = −4√{sin s}e^{i(π/4−s/2)} (left minus right), so
  J(ε) = S_left(iε) − S_right(iε) = −(4i/iε)∫₀^π e^{−(π/2+s)/ε}√{sin s}e^{i(π/4−s/2)}ds·(1+O(e^{−π/ε}))
       = 2√π e^{−3iπ/4} ε^{1/2} e^{−π/(2ε)} (1 − 3iε/4 + O(ε²)).
D13c's dictionary entry confirmed in every slot: |S| = 2√π, exponent b = 1/2 (not 3/2 — the Borel–Laplace 1/η), rate A = π/2.

Measured (steps_24–27): |J|/|pred| = 0.9803, 0.9889, 0.9950, 1.0003 at ε = 0.2, 0.15, 0.1, 0.075 (→ 1); phase deviations −0.147, −0.111, −0.075, −0.052 vs predicted −3ε/4 = −0.15, −0.1125, −0.075, −0.056 (the (1−3iε/4) correction matches); R = B−H jump measured ≤ 1.7e-7 … 4.3e-12 ≈ e^{−π/ε}, negligible; the right lateral sum equals I(iε) (direct quadrature at complex η) to 2e-7 … 4e-12.

---

## Evidence table

| check | result | residual/error | digits | step |
|---|---|---|---|---|
| I = ∫e^{−u}H(ηu)K du, 10 (η,q) pairs | ✓ | 5e-44 … 8e-31 | 30–44 | 02 |
| H(z) = √2√{1+e^{−2z}} (values + formal series) | ✓ EXACT | 0 … 2e-81 | 81/exact | 03, verify[1] |
| H^{[0]}(x) = √2Σc_kJ_0(2√{2kx}), x ≥ 0 | ✓ | 2.4e-10 | 9–10 | 03 |
| H^{[0]}(x) = √2Σc_kJ_0(2√{(1−2k)x}), x ≤ 0 | ✓ | 2e-10 … 7e-10 | 9–10 | 38 |
| E* = max|E(is)| on slits = 5/24 | CERTIFIED (grid) | ≤ 0.20833334 | — | 39–41 |
| f_n = Σ_j H_{n−j}M_{n,j}, n = 1…34 | ✓ | 3.5e-18 … 3.4e-39 | 18–39 | 04 |
| M_{n,0} = n!η_D(n); tail bound ratio ≤ 1 | ✓ | 2e-18 … 2e-43 | 18–43 | 04 |
| f_n/law → 1, n ≤ 180, 4 classes | ✓ | Richardson const 0.9999962 | ~6e-6 | 05, verify[3] |
| f_n/(n!H_n) = 1−O(2^{−n}) | ✓ | 9.5e-13 (n=40) … 6.6e-55 (n=180) | 12–55 | 05 |
| Darboux remainder |R_n| ≤ 4.45ρ^{−n}n^{−7/2}+135(2.9ρ)^{−n} | PROVED constant | — | 3.2 |
| two-term residual, square | ≤ 13.07 (2/π)ⁿn^{−5/2} | O(n^{−2}) const ≈ 4 | — | 11 |
| rhombus two-term law, 3 cases | worst 4.8e-4 / 1.8e-3 / 4.9e-3 | ×n² ≈ 0.5–1.3 | — | 06–11, 42–43 |
| Padé [90/90] reproduces b_0..b_180 | ✓ | 2.25e-358 | 358 | 12 |
| references I(η), 8 points | CERTIFIED (arb) | radius ≤ 5.4e-62 | 61 | 34 |
| Padé–Borel digits at 0.1/0.2/0.3/−0.2 | 44.6/31.3/25.4/23.8 (direct), 57.0/43.7/35.4/42.7 (H-split) | paper: 44/31/25/30; see decision §2 | — | 18–21, 44, 45 |
| direct Padé digits at −0.1/−0.3 | 39.5 / 15.6 | — | — | 18, 45 |
| S(η) = I(η) on (−1,∞) | ✓ | 5.4–26.9 digits at ±0.5…±0.9 | 5–27 | 29 |
| B(ζ) = Σ𝔅_j vs Padé/series | ✓ | 5e-18 … 2e-9 | 9–18 | 23 |
| Padé poles all |z| > π/2 | ✓ innermost 1.57175081 | none near ℝ | — | 28 |
| Stokes jump |J|/|pred| | 0.9803, 0.9889, 0.9950, 1.0003 | → 1 | 25–27 |
| jump phase vs −3π/4 | Δ = −0.147, −0.111, −0.075, −0.052 | → 0 (−3ε/4 law ✓) | — | 25–27 |
| S_right(iε) = I(iε) | ✓ | 2.3e-7 … 3.7e-12 | 7–12 | 25–27 |
| f450 vs D11 exact ζ-forms (n ≤ 6) | ✓ | 0 … 1.2e-91 | 91 | 35 |

## Files

- `proof.md` — the full proof with lemma ledger (PROVED/CERTIFIED/OPEN-with-gap; the only remaining gap: possible off-axis R₂-singularities on |Im ζ| = (j+1)π/2, j ≥ 1, are not pinned down — stated as such).
- `verify_grw1.py` — consolidated rerun (≈4 s): exact H_n and H² identity, Darboux laws + explicit constants, f_n ratios and Richardson constants, tail bound, Laplace identity, Padé–Borel [30/30] vs certified reference (22.9 digits), singular directions.
- step files: `step_02` (identity), `step_03` (J_0 form ℝ⁺), `step_04` (f_n identity + tail), `step_05` (growth ratios), `step_06/07` (rhombus engine), `step_08–11` (Darboux 1/n term, two-term laws), `step_12` (Padé build), `step_14–17` (references), `step_18` (direct Padé–Borel), `step_20–21` (H-split), `step_23` (B = Σ𝔅_j cross-check), `step_25–27` (Stokes jump), `step_28` (poles), `step_29` (S = I on (−1,∞)), `step_34` (arb certified refs), `step_35` (f450 vs ζ forms + ℝ⁻ growth), `step_36–37` (contour diagnostics), `step_38` (J_0 form ℝ⁻), `step_39–41` (E* certified), `step_42–43` (rhombus cases q = 5/13, −4/5), `step_44–45` (paper-30-digit reconciliation).
- data: `pade9090.json`, `padeR9090.json`, `pade7070.json`, `iref_de160.json`, `iref_certified.json`, `pb_sums.json`, `pb_split_sums.json`, `stokes_jump.json`.
- reused: D11b2 `step2_num450_N180.json` (f₀…f₁₈₀ at 450 digits; cross-checked against D11's exact ζ-forms to 91 digits, step_35).
'''
with open(wd+"/result.md","w") as f:
    f.write(result)
import shutil
shutil.copy(wd+"/result.md", wd+"/result_nofinal.md")
print("result.md updated and copied to result_nofinal.md,", len(result.splitlines()), "lines")
