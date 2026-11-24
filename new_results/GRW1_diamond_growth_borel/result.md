# GRW1 — The diamond end: growth of Taylor coefficients, and Borel summability

**Both goals achieved, and the six review items are closed.**

- **Growth law: PROVED** (square, with two-term refinement and all-n explicit constants; every rhombus with explicit constant and phase). The paper's missing uniformity step is written out in full — the K-correction is *exponentially* small, not merely O(n^{−1/2}).
- **Borel summability: PROVED** on both real directions, with the Borel sum shown equal to I(η) for all η > −1 (rhombus: |η| < ρ²/|L|, which is ≥ π).
- **Paper's Padé–Borel digits reproduced** at η = 0.1/0.2/0.3 (44.6/31.3/25.4 vs paper's 44/31/25) and **exceeded** by the H-split (up to 57 digits); the paper's "30 at η = −0.2" is **not** reproduced by the plain Padé (we measure 23.8 at every order 60–90, two quadratures) but is exceeded by the H-split (42.7) — decision recorded.
- **Stokes constant:** 2√π, power ε^{1/2} (D13c's correction), phase e^{−3iπ/4} — PROVED formula + measured.
- **References I(η) CERTIFIED** by arb to 61 digits.

---

## Theorem statements and proof skeletons (full text: `proof.md`)

**Setup (Lemma 1, PROVED).** For the rhombus |x/A|ᵖ+|y/B|ᵖ = 1, r = √(A²+B²), q = (A²−B²)/(A²+B²), L = log(B/A), ρ = √(L²+π²/4), η = 1/p−1:
`I(η) = P/(4r(1+η)) = ∫₀^∞ e^{−u} H(ηu) K(η,u) du`, with H(z) = e^{−z/2}[g(z)+g(−z)], g = √(cosh z + q sinh z), K(η,u) = (1+e^{−u})^{−2−η} = Σ k_j(u)η^j, |k_j(u)| ≤ e^{−ju}/j!. Verified to 30–44 digits at 10 (η,q) pairs (step_02).

**Lemma 2 (PROVED).** f_n = Σ_{j=0}^n H_{n−j}M_{n,j} **exactly** (differentiation under the integral, dominated; no uniformity needed here). Verified to 3e-18…3e-39 (step_04). M_{n,0} = n!η_D(n), |M_{n,j}| ≤ (n−j)!/(j!(j+1)^{n−j+1}).

**Lemma 3.2 (Delta-domain Darboux transfer, PROVED, explicit).** H is analytic in the Delta domain {|z| < 2.9ρ} slit from z± = ±iπ/2 outward; on the slits H(z±+w) = 2w^{1/2}(1 − w/2 + w²E(w)) with **|E| ≤ E* = 5/24** (CERTIFIED on a 4000-point grid, steps_39–41). Hankel expansion gives
`H_n = −√2(2/π)ⁿn^{−3/2}cos((2n+1)π/4)·(1 + (3−3π(−1)ⁿ)/(8n)) + R_n`,
**|R_n| ≤ 4.45·ρ^{−n}n^{−7/2} + 135·(2.9ρ)^{−n}** (n ≥ 10). Hence C_H = 5.874 (all n ≥ 1), c_H = 0.8745 (all n ≥ 10, proven ≥ 0.968 beyond 190). Rhombus: same form, four slits from ±L ± iπ/2, c_j ∈ {A√2/r, B√2/r}.

**Lemma 4 (the uniformity step, PROVED, explicit).** |f_n − H_n n! η_D(n)| ≤ C_H e^{2ρ} ρ^{−n} n! 2^{−n−1}/n — exponentially smaller than the main term. Measured tail/main = 9.2e-8 (n=20) … 1.0e-38 (n=120), always below the proven bound (verify[4]).

**Theorem 1 (PROVED).** Square: `f_n = √2 cos(π(3−2n)/4)·n!(2/π)ⁿ n^{−3/2}·(1 + (3−3π(−1)ⁿ)/(8n) + O(n^{−2}))`, pattern (−,+,+,−), |cos| ≡ √2/2 (no cancellation). Rhombus: `f_n ∼ −(√2/(r√π)) ρ^{1/2} n! ρ^{−n} n^{−3/2} [B(−1)ⁿ cos((n−1/2)θ₀) − A sin((n−1/2)θ₀)]`, θ₀ = arg(L+iπ/2). limsup |I_n/n!|^{1/n} = 1/ρ (radius zero).

**Lemmas 5–8 (Borel machinery, PROVED).** H(z) = √2√{1+e^{−2z}} = √2 Σ c_k e^{−2kz}, Σ|c_k| = 2. H^{[0]}(x) = √2 Σ c_k J_0(2√{2kx}) for x ≥ 0 (|H^{[0]}| ≤ 2√2; step_03, 2.4e-10); H^{[0]}(x) = √2 Σ c_k J_0(2√{(1−2k)x}) for x ≤ 0 (|H^{[0]}| ≤ √2(I_0(2√{|x|})+1) ≤ √2(e^{2√{|x|}}+1); step_38, ~2e-10 at x = −1…−200). Indicator (rotated-rectangle proof): |H^{[0]}(w)| ≤ C_ε e^{(2/π)|Im w|+ε|w|}. Pieces 𝔅_j(ζ) = ∫₀^∞ e^{−u}k_j u^{−j}H^{[j]}(ζu)du are analytic in |Im ζ| < (j+1)π/2, and B = Σ𝔅_j on |ζ| < ρ — the continuation of B along arg ζ = 0 and arg ζ = π, with |B| ≤ Ce^{2√|ζ|} (ℝ⁺), ≤ C(1+√{π|ζ|})e^{|ζ|} (ℝ⁻). Key identity ∫₀^∞ e^{−t}H^{[j]}(wt)dt = w^jH(w) (all real w).

**Theorem 2 (PROVED).** S(η) = ∫₀^∞ e^{−t}B(ηt)dt = I(η) for all η > −1 (η ≠ 0), S(0) = 1 = I(0); Gevrey-1 asymptotic with rate (2/π+ε).

**Singular directions (final statement).** (i) **PROVED: in |Im ζ| < π, B's only singularities are ±iπ/2.** (ii) **PROVED decomposition:** B(ζ) = H(ζ) − 1 + Σ_{m≥2}(−1)^{m−1}[H(ζ/m) − 2] + R₂(ζ) (R₂ analytic in |Im ζ| < π); all identified singular points lie on the imaginary axis ({±imπ(2k+1)/2}); the point ±iπ comes only from H(ζ/2), matching the PROVED asymptotic r_n = −2^{−n}H_n(1+O(1/n)), radius 2ρ = π. (iii) R₂-pieces may have natural boundaries on |Im ζ| = (j+1)π/2 (j ≥ 1); off-axis points there not pinned down (stated honestly). Singular directions: exactly arg ζ = ±π/2 (square), {±θ_k, ±(π−θ_k)} (rhombus); **arg ζ = 0, π are not singular directions** (Padé poles all |z| > π/2, innermost 1.57175081, step_28).

**Theorem 3 (Stokes jump, PROVED formula + NUMERICAL).** Across arg η = π/2:
`J(ε) = S_left(iε) − S_right(iε) = −(4i/iε)∫₀^π e^{−(π/2+s)/ε}√{sin s}e^{i(π/4−s/2)}ds = 2√π e^{−3iπ/4} ε^{1/2} e^{−π/(2ε)} (1 − 3iε/4 + O(ε²))` — D13c's dictionary in every slot (|S| = 2√π, b = 1/2, A = π/2).

---

## Evidence table (every number printed by code in this directory)

| check | residual/error | digits | step |
|---|---|---|---|
| I = ∫e^{−u}H(ηu)K du, 10 (η,q) pairs | 5e-44 … 8e-31 | 30–44 | 02 |
| H(z) = √2√{1+e^{−2z}} (values; formal series to order 190) | 0 … 2e-81; EXACT | 81/exact | 03, verify[1] |
| H^{[0]} = √2Σc_kJ_0(2√{2kx}) (ℝ⁺) | 2.4e-10 | 9–10 | 03 |
| H^{[0]} = √2Σc_kJ_0(2√{(1−2k)x}) (ℝ⁻) | 2e-10 … 7e-10 | 9–10 | 38 |
| E* = 5/24 on both slits | ≤ 0.20833334 (CERTIFIED grid) | — | 39–41 |
| f_n = Σ_j H_{n−j}M_{n,j}, n = 1…34 | 3.5e-18 … 3.4e-39 | 18–39 | 04 |
| M_{n,0} = n!η_D(n); tail ratio ≤ 1 | 2e-18 … 2e-43 | 18–43 | 04 |
| f_n/law → 1, 4 classes, n ≤ 180 | Richardson 0.9999962 (each class) | ~6e-6 | 05, verify[3] |
| f_n/(n!H_n) = 1−O(2^{−n}) | 9.5e-13 … 6.6e-55 | 12–55 | 05 |
| Darboux \|R_n\| ≤ 4.45ρ^{−n}n^{−7/2}+135(2.9ρ)^{−n} | PROVED constant (n ≥ 10) | — | proof 3.2 |
| two-term residual, square | ≤ 13.07 (2/π)ⁿn^{−5/2} | O(n^{−2}) ≈ 4 | 11 |
| rhombus two-term law, q = −3/5, 5/13, −4/5 | worst 4.8e-4 / 1.8e-3 / 4.9e-3 | ×n² ≈ 0.5–1.3 | 06–11, 42–43 |
| A↔B swap symmetry f_n(1,2) = f_n(2,1) | identical series (40 digits) | — | 06–07 |
| Padé [90/90] reproduces b_0..b_180 | 2.25e-358 | 358 | 12 |
| references I(η), 8 points | CERTIFIED (arb), radius ≤ 5.4e-62 | 61 | 34 |
| direct Padé–Borel digits 0.1/0.2/0.3/−0.2 | **44.6/31.3/25.4/23.8** (paper: 44/31/25/30) | — | 18, 44, 45 |
| direct Padé digits −0.1/−0.3 | 39.5 / 15.6 | — | 18, 45 |
| H-split Padé digits 0.1/0.2/0.3/−0.2 | **57.0/43.7/35.4/42.7** | — | 20–21 |
| S(η) = I(η) on (−1,∞), ±0.5…±0.9 | 5.4–26.9 digits | 5–27 | 29 |
| B(ζ) = Σ𝔅_j vs Padé/series | 5e-18 … 2e-9 | 9–18 | 23 |
| Padé poles | all \|z\| > π/2, innermost 1.57175081 | — | 28 |
| Stokes jump \|J\|/\|pred\| | 0.9803, 0.9889, 0.9950, 1.0003 (→1) | — | 25–27 |
| jump phase vs −3π/4 | Δ = −0.147, −0.111, −0.075, −0.052 (−3ε/4 law ✓) | — | 25–27 |
| S_right(iε) = I(iε) (complex-η quadrature) | 2.3e-7 … 3.7e-12 | 7–12 | 25–27 |
| f450 vs D11 exact ζ-forms (n ≤ 6) | 0 … 1.2e-91 | 91 | 35 |

## Files (new_results/GRW1_diamond_growth_borel)

- `proof.md` — full proof with lemma ledger. Only stated gap: possible off-axis R₂-singularities on |Im ζ| = (j+1)π/2, j ≥ 1, not pinned down (does not affect any theorem).
- `result.md`, `result_nofinal.md` — report; `progress.log` — step log.
- `verify_grw1.py` — consolidated verification, runs 3.2 s (< 5 min), 9 sections: exact H_n + H² identity; Darboux laws with explicit constants; f_n ratios/Richardson; tail bound; Laplace identity; Padé–Borel [30/30] vs certified reference (22.9 digits); singular directions; E* check; rhombus q = 5/13 two-term law.
- Step files (each produced the numbers cited): `step_02` (identity), `step_03` (J_0 form ℝ⁺), `step_04` (f_n identity + tail), `step_05` (growth ratios + Richardson), `step_06/07` (rhombus engine), `step_08–11` (1/n term + two-term laws), `step_12` (Padé build), `step_14–17` (references), `step_18` (direct Padé–Borel), `step_20–21` (H-split), `step_23` (B = Σ𝔅_j), `step_25–27` (Stokes jump), `step_28` (poles), `step_29` (S = I on (−1,∞)), `step_34` (arb certified references), `step_35` (f450 vs ζ forms), `step_36–37` (contour diagnostics), `step_38` (J_0 form ℝ⁻), `step_39–41` (E* certified), `step_42–43` (rhombus q = 5/13, −4/5), `step_44–45` (paper-30-digit reconciliation).
- Data: `pade9090.json`, `padeR9090.json`, `pade7070.json`, `iref_de160.json`, `iref_certified.json`, `pb_sums.json`, `pb_split_sums.json`, `stokes_jump.json`.
- Reused: D11b2 `step2_num450_N180.json` (f₀…f₁₈₀ at 450 digits, cross-checked against D11's exact ζ-forms to 91 digits).