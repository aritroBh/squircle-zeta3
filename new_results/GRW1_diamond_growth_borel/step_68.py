wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
result = r'''# GRW1 — the diamond end: growth law of Taylor coefficients, and Borel summability

**Verdict.** Both goals achieved. The growth law is **PROVED** (square, with the two-term refinement, and every rhombus with explicit constant and phase; the paper's missing uniformity step is written out in full and gives an *exponentially* small correction), and Borel summability on both real directions is **PROVED**, with the Borel sum shown equal to I(η) for all η > −1. The paper's Padé–Borel digits are reproduced (direct Padé) and exceeded (H-split Padé). The Stokes constant across arg η = π/2 is 2√π with power ε^{1/2} and phase e^{−3iπ/4}, matching D13c's dictionary; it is also computed as an exact Laplace integral. All reference values of I(η) are **CERTIFIED** (arb) to 61 digits.

Labels: PROVED / CERTIFIED / NUMERICAL / EXACT (exact rational arithmetic). Full proofs in `proof.md` (lemma ledger at the end). Reproduce with `verify_grw1.py` (3.5 s) and the step files.

---

## 1. The growth law (Goal 1) — PROVED

Setup (Lemma 1, PROVED; verified to 30–44 digits at 10 (η,q) pairs, step_02):
  I(η) = ∫₀¹[t(1−t)]^{η/2} g(ηl)dt = ∫₀^∞ e^{−u} H(ηu) K(η,u) du,
  g(z) = √(cosh z + q sinh z), H(z) = e^{−z/2}[g(z)+g(−z)], K(η,u) = (1+e^{−u})^{−2−η}.

Exact coefficient formula (Lemma 2, PROVED; verified 3e-18…3e-39, step_04): with M_{n,j} = ∫₀^∞ e^{−u}u^{n−j}k_j(u)du, M_{n,0} = n!η_D(n), |M_{n,j}| ≤ (n−j)!/(j!(j+1)^{n−j+1}) (Lemma 1(b), PROVED),
  f_n = Σ_{j=0}^n H_{n−j} M_{n,j}   (exact, no uniformity needed at this step).

Darboux at the four branch points ±L ± iπ/2 (Lemma 3, PROVED): local coefficients c = A√2/r at L±iπ/2, c = B√2/r at −L±iπ/2, all real; two-term transfer with remainder O(ρ^{−n}n^{−7/2}) + O(ρ̃^{−n}), ρ̃ < √(L²+9π²/4).

**Uniformity step (Lemma 4, PROVED, explicit constants)** — the one the paper left unwritten:
  |f_n − H_n n! η_D(n)| ≤ C_H e^{2ρ} ρ^{−n} n! 2^{−n−1}/n   (square: C_H = 5.874, EXACT for n ≤ 190),
so the K-correction is O(n^{1/2}2^{−n}) relative to the main term — exponentially small, far below the required O(n^{−1/2}).

**Theorem 1 (PROVED).**
- Square:  f_n = √2 cos(π(3−2n)/4)·n!(2/π)ⁿ n^{−3/2}·(1 + (3−3π(−1)ⁿ)/(8n) + O(n^{−2})).
  Sign pattern (−,+,+,−) ✓. |cos| ≡ √2/2, so no cancellation.
- Rhombus (A ≠ B): f_n ∼ −(√2/(r√π)) ρ^{1/2} n! ρ^{−n} n^{−3/2}·[B(−1)ⁿ cos((n−1/2)θ₀) − A sin((n−1/2)θ₀)], θ₀ = arg(L+iπ/2).
- Radius zero: limsup |I_n/n!|^{1/n} = 1/ρ, ρ = √(log²(B/A)+π²/4) (paper's thm:zero).

Numerical confirmation:
- ratios f_n/law, all four classes mod 4, n ≤ 180: → 1 (e.g. n=180: 0.995415); Richardson constants 0.9999962 in every class (paper had ~6e-6 at n ≤ 64; we match at n ≤ 180);
- f_n/(n!H_n) = 1 − O(2^{−n}) measured down to 6.6e-55 (n=180);
- two-term law for H_n (EXACT rationals): residual ≤ 13.07·(2/π)ⁿn^{−5/2} (n ≤ 190), i.e. O(n^{−2}) with const ≈ 4;
- rhombus (q = −3/5, A = 1, B = 2): ratios → 1 at all n with |bracket| ≥ 1 (99/141 of n ∈ [100,240]); two-term law worst residual 4.8e-4.

## 2. Borel summability (Goal 2) — PROVED

B(ζ) = Σ f_n ζⁿ/n! (radius ρ by Theorem 1). Write H^{[j]}(w) = Σ_m H_m w^{m+j}/(m+j)!.

- Lemma 5 (PROVED): H(z) = √2√{1+e^{−2z}} = √2 Σ c_k e^{−2kz}, c_k = binom(1/2,k), Σ|c_k| = 2; for x ≥ 0, H^{[0]}(x) = √2 Σ c_k J_0(2√{2kx}), |H^{[0]}| ≤ 2√2. (Confirmed to 2.4e-10, step_03.)
- Lemma 6 (PROVED): for x ≤ 0, H^{[0]}(x) = √2 Σ c_k J_0(2√{(1−2k)x}) = √2[I_0(2√{|x|}) + Σ_{k≥1}c_k J_0(2√{(2k−1)|x|})], so |H^{[0]}(x)| ≤ √2(I_0(2√{|x|})+1) ≤ √2(e^{2√{|x|}}+1). (Confirmed to ~2e-10 … 7e-10 at x = −1…−200, step_38.)
- Lemma 7 (PROVED): 𝔅_j(ζ) = ∫₀^∞ e^{−u}k_j(u)u^{−j}H^{[j]}(ζu)du are analytic in strips around ℝ; B = Σ_{j≥0}𝔅_j on |ζ| < ρ, so the Σ𝔅_j is the analytic continuation of B along arg ζ = 0 and arg ζ = π, with |B(ζ)| ≤ C_ε e^{(1+ε)|ζ|} (ℝ⁻) and ≤ Ce^{2√|ζ|} (ℝ⁺) — at most exponential growth.
- Lemma 8 (PROVED): ∫₀^∞ e^{−t}H^{[j]}(wt)dt = w^j H(w) for all real w.
- **Theorem 2 (PROVED):** S(η) = ∫₀^∞ e^{−t}B(ηt)dt = I(η) for all η > −1 (η ≠ 0), S(0) = 1 = I(0); Gevrey-1 asymptotic Σ f_n ηⁿ with rate (2/π+ε); singular directions are arg ζ ∈ {±θ_k, ±(π−θ_k)}, θ_k = arctan(π(2k+1)/(2L)) — for the square exactly arg ζ = ±π/2; **the real axis is not a Stokes line** (Padé poles confirmed beyond π/2, innermost 1.57175081, none within 0.05 of ℝ, step_28).

Numerical confirmation:
- references I(η) **CERTIFIED** by arb at η = ±0.05, ±0.1, ±0.2, ±0.3, total radius ≤ 5.4e-62 (step_34), matching independent DE quadrature to 55+ digits;
- direct Padé [90/90] (build certified 2.25e-358): 44.6 / 31.3 / 25.4 / 23.8 digits at η = 0.1 / 0.2 / 0.3 / −0.2 — reproduces the paper's 44 / 31 / 25 at the positive points (paper claims 30 at −0.2; we measure 23.8 with the direct Padé, and 42.7 with the H-split);
- H-split Padé (R = B − H has radius 2ρ = π, measured 0.305→0.318): **57.0 / 43.7 / 35.4 / 42.7** digits at 0.1 / 0.2 / 0.3 / −0.2, and 55.6 / 54.3 / 34.0 at ±0.05 / −0.1 / −0.3 (steps_20–21);
- Borel sum converges to I(η) on the whole range η ∈ (−1,∞): ±0.5: 26.9/24.6 digits, ±0.7: 22.4/18.8, ±0.9: 19.4/5.4 (step_29);
- B(ζ) on the real axis via the independent Σ𝔅_j integrals matches series/Padé to 5e-18 … 2e-9 at ζ = ±0.5, ±1, ±2, ±3 (step_23).

## 3. Stokes constant across arg η = π/2 (Goal 3) — PROVED formula + NUMERICAL

Theorem 3 (proof.md Part V): with the cut going up from iπ/2, ζ = i(π/2+s): Disc B = −2H = −4√{sin s}e^{i(π/4−s/2)} (left minus right), so
  J(ε) = S_left(iε) − S_right(iε) = −(4i/iε)∫₀^π e^{−(π/2+s)/ε}√{sin s}e^{i(π/4−s/2)}ds·(1+O(e^{−π/ε}))
       = 2√π e^{−3iπ/4} ε^{1/2} e^{−π/(2ε)} (1 − 3iε/4 + O(ε²)).
D13c's dictionary entry confirmed in every slot: |S| = 2√π, exponent b = 1/2 (not 3/2 — the Borel–Laplace 1/η), rate A = π/2.

Measured (steps_24–27): |J|/|pred| = 0.9803, 0.9889, 0.9950, 1.0003 at ε = 0.2, 0.15, 0.1, 0.075 (→ 1 with the (1−3iε/4) correction matching: predicted phase −3ε/4 = −0.15, −0.1125, −0.075, −0.056 vs measured −0.147, −0.111, −0.075, −0.052); the R = B−H jump measured ≤ 1.7e-7 … 4.3e-12 ≈ e^{−π/ε}, negligible as predicted; the right lateral sum equals I(iε) (direct quadrature at complex η) to 2e-7 … 4e-12.

---

## Evidence table

| check | result | residual/error | digits | step |
|---|---|---|---|---|
| I = ∫e^{−u}H(ηu)K du, 10 (η,q) pairs | ✓ | 5e-44 … 8e-31 | 30–44 | 02 |
| H(z) = √2√{1+e^{−2z}} (values) | ✓ | 0 … 2e-81 | 81 | 03 |
| H² ≡ 2(1+e^{−2z}) as formal series to order 190 | ✓ EXACT | 0 | exact | verify[1] |
| H^{[0]}(x) = √2Σc_kJ_0(2√{2kx}), x ≥ 0 | ✓ | 2.4e-10 | 9–10 | 03 |
| H^{[0]}(x) = √2Σc_kJ_0(2√{(1−2k)x}), x ≤ 0 | ✓ | 2e-10 … 7e-10 | 9–10 | 38 |
| f_n = Σ_j H_{n−j}M_{n,j}, n = 1…34 | ✓ | 3.5e-18 … 3.4e-39 | 18–39 | 04 |
| M_{n,0} = n!η_D(n) | ✓ | 2e-18 … 2e-43 | 18–43 | 04 |
| tail bound |M_{n,j}| ≤ (n−j)!/(j!(j+1)^{n−j+1}) | ✓ (ratio ≤ 1) | — | 04 |
| f_n/law → 1, n ≤ 180, 4 classes | ✓ | Richardson const 0.9999962 | ~6e-6 | 05, verify[3] |
| f_n/(n!H_n) = 1−O(2^{−n}) | ✓ | 9.5e-13 (n=40) … 6.6e-55 (n=180) | 12–55 | 05 |
| two-term Darboux residual | ≤ 13.07 (2/π)ⁿn^{−5/2} | O(n^{−2}) const ≈ 4 | — | 11 |
| rhombus law (q=−3/5) | ✓ at |bracket|>1 | worst 2% (1-term), 4.8e-4 (2-term) | — | 06–11 |
| Padé [90/90] reproduces b_0..b_180 | ✓ | 2.25e-358 | 358 | 12 |
| references I(η), 8 points | CERTIFIED (arb) | radius ≤ 5.4e-62 | 61 | 34 |
| Padé–Borel digits at 0.1/0.2/0.3/−0.2 | 44.6/31.3/25.4/23.8 (direct), 57.0/43.7/35.4/42.7 (H-split) | paper: 44/31/25/30 | — | 18–21 |
| S(η) = I(η) on (−1,∞) | ✓ | 5.4–26.9 digits at ±0.5…±0.9 | 5–27 | 29 |
| B(ζ) = Σ𝔅_j vs Padé/series | ✓ | 5e-18 … 2e-9 | 9–18 | 23 |
| Padé poles all |z| > π/2 | ✓ innermost 1.57175081 | none near ℝ | — | 28 |
| Stokes jump |J|/|pred| | 0.9803, 0.9889, 0.9950, 1.0003 | → 1 | 25–27 |
| jump phase vs −3π/4 | Δ = −0.147, −0.111, −0.075, −0.052 | → 0 (−3ε/4 law ✓) | — | 25–27 |
| S_right(iε) = I(iε) | ✓ | 2.3e-7 … 3.7e-12 | 7–12 | 25–27 |
| f450 vs D11 exact ζ-forms (n ≤ 6) | ✓ | 0 … 1.2e-91 | 91 | 35 |

## Files

- `proof.md` — the full proof with lemma ledger (PROVED/CERTIFIED/OPEN; nothing OPEN).
- `verify_grw1.py` — consolidated rerun (< 5 s): exact H_n and H² identity, Darboux laws + explicit constants, f_n ratios and Richardson constants, tail bound, Laplace identity, Padé–Borel [30/30] vs certified reference (22.9 digits), singular directions.
- step files: `step_02` (identity), `step_03` (J_0 form ℝ⁺), `step_04` (f_n identity + tail), `step_05` (growth ratios), `step_06/07` (rhombus engine), `step_08–11` (Darboux 1/n term, two-term laws), `step_12` (Padé build), `step_14–17` (references), `step_18` (direct Padé–Borel), `step_20–21` (H-split), `step_23` (B = Σ𝔅_j cross-check), `step_25–27` (Stokes jump), `step_28` (poles), `step_29` (S = I on (−1,∞)), `step_34` (arb certified refs), `step_35` (f450 vs ζ forms + ℝ⁻ growth), `step_38` (J_0 form ℝ⁻).
- data: `pade9090.json`, `padeR9090.json`, `pade7070.json`, `iref_de160.json`, `iref_certified.json`, `pb_sums.json`, `pb_split_sums.json`, `stokes_jump.json`.
- reused: D11b2 `step2_num450_N180.json` (f₀…f₁₈₀ at 450 digits; cross-checked against D11's exact ζ-forms to 91 digits, step_35).
'''
with open(wd+"/result.md","w") as f:
    f.write(result)
print("result.md written,", len(result.splitlines()), "lines")
