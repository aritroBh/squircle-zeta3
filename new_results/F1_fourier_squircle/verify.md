# F1 (Fourier meets the squircle): adversarial verification

Verifier scratch folder: `<scratch>/verify_F1_fourier_squircle/`. Copies of every verifier script and
its output are in `F1_fourier_squircle/verify/`. Python: the project `.venv`. Every number below is printed by one of
those scripts.

**Verdict: CONFIRMED WITH CORRECTIONS.**
- Every mathematical statement that result.md labels PROVED holds.
- Every NUMERICAL statement reproduces, and several were extended.
- The corrections are one sign typo, one overstated "PROVED given 3b", and two short gaps that are easy to close. One
  of those gaps, the cubic Parseval step, I close below with a cleaner argument. No claimed constant is wrong.

---

## 1. Re-runs of the solver's scripts

I copied all six scripts f1_01 … f1_06 and ran them unchanged. The only edit was f1_06's figure path, so that the
project PNG is not overwritten. **All six outputs are identical, byte for byte, to the stored `.out` files.** So
every residual quoted in result.md is reproduced: 1e-50 for the duality forms, ≤4e-71 for the endpoint coefficients,
c₂ and c₃ to 1.6e-61, the τ⁴/τ⁵ PSLQ relations, D11 weight homogeneity for n ≤ 10, the cubic Parseval to 1.7e-26 with
tail-fit stability 2.8e-22, ζ(3…11) to 1.6e-61, and the comb residual/τ⁴ ≈ −3.9.

## 2. Independent checks (different method or parametrisation, ≥ 40 digits)

| claim | independent method | result | script |
|---|---|---|---|
| Diamond theorem: c₂ = ζ(2), c₃ = ζ(3), c₄ = −(33/16+7k²)ζ(4), normalised by (−1)^kπ/√2 | Direct θ-quadrature of (8/π)∫₀^{π/4} cosθ cos4kθ·expm1(τ log1p(tan^qθ)), with no s-substitution and no Watson coefficients. Then 16-level Richardson/Neville extrapolation τ → 0 (τ₀ = 0.02/(1+k), down to 8e-8). dps 90. | For k = 0, 1, 2, 3, 5, 7: \|c₂ − ζ(2)\| ≤ 6.4e-53, \|c₃ − ζ(3)\| ≤ 2.2e-52, \|c₄ − claim\| ≤ 8.5e-51 (extrapolation stability ≤ 1.1e-50). | v1 |
| c₅(k) = −[(3+16k²)ζ(5) + (7+16k²)ζ(2)ζ(3)]/4 (NUMERICAL in result.md, k ≤ 3) | Same Richardson run. | \|c₅ − claim\| ≤ 2.7e-50 for k = 0, 1, 2, 3, **and also for k = 5, 7**, which were not tested before. | v1 |
| k²-part of c₅ ⇔ ∫₀¹ log²x log²(1+x)/x = 2ζ(2)ζ(3) − (29/8)ζ(5) | quadrature | diff 5.7e-42 | v6 |
| Classical I_{0,3}, I_{1,2} quoted in the τ⁴ proof | quadrature vs the quoted closed forms | 1.5e-50, 1.3e-51. c₄ assembled from them: ≤ 1.7e-49 | v4 |
| Circle end P₀…P₃ (support-function route) | **Polar arclength** P = 8∫₀^{π/4}√(ρ² + ρ′²), at 44 Chebyshev nodes in ε ∈ [−0.1, 0.1], Taylor coefficients taken from the interpolant. dps 60. | P₀…P₃ match the closed forms to ≤ 1.1e-52. P₄ and P₅ match circle_end's stored values to 7.6e-51 and 9.1e-49. | v2 |
| D₂/(2π²) = 7π²/48 − 1 and D₃/(4π²) = −(15/32)ζ(3) + 43π²/384 − (7/48)π²log2 + log2 − 5/16 | P from v2 together with the exact Gamma area | diffs 4.4e-55 and 3.5e-53 | v2 |
| Cubic Parseval Σ(16k²−1)α_kβ_k = D₃/(4π²) | **Single integral, no series and no tail fit**: Σ(16k²−1)α_kβ_k = 2·avg(f′g′ − fg) + α₀β₀/2, with f = κ₁, g = κ₂ + κ₁²/2, and f′, g′ written out by hand and checked against mp.diff. | 1.67e-51 against the closed form. PSLQ on {S, ζ(3), π², π²log2, log2, 1} at 50 digits gives [−384, −180, 43, −56, 384, −120], residual 3e-49 (max coefficient 384, 6 elements: acceptable). The ε² analogue matches 7π²/48 − 1 to 4.0e-51. | v3 |
| β_k explicit formula (solver checked k ≤ 6 directly) | oscillatory quadrature at k = 50 and 200 | 9.3e-43, 6.2e-45 | v6 |
| Euler sums ΣH_k/k² = 2ζ(3), ΣO_{2k}/k² = (7/4)ζ(3) | integral representations | 5e-40, 3e-40 | v7 |
| ζ(3) = (4/π)∫₀^{π/2}u[(4θ−π)² + 8 − π²/3]dθ and ζ(5) = 8⟨u, β₄ − 4β₂⟩ | quadrature in u = tanθ, split at 1 and mirrored | 2.7e-51, 2.7e-51 | v4 |
| a_k(v), k = 1..4 (PSLQ values) and the w-integral formula | θ-quadrature vs the w-formula vs the closed forms | w-formula diff 0 at 50 digits; closed forms ≤ 9e-52; a₀(v) exact | v4 |
| Prop. 3b (uniform 1/sinh form) | Composite Gauss–Legendre in θ (double precision) over **all** k ≤ 400 (the solver sampled only 8 values of k) | max error / τ³ = 0.5185, 0.5304, 0.5359 at τ = .02, .01, .005, always at k = 1. This reproduces the solver's numbers. | v5 |
| Lemma ∫₀^∞(1/y − csch y)² = 2log2 − 1 | the sine transform used in its proof, checked by quadosc | ≤ 1.4e-42 | v6 |
| §6 Minkowski mix: P(t) = 4√2(1−t) + 8t, A(t) = 2 + 4t − 2t² | support-function quadrature | exact to 1e-50 | v4 |

## 3. Proof reading

- **Theorem 1 (P = πa₀(h_p), three integral forms).** Correct. I re-derived the substitution u = tanθ = e^{−τs}; the
  k = 0 kernel is e^{−x}(1+e^{−2x})^{−3/2} = cos²θ sinθ. The dual-norm/Cauchy step is classical. The M₀ "three
  lines" are correct: I re-derived ∫(c+s)log(c+s) = 2√2L − 2.
- **Theorem 3 (diamond, every k). PROVED, and I agree.**
  - The representation is right: π/4 − θ = gd(x)/2, so cos4kθ = (−1)^k cos(2k gd x) and
    cos²θ sinθ = 2^{−3/2}e^{x/2}sech^{3/2}x.
  - Hand checks: g₀ = (−1)^k√2/4, g₁ = (−1)^k√2/8, g₂ = −(−1)^k(5√2/32 + √2k²/2).
  - c₂ and c₃ are as stated.
  - For c₄, I redid the cancellation by hand. The Li₄(½), λ⁴, ζ(2)λ² and ζ(3)λ terms cancel and leave
    ζ(4)(1/64 − 35/128 − 7k²/8), which is −(33/16 + 7k²)ζ(4)·(√2/π)(−1)^k.
  - The remainder argument is fine:
    - 0 ≤ (1+e^{−s})^τ − 1 ≤ 2τe^{−s};
    - |G_k − Taylor_M| ≤ C_M x^M on all of [0, ∞), because G_k is bounded and analytic at 0.
  - Pure weight is PROVED: g_{k,m} ∈ (−1)^k√2·ℚ[k²] with deg ≤ ⌊m/2⌋, and I_{m,j} is an alternating MZV of weight m+j+1.
- **Corollary 3a (Dirac comb).** Correct. I redid the two-sided computation: the φ′ terms cancel by the reflection
  θ ↦ π/2 − θ about each normal, and each side contributes (√2/4)φ(θ_j)(ζ(2)τ² + ζ(3)τ³)/2. Only φ ∈ C² is needed.
- **Re-expansion in d.** τ = d − d² + …, which gives P/√2 = 4 + ζ(2)d² + (ζ(3) − 2ζ(2))d³. That agrees with the ledger
  (R2/D11) P‴(1) = 6√2ζ(3) − 2√2π².
- **Section 4, cumulants.** log h = κ₁ε + κ₂(ε²+ε³) + ⅔κ₃ε³ + O(ε⁴) is correct. I re-derived it from the CGF of the
  two-point law with t = ε/(1−ε).
- **u, α_k.** The coefficients −1/(4k(4k²−1)) are correct (re-derived from the log-cos series).
- **ζ(2j+1) identities.** The algebra Z_j = 4Z_{j−1} − ζ(2j+1) is correct, and so is β₂ = ((4θ−π)² − π²/3)/4.
- **Area lemma.** ∫(1/y − csch y)² = 2log2 − 1: the proof chain is correct (sine transform, Parseval,
  ∫₀^∞(1 − tanh x)² dx = 2log2 − 1).
- **η-quotient asymptotics.** Σ1/(k sinh 2πkτ) = −2 log Π_{m odd}(1 − q^m), and the modular transformation of
  (q;q)_∞/(q²;q²)_∞ gives π/(12τ) − log2 + πτ/6 + O(e^{−c/τ}). Correct.

## 4. Corrections

1. **Sign typo (§2, "Derivative at p = 1 of log ρ_p").**
   - result.md says a_k(v) ∈ (−1)^kπ/4 + ℚ + ℚ/π. It should be **(−1)^{k+1}π/4**.
   - result.md's own table shows this: a₁ = π/4 − 1, a₂ = −π/4 + …. Checked in v4.
2. **"Consistency check (PROVED given 3b)" for the deficit slope −4/π is overstated.**
   - 3b only gives an error O(τ³) uniform in k.
   - Inserted into Σ(16k²−1)(a_k² − A_k²), the cross term Σ 32k²|A_k|·O(τ³) is a divergent sum of O(τ³) terms, so it
     is not controlled without an error bound that decays in k.
   - The computation reproduces −4/π exactly, but it is a **heuristic consistency check**, not a proof.
   - The value −4/π itself is trivially exact from P and A.
   - The lemma ∫(1/y − csch y)² = 2log2 − 1 used in it is PROVED.
3. **Cubic Hurwitz–Parseval: the "modulo dominated convergence" gap is closable, so it can be labelled PROVED.**
   Cleaner route, avoiding termwise differentiation of the Fourier series:
   - Write h = 1 + εf + ε²g + ε³k₃ + …. Differentiate P = ∫h and A = ½∫(h² − h′²) under the integral sign. This is
     justified because h and h′ are smooth in ε near 0, with dominating functions polynomial in |log|c||, |log|s||.
   - That gives D₃ = 2∫f∫g − 4π∫(fg − f′g′) (integrals over [0, 2π]).
   - Parseval for f and g, which holds since f′, g′ ∈ L², gives exactly 4π²Σ(16k²−1)α_kβ_k = D₃.
   - This is the integral form evaluated in v3, which agrees to 1.7e-51.
4. **Derivatives vs asymptotic expansion (minor gap).**
   - Statements like d²a_k/dp² = (−1)^k√2π/3 and P‴(1) are read off from an asymptotic expansion in τ. That needs
     a_k(τ) ∈ C^∞[0, τ₀).
   - This follows by differentiating the s-representation in τ under the integral sign; each derivative has its own
     Watson expansion. One line, not written in result.md.
   - For k = 0, two-sided smoothness at p = 1 is already known from R2/D11.
5. **Prop. 3b wording.** "sup over k ≤ 200" was really a sup over 8 sampled k. I confirm the statement over all
   k ≤ 400 (v5). The O(τ³) error is exactly the missing ζ(3)τ³ term, since √2ζ(3)/π = 0.541 against the observed 0.53.
6. **Cosmetic.**
   - f1_02's printed "check" line for a₀(log ρ₁) has a sign typo; the solver flagged this itself. The PSLQ relation
     and the stated value log2 − 4G/π are correct.
   - f1_06 writes its figure to a hard-coded project path.

## 5. Novelty and literature

- The solver labels the classical ingredients correctly: the support function = dual norm, Cauchy's formula,
  Hurwitz 1901, and the Euler-type log-cos integrals behind ζ(2j+1) = 8⟨u, β_{2j} − 4β_{2j−2}⟩.
- For k = 0 the diamond theorem re-derives P″(1) and P‴(1), which are already in the ledger (R2/D11). Its value there
  is the pure-weight τ form and an independent derivation.
- What looks new:
  - the all-k statement, with the same ζ(2)τ² + ζ(3)τ³ in every support Fourier coefficient;
  - the Dirac-comb reading;
  - pure weight in τ;
  - the explicit cubic Parseval sum with −(15/32)ζ(3).
- I ran three more web searches: support function / superellipse / Fourier; ℓ_p unit ball perimeter near p = 1 with
  ζ(3); support function of ℓ_p ball / mean width / Fourier coefficients. They found nothing that preempts these
  results. The hits were Rodríguez–Montilla 2607.09048 (hypergeometric, no support function), a Blaschke–Santaló
  paper using P = 2πa₀ for generic convex sets, and constant-width papers.
- This is still a limited search, so "to our knowledge new" is the right wording. Mean widths of ℓ_p^n balls are
  standard in asymptotic convex geometry, but only as high-dimensional estimates, not as this 2-D expansion.

## 6. Status table after verification

| claim | solver label | verifier |
|---|---|---|
| P = πa₀(h_p); M₀ in three lines; P′(2) | PROVED | CONFIRMED |
| Endpoint closed forms (h₁, h_∞, ρ₁, log ρ₁, u, γ_m, ∂_pγ_m) | PROVED | CONFIRMED (re-run ≤ 4e-71) |
| a_k(v) structure / values | structure PROVED / NUMERICAL | CONFIRMED, with a sign correction: (−1)^{k+1}π/4 |
| Diamond theorem c₂, c₃, c₄ (all k) | PROVED | CONFIRMED (proof re-derived; independent Richardson ≤ 8.5e-51, k ≤ 7) |
| c₅(k) | NUMERICAL | CONFIRMED NUMERICAL, extended to k = 5, 7 (2.7e-50) |
| Pure weight in τ; D11 homogeneity n ≤ 10 | PROVED | CONFIRMED |
| Dirac comb corollary | PROVED (sketch) | CONFIRMED (needs only C²) |
| Uniform 1/sinh form (3b) | NUMERICAL | CONFIRMED over all k ≤ 400 |
| Deficit slope −4/π from 3b | "PROVED given 3b" | downgrade to heuristic consistency check (exact value matches) |
| P₀…P₃ at the circle | PROVED | CONFIRMED (polar-arclength Chebyshev, ≤ 1.1e-52) |
| Σ(16k²−1)α_k² = 7π²/48 − 1 | PROVED | CONFIRMED |
| Cubic Parseval = −(15/32)ζ(3) + … | PROVED mod DCT | upgrade to PROVED (argument in §4.3); integral form 1.7e-51 |
| ζ(2j+1) = 8⟨u, β_{2j} − 4β_{2j−2}⟩ | PROVED (Euler-type) | CONFIRMED |
| Linear (Minkowski) homotopy stays locked | (derivation) | CONFIRMED |

Bottom line: F1's three questions are answered. (a) The endpoint coefficients are in closed form. (b) ζ(3) appears
through explicit Fourier identities: in every support coefficient at the diamond end, and in the cubic
Hurwitz–Parseval sum at the circle. (c) P″(1) and P‴(1) are read off a₀. The supporting mathematics survives
independent checking at ≥ 49 digits.
