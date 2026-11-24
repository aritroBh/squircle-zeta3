# Verification of H1 (octahedron end): adversarial review

Verifier scripts and outputs are in `<scratch>/verify_H1_octahedron_end/`. The solver's scripts were
re-run in its `rerun/` subfolder, and every number below was printed by one of those scripts. Python is
`.venv/bin/python`. I modified none of the solver's files; this file is the only one I added.

**Verdict: CONFIRMED, with minor corrections to labels and wording (§6).** Every exact statement checked out against
my own derivation and code at ≥ 70 digits, using a different chart from the solver's: T1–T7, the S₁/S₂/S₃ values,
the ζ(3) coefficient, the log form, the general-n η²/η³/η⁴ law, the c₁ formula and the n = 2 diamond law. The
large-order claim T8 is correctly labelled NUMERICAL. I reproduced its singularity location only to about 2–3 digits.

## 1. Re-runs of the solver's scripts

The scripts were copied to `rerun/` and run unchanged. For each pair of files below the result is byte-identical,
except where noted.

* JSON files: `routeA_n2_K12.json`, `routeA_n3_K12.json`, `routeB_n3_K10.json`, `routeA_n4-5-6-7-8_K8.json`,
  `hi_coeffs_n{2,3,4,5}_K80.json` and `direct_S_dps40.json` (via `direct_dps40.out`). The n = 9–12 Route A data
  were regenerated in a single file. Its content equals the solver's per-n files.
* Outputs: `compare_n3_K12_dps40.out`, `pslq_independent.out`, `logform_n3.out`, `logform_general_n.out`,
  `logB_coeffs.json`, `final_numbers.out`, `borel_growth.out`, `growth_exponent.out` and `singular_curve.out`.
* Printed residuals I confirmed:
  * |N1 − N2| ≤ 1.0e−39;
  * series(K = 12) minus direct quadrature: 1.1e−33 at d = ±0.001;
  * Route A ≡ Route B: every difference is 0, through d¹⁰;
  * n = 2 against D11's a₀…a₁₀: every difference is 0;
  * PSLQ (9c₂ − 9 + π² = 0; −9c₃ + 9 + π² − 16ζ(3) = 0): residual 0 at 50 and 60 digits.

## 2. Independent checks (different method and parametrisation)

**V1 (`V1_dirichlet_logmoments.py`): exact and symbolic.**
* Method: expand the simplex integrand pointwise in η as a polynomial in X_i = log t_i, then evaluate every moment
  E[ΠX_i^{k_i}] exactly from the Dirichlet MGF (n−1)!ΠΓ(1+u_i)/Γ(n+Σu_i). This is uncentred, with no
  Gamma-variable or Lukacs step and no finite differences, and it is my own code.
* Coverage: n = 2, 3 to order 10; n = 4, 5 to order 6; n = 6, 7 to order 5.
* Euler's γ cancels in every coefficient.
* `V3_compare_exact.py` compares V1 with the solver's Route A files (original and re-run) and with Route B. Every
  η- and d-coefficient difference is exactly 0.

**V2 (`V2_central_projection.py`): numerical, new chart.** I used radial projection onto the plane u₁ = 1:

  S(p) = 48 ∫₀¹∫₀¹ y (1 + y^p + (yw)^p)^{−(2+p)/p} (1 + y^{2(p−1)} + (yw)^{2(p−1)})^{1/2} dw dy

* The integrand comes from dA = |∇‖u‖_p| ‖u‖_p^{−3} dΩ. I checked it at p = 1 (4√3) and p = 2 (4π, to 1.7e−80).
* Each node carries truncated power series in d. The nested tanh–sinh rule uses step 2⁻⁶, 80-digit arithmetic, and
  step-halving error ≤ 1e−70 for c₀…c₁₀.

| quantity | V2 − claimed exact | V2 error estimate |
|---|---|---|
| c₂ = 1 − π²/9 | −9.2e−82 | 1.7e−80 |
| c₃ = 1 + π²/9 − 16ζ(3)/9 | −9.9e−82 | 5.5e−79 |
| c₄ = 1 − 2π²/9 + 32ζ(3)/9 − 16π⁴/405 | −1.1e−81 | 9.7e−78 |
| c₅ (claimed closed form) | 5.9e−80 | 9.8e−76 |
| c₆…c₁₀ against V1 exact | ≤ 6.9e−77 | ≤ 2.7e−70 |
| c₆, c₇, c₈ against the decimals printed in result.md | 4e−39, 6e−30, 3e−27 (agree to every printed digit) | |
| S(1±0.001), S(1±0.01), S(1±0.1), S(2) against the solver's N1 | ≤ 1.3e−39 (N1 has 40 digits) | ≤ 2e−79 |

* S₂ = 4√3c₂ = −0.669421780076565987998202567942442659218.
* S₃ = −0.2796731272751286206150690832602315974454.
* The ζ(3) coefficient of S at d³ is −64√3/9 (V3, exact).

**PSLQ on the V2 numbers (`V5_check_V2_numbers.py`).** Each relation below appeared unchanged at 50, 60 and 70
digits. Residuals were evaluated at 80 digits.

| constant | relation | residual | max coefficient |
|---|---|---|---|
| c₂ | 9c₂ − 9 + π² = 0 | 0 | 9 |
| c₃ | −9c₃ + 9 + π² − 16ζ(3) = 0 | 0 | 16 |
| c₄ | (405, −405, 90, −1440, 16)·(c₄, 1, π², ζ(3), π⁴) = 0 | 2e−78 | 1440 |
| c₅ | basis {c₅, 1, π², ζ(3), π⁴, π²ζ(3), ζ(5)}, coefficients (−405, 405, 90, −2880, 48, 80, −2400) | 4e−77 | 2880 |

* Each relation is identical to the exact c_k.
* Enlarged basis at 75 digits, run for c₂, c₃, c₄. The extra constants were log 2, log 3, √2, √3, π√3, π²√3,
  √3·L(χ₋₃,2), Catalan, π log 2 and the dihedral angle arccos(−1/3). PSLQ returns the pure single-zeta relation
  each time.
* Note: my first try also included Cl₂(π/3). PSLQ then found the identity Cl₂(π/3) = (3√3/4)L(χ₋₃,2), so the two
  constants must not be in the same basis. The solver's basis contains only one of them, so it is fine.

**Closed-form claims (`V3_compare_exact.py`, `V7_eta5_z2z3.py`, exact, on V1 data).** All of these hold:
* **T5:** L₂…L₇ hold as stated, and L₈–L₁₀ also equal the solver's printed values.
* **Rat₃:** (1−2η)/(1−η) = 1 − d, which confirms Rat₃ = 1/(2−p).
* **T6, n = 2…7:** the η¹ remainder is 0. The η², η³ and η⁴ coefficients equal A₂ζ(2), A₃ζ(3) and W₄.
* **σ_k:** log Rat_n = Σσ_k d^k/k.
* **c₁:** c₁ = (n−1)(H_{n−1} − 1).
* **n = 4, 5 table:** the c₂ and c₃ entries match.
* **n = 2:** reproduces P/√2 = 4 + (π²/6)d² + (ζ(3) − π²/3)d³, and 4c_k = D11's a_k for k ≤ 10.
* **η⁵ observation:** the ζ(2)ζ(3) coefficient of log B_n at η⁵ equals 4(n−1)(n−3)/n³ for n = 2…7.
* **Rational parts:** for n = 3 the rational part of every c_k (k ≤ 10) is 1.

**High orders (`V4_dirichlet_highorder.py`, `V6_highorder_borel.py`).**
* V4 regenerates the 80-order coefficients independently. It uses the Dirichlet form, with denominators
  Γ(n + (2J−n)η) and an EGF in J, at 380 digits. It agrees with the solver's `hi_coeffs` to a relative 4e−120
  for n = 2…5, in both η and d.
* Padé–Borel on the V4 data reproduces the solver's poles. For n = 3, η-series, [40/40]: −0.347372 ± 1.576058i.
* My own estimator does not assume t*. It uses a two-term recurrence fit, so w_k → 1/t*. Results:
  * **n = 2:** Richardson converges to ±1.5707963i, i.e. ±iπ/2, to 8 digits.
  * **n = 3, 4, 5:** the raw t_k at k = 77 are −0.3536 + 1.6009i, −0.5609 + 1.6002i and −0.7082 + 1.5992i. The
    first-order Richardson extrapolations are −0.340 + 1.572i / −0.352 + 1.570i (η / d) for n = 3, −0.541 + 1.568i /
    −0.556 + 1.576i for n = 4, and −0.688 + 1.565i / −0.695 + 1.581i for n = 5. The predicted values are −0.3466,
    −0.5493 and −0.6931 (each + 1.5708i).
  * **Higher-order Richardson diverges for n ≥ 3.** This fits the solver's own remark that a cut starts at t_edge,
    so the 1/k expansion is not clean.
* The pattern Re t* = −log(n−1)/2 shows up consistently across n = 2…5. So T8's location is supported to about
  1e−2. Its exponent and amplitude are not determined, as the solver also says.

## 3. Proof review

* **T1 (PROVED, correct).**
  * I re-derived the octant formula (1) by hand: dA = |∇F|/|F_z|, with x = t^ε.
  * Lukacs step: E[h(G)] = E[R^η]·E_Dir[h] with R ~ Gamma(nε), so the denominator is Γ(1+(n−1)ε). Checked.
  * The p = 2 check gives 2π^{n/2}/Γ(n/2). Checked.
  * Literature: I read pp. 241–244 of Naor–Romik (Ann. IHP 39 (2003)), from the Numdam PDF. Lemma 1 (p. 244) states
    dσ_K/dμ_K = n·vol(K)/area(∂K)·‖∇‖·‖_K‖₂. Theorem 1 (Schechtman–Zinn) represents the cone measure. Together
    these give T1, so the attribution "classical in substance" is correct.
* **T2 (PROVED, correct).**
  * The binomial remainder bound on x ≥ −1 is valid: the Taylor remainder for |x| ≤ ½, and |x|^K ≤ 2|x|^{K+1} for |x| ≥ ½.
  * The mean-value bound |G^{2η} − 1| ≤ 2|η||log G| max(1, G^{2η}) is valid.
  * Moments E[G^c] are finite for c > −ε, which holds for small |η| on both sides.
  * Each E[x^k] is analytic.
* **T3 (PROVED, correct; two small gaps in the exposition).**
  * (a) The m_i are not independent. The "partition into cumulants" step should be stated for the underlying iid
    log G_j, of which the m_i are linear forms. Shift invariance then removes κ₁. The conclusion is unchanged.
  * (b) Route B's termwise expectation is not justified in the text. It is easy to justify, because log G has an
    exponential left tail, so E[e^{c|η||m|}|m|^K] < ∞. Alternatively, it follows from the Route A remainder together
    with uniqueness of asymptotic expansions.
  * The flat-face algebra checks out: Λ_n, Rat_n, the cancellation of γ, and the product form of Rat_n in d.
* **T6 (PROVED, correct for all n).**
  * I re-derived the η⁴ coefficient of log B_n from the standard sample-moment formulas:
    * E m₂ = ((n−1)/n)κ₂;
    * E m₃ = ((n−1)(n−2)/n²)κ₃;
    * Var m₂ = ((n−1)²/n³)κ₄ + (2(n−1)/n²)κ₂²;
    * E m₄ = ((n−1)(n²−3n+3)/n³)μ₄ + (3(n−1)(2n−3)/n³)μ₂².
  * I added the a-dependence: +3(n−1)³ζ(4)/n³ from ψ′ and −4(n−1)²(n−2)ζ(4)/n³ from ψ″.
  * Result: ζ(2)² coefficient −(n−1)/n², and ζ(4) coefficient (n−1)(n²−3n+4)/n³. Both are exactly as claimed.
  * W₄'s third term is Λ_n at j = 4. Checked.
* **Channel table (§5).** −64/27 + 24/27 − 8/27 = −16/9, and for n = 2: −1/4 + 1/2 + 0 = 1/4. Both checked.
  * E[m₁m₂m₃] = ((n−1)(n−2)/n²)ψ″(1) = −(4/9)ζ(3). Checked.
  * μ₃ ≡ 0 for n = 2. Checked.
  * The signs near edges (−) and vertices (+) are correct.
* **T8 heuristic.**
  * The Gamma shape depends on η, but this only multiplies the Borel integrand by an entire factor e^{(n−1)s/n}.
    So the singularity location is unaffected. Checked.
  * The edge equation e^{2s} = −1/(n−1) and the vertex value log 2 ± iπ are both correct.
  * |t_edge| = √(π² + log²2)/2 = 1.608575.

## 4. Literature

Searched 2026-10-21 (web): I found no expansion of the ℓ_p-sphere or superellipsoid surface area about p = 1 in any
dimension ≥ 3. That covers MathWorld "Supersphere", general ℓ_p surface-area and cone-measure work, and
Naor–Romik/Naor 2007. So the novelty statement "to our knowledge new" stands. Two remarks on that claim:
* Once T1 is in hand, T2/T3 are elementary consequences of standard Gamma/polygamma expansions.
* The residual risk that it appears in a paper whose full text I did not read is real but modest.

## 5. Numbers confirmed

| quantity | value |
|---|---|
| S₁ = [d¹]S | 4√3 |
| S₂ = [d²]S | 4√3(1 − π²/9) = −0.6694217800765659879982… |
| S₃ = [d³]S | 4√3(1 + π²/9 − 16ζ(3)/9) = −0.2796731272751286206150… |
| ζ(3) first appears at | d³ |
| coefficient of ζ(3) in S | −64√3/9 |
| S′(1) | S(1) |
| log form T5 | L₂ = −π²/9, L₃ = −16ζ(3)/9, L₄ = −37π⁴/810, L₅ = −160ζ(5)/27, L₆ = 16ζ(3)²/81 − 572π⁶/76545 |
| T6 | A₂, A₃, W₄ |
| n = 2 | diamond law |
| sign comparison | 3D −π²/9 and −16ζ(3)/9, against 2D +π²/24 and +ζ(3)/4 |

## 6. Corrections and caveats (none changes a result)

1. **Notation.** S_k means the Taylor coefficient [d^k]S(1+d), not the derivative. So S″(1) = 2S₂. "S₁ = 4√3" is the d¹
   coefficient, and it happens to equal S(1).
2. **"MZVs never appear".** This should read "no irreducible MZVs beyond products of single zeta values". The
   statement only has content from weight 8 (d⁸) on.
3. **Skewness channel.** The value −(4/9)ζ(3) and the sign facts are PROVED. The phrase "signed balance of edge
   regions against vertex regions" is interpretive: the integrand also changes sign in the interior, and the
   decomposition is not canonical (result.md concedes this in §5). Label that sentence as interpretation, not as
   part of the PROVED claim.
4. **T8.**
   * The edge-singularity location is NUMERICAL to about 2–3 digits. Higher-order extrapolation is unstable for
     n ≥ 3, and the exponent and amplitude are undetermined.
   * The statement "the vertex direction gives a farther singularity at log 2 ± iπ" is heuristic only. Nothing in
     the numerics resolves it.
   * Zero radius and Borel summability for n ≥ 3 remain CONJECTURES, as labelled.
5. **T3 write-up.** Add the two missing justifications from §3: cumulants of the underlying iid log G_j, and a
   dominated termwise expansion in Route B.
6. **PSLQ basis.** If Cl₂(π/3) is ever added to an enlarged basis, drop √3·L(χ₋₃,2): the two are rationally dependent.
