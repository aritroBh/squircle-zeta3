# Referee report on LAY1 `result.md` (Layer Theorem)

Referee: adversarial check, 2026-11-22. All code is my own and lives in `referee/`. I did not run or read the author's `lay1_*.py`. Interpreter: `Zeta(3)/.venv/bin/python3` (mpmath 1.3.0, sympy 1.14.0).

How each item was checked:
- **[P]** proof-reading by hand
- **[S]** symbolic computation (sympy)
- **[N]** numerics (mpmath)

| script | output | covers |
|---|---|---|
| `ref1_symbolic.py` | `ref1_symbolic_output.txt` | items 1 and 3 [S] |
| `ref2_hypotheses.py` | `ref2_hypotheses_output.txt` | items 2 and 3 [N] |
| `ref3_convexity_numeric.py` | `ref3_convexity_numeric_output.txt` | item 3 [N], dps 300 |
| `ref4_moments.py` | `ref4_moments_output.txt` | item 4 [N] |
| `ref_extra2.py` | `ref_extra2_output.txt` | item 4 [N], Lemma 5 through ₂F₁ |
| `ref5_endtoend.py` | `ref5_endtoend_output_dps60.txt`, `..._dps80.txt` | item 5 [N] |
| `ref_extra.py` | `ref_extra_output.txt` | 1/τ = 10 rows; τ = 0.01 convexity |

---

## Item 1. Proof of Theorem 1 (Section 2): **HOLDS** (three cosmetic fixes)

**[P] Lemma 1.**
- The majorant argument is correct: E_n is a polynomial in the y_k with nonnegative coefficients, the bound r^k|y_k| ≤ Y₀Ω^k holds, and e^X − 1 ≤ Xe^X.
- The tail bound needs τ/r ≤ 1, which is exactly τω ≤ Ω.
- The series in t is entire, so evaluating it at t = τ is legitimate.

**[P] Step 0.** The identity (2.1)–(2.2) is exact. A is bounded by (H1), so mA is integrable.

**[P] Step 1.** Correct. The bound is uniform in θ because Y₀ = max sup_T|a_k| and |exp(Σ τ^k a_k)| ≤ e^{Y₀N} for τ ≤ 1.

**[P] Step 2(a).** Correct.
- The Taylor expansion is in τ at fixed s, with the sign of s fixed.
- The N-th τ-derivative is bounded by (|H|+|S±T|)|s|^N ≤ (|H|+|S|+T)|s|^N.
- The bound |M_a(s)| ≤ (|H|+|S|+T)|s|^a/a! holds for **all** s, which Step 4 needs.

**[P] Step 2(b).** Correct.
- (H1) as stated gives what is used: log A = Σ_{k<N} τ^k a_k + O(τ^N), using sup|a_N| < ∞ and C₁τ^{N+1}. Taylor's theorem for a_k ∈ C^{N−k}(I_j) gives a remainder τ^k(τ|s|)^{N−k} = τ^N|s|^{N−k}.
- Lemma 1 applies with ω = 1+|s| and Ω = 1+W, because τω ≤ τ + W ≤ Ω when τ ≤ 1 and τ|s| ≤ W.
- The coefficient bound |P_n| ≤ Z(1+|s|)^n uses only the first claim of Lemma 1, so it also holds for all s.
- Y₀ here has to cover sup_{I_j}|a_k^{(i)}|; the text says so implicitly.

**[P] Step 2(c).** Correct. The s-dependent Y₀ = C₂e^{−c|s|} is applied pointwise, and |exp Σ| ≤ e^{C₂S₁}.

**[P] Steps 3–4 (range of each bound).**
- **Product-remainder identity:** r_m·A·E + m_<·r_A·E + m_<·A_<·r_E telescopes to mAE − m_<A_<E_<. Checked by expansion.
- **Range of the bounds:** every bound that needs |s| ≤ W/τ (r_m, r_A, r_E, and |A| ≤ C) is used only in (2.4), which is restricted to |s| ≤ W/τ. The tail estimate on |s| > W/τ uses only the coefficient bounds on M_a, P_b and G_c, which hold for all s. There is no out-of-range use.
- **Truncated terms with n ≥ N:** these use τ^n ≤ τ^N (τ ≤ 1) and a + b ≤ 2N − 2. Correct.

**[P] Step 5.** Correct. The outside region has measure ≤ 2π, and e^{−κ′/τ} ≤ (N+1)! τ^{N+1}/κ′^{N+1}.

**[S] Coefficients.** `ref1_symbolic.py` expands the generating functions with fully general symbols (H, S, T, a_k^{(i)}(θ_j), g₁, g₂, g₃), separately for s > 0 and s < 0. It verifies:
- B₂ and B₃ as stated;
- M₀ to M₃ as stated;
- P₂ = a₂ + a₁²/2 + a₁′s.

It then checks the claimed forms against the t^n coefficients of M·P·G:

| formula | claimed form | difference |
|---|---|---|
| Φ₀ | 0 | 0 |
| Φ₁ | HG₁ | 0 |
| Φ₂ | HG₂ + (M₁+HP₁)G₁ | 0 |
| Φ₃ (the c₄ formula) | HG₃ + (M₁+HP₁)G₂ + (M₂+M₁P₁+HP₂)G₁ | 0 |

So c₂, c₃ and c₄ in lines 53–55 are read off correctly.

**Fixes (cosmetic, none affects correctness):**
1. **Line 49.** For n = 0 the term Φ_{j,n−1} = Φ_{j,−1} is undefined. Add "Φ_{j,−1} := 0".
2. **Line 57 ("C depends only on …").** The list should also include:
   - J, the number of kinks;
   - sup_T|m|;
   - sup_T|a_k| for 1 ≤ k ≤ N (the bulk Y₀ in Step 1).
3. **Line 27 (non-decaying error remark).** Add the conclusion: (H2) then holds with c replaced by min(c, c′). This also needs κ − c′W > 0, which is what makes e^{−(κ−c′W)/τ} = O(τ^N).

## Item 2. Squircle hypotheses (Section 3): **HOLDS**

**[P] Algebra.** Each of the following is correct:
- h_p = ‖n‖_q, from Hölder duality with (q−1)p = q;
- h_p = m(1+r^q)^{1/q};
- the tent with H = T = 1/√2 and S = 0;
- r = tan(π/4 − x);
- log r = −2 artanh(tan x);
- Λ′ = sec²x/(1−tan²x) = sec 2x.

**[P] Bounds on Λ(x) − x.**
- sec y − 1 ≤ √2·y²/2 on [0, π/4], which gives Λ − x ≤ (2√2/3)x³.
- Monotonicity of (sec y − 1 − y²/2)/y⁴ comes from the nonnegative Taylor coefficients of sec. This gives K = 0.27802… < 0.28 and the remainder (32K/5)τ⁴|s|⁵.

**[P] N = 4 profile.** With F′ = −1/(1+e^y), the first-order term is τ·F′(2|s|)·(4/3)τ²|s|³ = −(4/3)τ³|s|³/(1+e^{2|s|}), so g₃ is correct. The remainder bound uses ξ ≥ 2|s| and |F″| ≤ e^{−ξ}.

**[P] Outside the windows.** By D₄ the outside set is |θ − jπ/2| ≤ π/8. There r ≤ √2 − 1 = e^{−log(1+√2)}, and Bernoulli's inequality applies. W = π/8 is strictly less than half the gap (π/4), as (H0) requires.

**[N] Checks in `ref2_hypotheses.py`, dps 120.** Every bound the text claims is ≤ 1 holds:

| quantity | max ratio | required |
|---|---|---|
| (Λ−x) / ((2√2/3)x³) | 0.8405 | ≤ 1 |
| (Λ−x−2x³/3) / ((16K/5)x⁵) | 0.9165 | ≤ 1 |
| \|L−τg₁\| / ((4√2/3)τ³\|s\|³e^{−2\|s\|}) | 0.7071 | ≤ 1 |
| \|L−τg₁−τ³g₃\| / (τ⁵(1+\|s\|)⁶e^{−2\|s\|}) | 0.693 | bounded |
| outside: (h−m) / (τ(√2−1)^{1/τ}) | 0.924 | ≤ 1 |

The grid covered τ from 10⁻¹ to 10⁻⁴. In addition:
- h_p = ‖n‖_q agrees with a brute-force sup over the boundary of the p-ball;
- log r + 2 artanh(tan x) = 10⁻¹²¹;
- Λ′ − sec 2x = 0.

## Item 3. Soft-max hypotheses and Proposition 4 (Section 4): **HOLDS WITH FIX** (one literal hypothesis mismatch; the convexity proof is correct)

**[P] Window formula and bounds.** All correct:
- the three exponents −2 sin u, −2 cos u and −2(cos u + sin u);
- b ≤ 2e^{−√2/τ};
- δ ≤ τ²|s|³/3 and δ ≤ 0.2056|s|, so the lower endpoint is ≥ 1.794|s|;
- the N = 4 expansion with x³/6 − x⁵/120 ≤ x − sin x;
- g₃ = |s|³/(3(1+e^{2|s|})) = (|s|³/6)(1 − tanh|s|).

**[N] Window formula.** The exact formula matches the raw vertex-sum definition to 10⁻¹²¹. The largest ratio to the claimed N = 3 bound is 0.7146, and the N = 4 remainder ratio is 0.025.

**[S] Proposition 4.** The derivative algebra is exact (c = cos, sn = sin, E_i, with D E_i = (ψ_i′/τ)E_i). The following all check:
- Φ + Φ″ = τfQ: difference 0.
- The displayed identity for τ²(1+ε)²Q: difference 0.
- The Cauchy–Schwarz step: RHS − (lower bound) = (1−τ)M[ε Σψ_i′²E_i − (Σψ_i′E_i)²] identically. The bracket is ≥ 0 by Cauchy–Schwarz with weights E_i, whose sum is ε.
- Φ′(0) = 0: verified by hand ([P]) with 1+ε = S and ε′ = −S/τ.

**Signs of the three terms on [0, π/4].**

| term | expression | sign claim | check |
|---|---|---|---|
| i = 2 | 2M′ψ′ − Mψ = −4c² + 6c·sn + 2sn² [S] | ≥ 4c²[M(1+τε) − τ(1+ε)] ≥ 4c²(1−τ) | correct: 6c·sn + 2sn² ≥ 0 and M ≥ 1 [P] |
| i = 3 | 2c² + 6c·sn − 4sn² [S] | > 0 | correct: 4(c − sn)sn ≥ 0 and 2Mc > 0 [P] |
| i = 4 | −2c² + 12c·sn − 2sn² = 2M² − 4M′² [S] | ≥ 2τ(1+ε)M² + 4M′²(1−τ) | correct [P] |

Q > 0 follows from the strictly positive i = 3 term. The reflection θ ↦ π/2 − θ and the C² gluing at the kinks are correct.

**[N] min over θ of h + h″, dps 300.** h is built from the raw vertex sum and h″ comes from `mpmath.diff`, so this is independent of the identity. Search: a grid of 401 points on [0, π/4] followed by golden-section refinement.

| τ | min(h + h″) | at θ |
|---|---|---|
| 1 | 0.541341 | 0 (the kink) |
| 0.5 | 0.867810 | π/4 |
| 0.2 | 0.0273855 | π/4 |
| 0.1 | 4.36913·10⁻⁵ | π/4 |
| 0.05 | 6.09525·10⁻¹¹ | π/4 |
| 0.02 | 5.60256·10⁻²⁹ | π/4 |

All values are positive, and the τ = 0.05 value matches the text's 6.1·10⁻¹¹. Extra check: at τ = 0.01 the true minimum is 2.17·10⁻⁵⁹, and at dps 55 it evaluates to 0. That is consistent with the text's diagnosis of AV5b's −2·10⁻⁵⁶ as rounding noise.

**Issue: line 201, "W = π/4. The four windows tile T".** (H0) (line 15) requires W *strictly smaller* than half the least gap. The gap is π/2, so W = π/4 is not allowed as written. Mathematically this is harmless: the tent formula holds for |u| < π/2, and the windows overlap only at endpoints, a set of measure zero. Either fix works:
- (a) In (H0), replace "smaller than half the least angular gap" with "at most half the least angular gap, with the tent formula required on all of I_j". Step 0 is unaffected because the overlaps are null.
- (b) Take W = π/4 − η and verify (H3) on η-neighbourhoods of π/4 + jπ/2. There 0 ≤ L ≤ τ(e^{−2 sin W/τ} + 2e^{−√2/τ}), so (H3) holds with κ = min(2 sin W, √2).

## Item 4. Lemma 5 and the moments (Section 5): **HOLDS** (one exposition fix)

**[P] Proof of the reflection identity.**
- The substitution u = e^{−2s} gives G = ½K(−λ/2, μ).
- The Beta split and the inversion u ↦ 1/u on (1, ∞) give K(−a−μ, μ) − 1/(a+μ), valid for Re(a+μ) < 0.
- The starting region {Re a > 0, Re(a+μ) < 0} meets the polydisc-type domain {|a| < 1, |a+μ| < 1} in a nonempty open set. That domain is convex, hence connected.
- The left side is analytic there, since K(a, μ) is analytic for Re a > −1.

**Gap in exposition: the pole cancellation is asserted but not shown.** It is subtle at (a, μ) = (0, 0), where both polar hyperplanes meet and 1/Γ(−μ) vanishes. **Fix: write it out.** Using Γ(z) = Γ(1+z)/z,

  Γ(a)Γ(−a−μ)/Γ(−μ) − 1/a + 1/(a+μ) = μ/(a(a+μ)) · [Γ(1+a)Γ(1−a−μ)/Γ(1−μ) − 1].

The bracket is analytic on the domain and vanishes on {a = 0} and on {a + μ = 0}. It is therefore divisible by a and then by a + μ (divide in a for fixed μ). So the right side is analytic on the whole domain, including (0, 0), and the identity theorem applies.

**[P] Wording at line 259.** The text compares left-hand coefficients 2(γ₀₂+γ₁₁) with "Σ ζ(k+1)μ^k" without saying that the right side carries the factor ½, so its μ^k coefficient is ζ(k+1)/2. The conclusions are correct; I checked them, and the μ¹ coefficient gives γ₀₁ = ζ(2)/4 as a consistency check. Say "the μ^k coefficient on the right is ζ(k+1)/2".

**[N] Lemma 5 at test points**, using direct s-quadrature for G and a second route G = (₂F₁(−μ, a; a+1; −1) − 1)/(2a) with a = −λ/2:

| (λ, μ) | \|LHS − RHS\| |
|---|---|
| (0.5, 0.2) | 5·10⁻⁵² |
| (−0.7, −0.3) | 0 |
| (0.25+0.1i, 0.4−0.2i) | 1·10⁻⁴⁹ (₂F₁ route: 7·10⁻⁶¹) |
| (1.5, 0.6) | 6·10⁻⁶¹ by ₂F₁ |
| (−3, −1.2), outside the stated domain | 5·10⁻⁴⁸ |

At (1.5, 0.6) the direct s-quadrature gives 1.5·10⁻¹⁷. That is quadrature error from the slowly decaying e^{−0.5s} integrand: the two quadratures of G(1.5, 0.6) differ from each other by 10⁻¹⁷. It is not a failure of the identity. The ₂F₁ route gives 6·10⁻⁶¹ there.

**[N] Moments by direct quadrature at dps 60.** All residuals are ≤ 8·10⁻⁶²:
- γ₀₁ = ζ(2)/4;
- γ₀₂ + γ₁₁ = ζ(3)/4;
- γ₂₁ = 7ζ(4)/64;
- J = γ₁₂ + γ₀₃ = ζ(4)/32;
- 2γ₀₃ + 2γ₁₂ + 4γ₂₁ = ζ(4)/2;
- ∫₀^∞ s³/(1+e^{2s}) ds = 21ζ(4)/64;
- the diagonal closed form G(τ, τ) = ¼[Γ(−τ/2)²/Γ(−τ) + 4/τ] at τ = 0.3 and τ = −0.4.

**[N + P] c₄ from the c₄ formula.** With A ≡ 1, S = 0 and T = H, the per-kink integrand is Φ₃ = H[g₃ + g³/6 + |s|g²/2 − s²g/2]. Integrating it directly, times 4 kinks:

| example | c₄ from quadrature | claimed | difference |
|---|---|---|---|
| squircle | −3.156937154220705612484910 | −(33√2/16)ζ(4) | 3·10⁻⁶¹ |
| soft-max | 0.270580808427784547879001 | ζ(4)/4 | 4·10⁻⁶² |

By hand, 8H[∫g₃ + J − γ₂₁] gives 4√2·(−28+2−7)/64·ζ(4) and 8·(7+2−7)/64·ζ(4). Correct.

## Item 5. End-to-end numerics: **HOLDS**

**Method.**
- I = c₀ + 8∫₀^{π/4}(h − m) dθ, with h − m = m·expm1(L) so that nothing cancels. The exact ∫m is added separately.
- Breakpoints sit at τ·{0, ¼, ½, 1, 1.5, 2, 3, 4, 6, …, 128} from the layer.
- Quadrature is tanh-sinh, cross-checked with Gauss–Legendre. The two agree to working precision: their difference is below one ulp at dps 60 and at dps 80, and quad's error estimates are ≤ 10⁻⁶².
- The dps 60 and dps 80 runs agree in all 50 printed digits of I.

So I(τ) is correct to at least 50 digits, beyond the ≥ 40 required.

**Results.**

| 1/τ | squircle R₃/τ⁴ | soft-max R₃/τ⁴ |
|---|---|---|
| 20 | −3.398529157816400091491 | 0.345898799375654791436 |
| 80 | −3.228470915453401717903 | 0.289222480556943123986 |
| 320 | −3.175460425951702591766 | 0.275237407198046956205 |
| 1280 | −3.161606709563698063477 | 0.271744839518098016656 |

I also ran 1/τ = 10, 40, 160, 640 and 2560. Every row of the author's tables in §6, including R₄/τ⁵, is reproduced digit for digit.

**Richardson extrapolation** (Neville, polynomial in τ, extrapolated to τ = 0). The difference from the claimed c₄:

| points used | squircle | soft-max |
|---|---|---|
| the 4 required (1/τ = 20, 80, 320, 1280) | 5.4·10⁻⁷ (6 digits) | 5.3·10⁻⁹ (8 digits) |
| 4 points, 1/τ = 320 … 2560 | 4.6·10⁻¹⁰ | 4.9·10⁻¹³ |
| 6 points, 1/τ = 80 … 2560 | 1.0·10⁻¹² | 2.2·10⁻¹⁵ |
| all 8 points, 1/τ = 20 … 2560 | 5.5·10⁻¹⁴ (**13 digits**) | 3.2·10⁻¹⁶ (**15 digits**) |

The steady gain with more points also supports the claim of a pure power series with no τ^k log τ terms.

Byproduct: c₅ ≈ −5.993407 for the squircle and ≈ 1.489923 for the soft-max, from Richardson on (R₃/τ⁴ − c₄)/τ. Both are consistent with the text's −5.99 and 1.4899.

---

## Status table (Section 7): statements I would downgrade

**None in substance.** Every PROVED entry survived the proof-reading, and every checkable formula survived the symbolic and numerical checks. I would add a qualifier to two rows, without downgrading them:

1. **"Soft-max: (H0)–(H3) for N = 3 and 4: PROVED".** This is literally inconsistent with (H0), because W = π/4 equals half the gap rather than being strictly less. Keep PROVED after either fix in Item 3, (a) or (b).
2. **"J = ζ(4)/32: PROVED (Lemma 5)".** This stays PROVED, but the pole-cancellation sentence in the proof of Lemma 5 should be replaced by the explicit formula in Item 4. As written, "the poles on the right cancel" is the one step a careful reader cannot check at (0, 0).

Minor edits that do not affect any status: Φ_{j,−1} := 0 (line 49); the full dependency list for C (line 57); "c replaced by min(c, c′)" (line 27); "the right side has coefficient ζ(k+1)/2" (line 259).
