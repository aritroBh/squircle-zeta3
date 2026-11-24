# Referee report on PUR1 (purity criterion for layer integrals)

Referee: independent Opus agent, 2026-11-22. Own code is in `referee/`: `r1_criterion.py` … `r5_scope.py`, each with its `_output.txt`. Tools: mpmath 1.3.0 and sympy 1.14.0 in a throwaway venv, deleted at the end. I read the author's scripts only after I had my own numbers, and I ran none of them. LAY1 §5 and §9 are taken as refereed.

**Verdict in one line.** The mathematics of (a), (b), (c) and (f) is right. (e) holds once one sentence of its justification is replaced. (g)(iii) contains a false statement: "τ⁴ is impure unless every H_j = 1". A rectangle soft-max with H ≠ 1 has a pure τ⁴, checked to 23 digits below.

There is no section called "surprises" in `result.md`. I read item 6's "surprises" as the Outcome bullets and the claims in §7, and I treat those below.

---

## 1. (a) Expansion with remainder: HOLDS WITH FIX (cosmetic)

**Proof-reading, line by line.**
- Step 1, the substitution v = τs: correct.
- Step 2: g ≤ min(log 2, e^{−2s}) gives B ≤ τg·2^τ ≤ 2τe^{−2s} for τ ≤ 1. The exponential-remainder bound is correct.
  - *Notation slip:* the sum should be Σ_{1≤b<L}, because B = e^{τg} − 1 has no b = 0 term.
- Step 3: the Taylor remainder of order K = N − 1 needs φ ∈ C^{N−1}[0, V₀], which is assumed. On [V₀, V), the bound |P_K| ≤ C(v/V₀)^K e^{βv} is fine, because deg P_K = K − 1 < K.
- Step 4: with βτ ≤ 1, ∫s^K e^{(βτ−2)s} ds ≤ K!, so the bound 2C_K K! τ^{N+1} is right. The constants depend on (N, φ, V₀, β) only, so the remainder is uniform in τ ∈ (0, min(1, 1/β)].
- Step 5:
  - The tail is ≤ 2τ∫_{V/τ}^∞ s^k e^{−2s} ds = O(τ^{−k} e^{−2V/τ}) ⊂ O(e^{−V/τ}). Correct.
  - With L = N − k, the b-sum runs to N − k − 1 and the error is O(τ^{N+1}) after the factor τ^{k+1}. Correct.
- **Interchange of the two expansions.** There is none to justify. Both the v-Taylor polynomial and the exponent series are *finite*, each with an explicit remainder, and the cross terms are bounded separately in steps 4 and 5. This is rigorous.
- The bound γ_kb ≤ 1/(b!(2b)^{k+1}) is correct.

**Fix (hypothesis).** "|φ| ≤ Ce^{βv} … automatic for V < ∞" is false if φ is integrable but unbounded near V. Two repairs work:
- assume φ bounded on [V₀, V); or
- keep φ ∈ L¹(V₀, V) and bound that piece directly: τ∫_{V₀/τ}^{V/τ}|φ(τs)|B ds ≤ 2τe^{−2V₀/τ}‖φ‖_{L¹}/τ, which is exponentially small.

None of the paper's examples is affected. The soft-max φ = 1 + v/√(1−v²) is bounded on [0, 1/√2].

**Symbolic / numerical.** Route (a), with γ_ab by quadrature, agrees with the criterion (B3) for the squircle F₂…F₇ to ≤ 9·10⁻⁴⁷ (`r1_criterion_output.txt`).

## 2. (b) The criterion: HOLDS

**Proof-reading.**
- **(B1) as formal power series.**
  - Write the left side as Σ_{i,j,b} ψ_i (i+j)!/j! · γ_{i+j,b} τ^{i+j+b}.
  - Write the right side, ψ(τ∂)G|_{λ=τ} = Σ_i ψ_i τ^i Σ_{a≥i,b} γ_ab a!/(a−i)! τ^{a−i+b}, then set a = i + j.
  - The two are the same sum, term by term. For each n the index set {i+j+b = n, b ≥ 1} is finite, so the identity is purely algebraic. No convergence of Taylor's theorem in λ is needed: the author's appeal to it is harmless but superfluous.
  - Each F_n uses only φ_k, k ≤ n − 2. This matches (a).
- **(B2) and the poles of S.**
  - (R) is an identity of functions analytic near (x, μ) = (0, 0), by LAY1 Lemma 5 with referee fix 2 of LAY1 §9.
  - The apparent poles of S lie on λ = 0 (x = −μ) and λ = 2μ (x = μ). Both lines pass through the origin, so S is analytic at the diagonal point only because they are removable.
  - Removability holds: Q(±μ, μ) = Γ(1−μ)Γ(1)/Γ(1−μ) = 1. A function analytic near (0,0) that vanishes on the two distinct lines x = ±μ is divisible by μ² − x² (divide by μ − x; the quotient vanishes on x = −μ off the origin, hence everywhere, by continuity). So (1−Q)/(μ²−x²) is a convergent power series at (0,0).
  - I verified this divisibility exactly, degree by degree, through total degree 8 (sympy `div`, zero remainder).
  - Differentiating (R) 2k times in x at x = 0 then gives 2∂_λ^{2k}G(μ,μ) = ∂_x^{2k}S|_{x=0} as power series in μ. So the τ-Taylor coefficients are what is claimed.
  - For odd derivatives the identity is 0 = 0, as the text says.
- **Weight and homogeneity.**
  - In log Q, Euler's γ cancels: (μ+x)/2 + (μ−x)/2 − μ = 0.
  - Every term of log Q is ζ(k) times a degree-k form, so it has weight 0. Hence Q has weight 0, and R := (1−Q)/(μ²−x²) has weight 2.
  - The τⁿ coefficient of τ^{2+2k}ψ_{2k}∂_x^{2k}R|_{x=0} comes from monomials of R of degree n − 2 in (x, τ), hence of ζ-weight n. The claim is correct.

**Symbolic.**
- (B3) built independently. Every coefficient through n = 7 is weight-homogeneous (automatic check).
- (B3) agrees with route (a) for the squircle, where ψ₂, ψ₄ and ψ₆ are all nonzero, so the even-derivative identity is tested beyond ψ ≡ 1.
- The special case ψ ≡ 1: the contour-extracted Taylor coefficients of the direct integral equal those of 1 − Γ(1−τ/2)²/Γ(1−τ) to ≤ 9·10⁻³⁸, for n = 0…9 (`r3_first_failure_output.txt`).

## 3. (c) Squircle: HOLDS

**Proof-reading.** I re-derived every step:
- tan(π/4−u) = e^{−2v} from tan u = tanh v;
- cos 2u = 1/cosh 2v;
- cos u + sin u = e^v/√cosh 2v;
- du/dv = 1/cosh 2v;
- hence φ = e^v/(√2 cosh^{3/2} 2v), ∫φ = 1/√2, and P = 4√2 + 8F.

D₄ gives 8 half-windows, and V = ∞, so there are no exponentially small terms.

**Symbolic.** ψ_k = √2/2, 0, −3√2/2, 0, 13√2/4, 0, −379√2/60, as stated. My (B3) gives:

| n | c_n from (B3) | equals claimed |
|---|---|---|
| 2 | √2ζ(2) | yes |
| 3 | √2ζ(3) | yes |
| 4 | √2(ζ(4) − 7ζ(2)²)/8 = −(33√2/16)ζ(4) | yes, after ζ(2)² = (5/2)ζ(4) |
| 5 | −√2(¾ζ(5) + (7/4)ζ(2)ζ(3)) | yes |
| 6 | √2(4039ζ(6)/256 − 7ζ(3)²/8) = 20.91147814643259204894245555535839107867 | yes, after ζ(2k) → π^{2k} |
| 7 | √2(129ζ(7)/16 + 141ζ(2)ζ(5)/16 + 1071ζ(3)ζ(4)/64) = 63.54451807765424743136703150673447368701 | yes |

The c₆ and c₇ rows differ from the text's forms symbolically only by ζ(2k) relations, and the difference is 0 exactly.

**Numerical** (`r2_squircle_quad_output.txt`), independent of the change of variables:
- direct θ-quadrature of 8∫₀^{π/4} cos θ(1+tan^{1/τ}θ)^τ dθ;
- the v-form equals θ-quadrature at τ = 0.5, 0.3, 0.1 and 0.03, to ≤ 4·10⁻⁴⁸ at dps 50;
- P(½) = 2π.
- Neville extrapolation of θ-quadrature alone (dps 90, τ = 1/(20·2^i), i < 14) gives c₆ − exact = 4.7·10⁻³³ and c₇ − exact = −2.5·10⁻³².

*Trivial inconsistency.* The Outcome says "37 digits" and §3 says "≥ 36". Make them agree.

## 4. (d) Soft-max (spot check): HOLDS

- ψ = 1, 0, −½, 5/6, −5/8, 79/120 (sympy).
- c₂…c₇ from (B1), with the even part by (B3) and the odd part by numerical γ's, equal route (a) to ≤ 1.4·10⁻⁴⁵.
- c₆ = 40γ₃₂ − ζ(3)²/2 + (11/16)ζ(6) = 0.042950827880284033929433804835135266 agrees with route (a) to 7·10⁻⁴⁷.
- γ₃₁ = 15ζ(5)/256 and γ₄₁ = η(6)/32 are confirmed.

## 5. (e) First failure: HOLDS WITH FIX

**Proof-reading.**
- The bookkeeping of Proposition E is right. Odd k contributes ψ_k[k!γ_{k1}τ^{k+2} + ((k+1)!γ_{k+1,1} + k!γ_{k2})τ^{k+3} + …], and odd k ≥ 2j+3 starts at τ^{2j+5}.
- The middle term ω(2j+2)!η(2j+4)/2^{2j+3} is genuinely needed. Without it the formula is off by 0.2368, 0.7392 and 5.604 per unit ω for j = 0, 1, 2.
- Implicit hypothesis: purity through τ^{2j+3} needs ψ_even and ω rational, as in Theorem B. State it.

**FIX (interchange in the Euler-sum lemma).** The sentence "the series converges absolutely, as Σ H_{n−1}/n < ∞" is false: Σ H_{n−1}/n diverges, and absolute convergence fails at x = 1, that is, s = 0. Replace it with Tonelli:

  Σ_n (2H_{n−1}/n) ∫₀^∞ s^a e^{−2ns} ds = Σ_n 2H_{n−1}·a!/(n(2n)^{a+1}) < ∞ for every a ≥ 0,

so the termwise integration is justified. The formula itself is right.

**Numerical** (`r3_first_failure_output.txt`).
- *Direct test of Prop E.* Take φ = e^v(1 + ωv^{2j+1}), V = ∞, for j = 0, 1, 2 and ω = 1 and 3/7.
  - F was evaluated by quadrature of its defining integral in the s-form. This is the analytic continuation, analytic for Re τ < 2.
  - The integral was taken on the circle |τ| = 0.75 with 96 nodes, and the Taylor coefficients extracted by discrete Cauchy integral. No γ's were used on this side.
  - |F_{2j+4}(contour) − Prop E| = 5·10⁻³⁷, 2·10⁻³³ and 9·10⁻³⁰ for j = 0, 1, 2. That is ≥ 29 digits.
  - Lower orders: the coefficients of v^{2j+1}e^v vanish below τ^{2j+3}, and equal (2j+1)!η(2j+3)/2^{2j+2} at τ^{2j+3}, to the same accuracy.
- *Euler sums.* γ_{a,2} as an Euler sum (Cohen–Villegas–Zagier acceleration at dps 110) agrees with quadrature for a = 0…6 to ≤ 1.5·10⁻⁶³, which is the quadrature precision.
- *Classical weight-4 value.* Σ(−1)^{n+1}H_n/n³ − [(11/4)ζ(4) − (7/4)ζ(3)log 2 + ½ζ(2)log²2 − log⁴2/12 − 2Li₄(½)] = 8·10⁻¹¹². So the γ₁₂ closed form holds to 2·10⁻¹¹².
- *PSLQ, dps 100.*
  - γ₁₂ against the alternating basis returns exactly the stated relation.
  - γ_{1,2}, γ_{3,2} and γ_{5,2} against the pure bases give no relation.
  - γ₂₂ = (16ζ(2)ζ(3) − 29ζ(5))/256 and γ₄₂ = (−251ζ(7) + 64ζ(2)ζ(5) + 112ζ(3)ζ(4))/4096 are found, as stated.
- *The text's test ψ = (1+v)e^{−v²}.* Route (a) by hand gives F₄ = γ₀₃ + 2γ₁₂ + γ₂₁ = 9ζ(4)/64 + γ₁₂ = γ₁₂ + π⁴/640. This agrees with the text and with Prop E: ψ₂ = −1 contributes −2γ₂₁, which cancels the middle term 7ζ(4)/32.

**Is non-purity claimed only conditionally?** Mostly yes: the status table, the Outcome and §5 say CONDITIONAL and NOT KNOWN correctly. Three places overstate it:
1. §5, "PROVED (the sums are not *formally* pure): … exactly the antisymmetric data that (R) leaves free." This is a heuristic count, not a theorem. That (R) does not determine a number says nothing about its arithmetic. Reword it as a remark.
2. §7(iii), "the τ⁶ failure is unavoidable" for H = 1. This is unconditional in tone. It holds only conditionally on γ₃₂ ∉ ℚ[ζ], and modulo the unrefereed polygon window reduction.
3. §10, item 4, "γ₁₂ contains Li₄(½)". This is true of the displayed expression. It is an irreducibility statement only under the same independence conjecture.

## 6. (f) Designed square: HOLDS (closed form, series and convexity all proved)

**Closed form, numerical** (`r4_designed_output.txt`, dps 60). Full-circle quadrature of h(θ) as defined, with u the signed angle from the nearest jπ/2, no symmetry reduction:

| τ | ∫₀^{2π} h | closed form − quadrature |
|---|---|---|
| 1 | 12 | 0 |
| 0.7 | 10.04606976218073610135388875514462992446 | −5·10⁻⁵⁷ |
| 0.37 | 8.568140896358803953502010031659757790131 | 1·10⁻⁶⁰ |
| 0.1 | 8.035430423800770449801944314832356216326 | 1·10⁻⁶⁰ |
| 0.04 | 8.00542082853631133638892314262757217917 | 0 |

That is ≥ 56 digits at five values of τ.

**Derivation of the exponentially small series.**
- (e^Λ)′ = cos u + sin u = M. So M du = e^v dv, with v = Λ ∈ [0, log 2], and Λ(π/4) = log(2 + √2/2 − √2/2) = log 2.
- ∫h = 8 + 8F_V, and F_V = τG(τ,τ) − T(V).
- T(V) by the binomial series: for v ≥ V the series converges uniformly, since x = e^{−2v/τ} ≤ 4^{−1/τ} < 1. ∫_V^∞ e^{v(1−2n/τ)} dv = τe^{V}e^{−2nV/τ}/(2n−τ) for τ < 2.
- So 8T(log 2) = 16τΣC(τ,n)4^{−n/τ}/(2n−τ). Correct. Numerically, at τ = 0.6, series − quadrature = 7·10⁻⁴⁸.
- At τ = 1, 1/Γ(0) = 0 and only n = 1 survives: 16 − 4 = 12. Correct.

**Series through τ⁷.** mpmath Taylor coefficients of 16 − 8Γ(1−τ/2)²/Γ(1−τ) equal every claimed coefficient to ≤ 2·10⁻⁶¹, including 79π⁶/120960 = (237/384)ζ(6).

**Convexity: an actual inequality chain, valid for every τ ∈ (0,1]; no grid is used.** I re-derived the identity

  h + h″ = 2σfMe^{−Λ}[−3M′ + M²e^{−Λ}(1+2K)]

from M″ = −M, f′ = −2fσΛ′, f″ = fσ(4Λ′²K − 2Λ″), Λ′ = Me^{−Λ} and Λ″ = M′e^{−Λ} − M²e^{−2Λ}. Then:
1. K = (1−σ)/τ + σ ≥ 1 exactly when τ ≤ 1.
2. So the bracket is ≥ 3e^{−Λ}(M² − M′e^Λ).
3. Since e^Λ = 2 − M′ and M² + M′² = 2, we get M² − M′e^Λ = 2(1 − M′) = 2(1 − cos u + sin u). This is > 0 on (0, π/4].
4. The prefactor is > 0.

Hence h + h″ > 0 on (0, π/4) for all τ ∈ (0,1].

The identity also agrees with sympy's own second derivative of h to relative 10⁻¹⁸⁴ … 10⁻²¹⁰ at dps 210.

- **u = 0.**
  - h′(0⁺) = 2^τ(1 − Λ′(0)) = 0, since Λ′(0) = M(0)e^0 = 1. Numerically ≤ 5·10⁻⁹⁹.
  - h is even about 0, so h ∈ C² there and there is no kink.
  - (h+h″)(0⁺) = 2^τ(1/τ − 1) ≥ 0, matched numerically. It vanishes only at τ = 1.
- **u = π/4.**
  - M′(π/4) = 0 and f′ < 0, so h′(π/4⁻) < 0.
  - The reflection θ ↦ π/2 − θ gives jump = −2h′(π/4⁻) > 0 for every τ ∈ (0,1]. A positive jump is a positive atom of h + h″, that is, a short edge, not a concave corner.
  - Exact jump formula: −2h′(π/4⁻) = 2√2·2fσΛ′(π/4) = **4(1 + 4^{−1/τ})^τ/(1 + 4^{1/τ})**. Numerically it is 1.0, 0.24253563, 3.9032009·10⁻³, 3.814694·10⁻⁶, 3.6379788·10⁻¹², 3.1554436·10⁻³⁰ at τ = 1, ½, 0.2, 0.1, 0.05, 0.02, equal to −2h′(π/4⁻) computed from the sympy derivative.

**min(h + h″) on (0, π/4], dps 210, 801-point grid**, with h + h″ computed by sympy differentiation:

| τ | min | at u |
|---|---|---|
| 1 | 1.2·10⁻⁵⁹ (→ 0) | 0⁺ |
| 0.5 | 0.418658 | π/4 |
| 0.2 | 0.0151691 | π/4 |
| 0.1 | 2.83226·10⁻⁵ | π/4 |
| 0.05 | 5.2735·10⁻¹¹ | π/4 |
| 0.02 | 1.12677·10⁻²⁸ | π/4 |

All grid values are positive. These are illustrations only; the proof above does not depend on them.

Outside the claimed range the proof's step K ≥ 1 fails, and so does convexity: (h+h″)(0⁺) = −0.383 at τ = 1.2. So "(0,1]" is sharp at the edge normals.

The support-function step is standard: h is piecewise C², with h + h″ ≥ 0 on each piece and positive jumps of h′, so h + h″ ≥ 0 as a measure (Schneider §1.7 / §2.5). It would help to cite that statement explicitly, not just "§1.7".

## 7. (g) Scope

**(i) Linear images of B_p: HOLDS WITH FIX.**
- I re-derived the computation: x dy − y dx = det A dθ, dt = det A dθ/x², |x| = Q(1,t)^{−1/2}. The other side uses Q(σa, b) = Q(a, σb), since Q is even. So φ = e^v[q(v)+q(−v)], with each half-window v ∈ [0, ∞). This is correct and exact.
- **My test: A = [[1,1],[0,2]]**, a shear that is not symmetric under any reflection of the axes. Here Q = 5ξ₁²/4 − ξ₁ξ₂/2 + ξ₂²/4, and the kinks are at θ = 0, 3π/4, π, 7π/4.
  - The kink at x = y has ψ = 2 + 6v² − 28v⁴ + …. The kink at x = −y has ψ = √2(½ − (3/16)v² − (181/256)v⁴ + …). Both are even, as predicted.
  - From (B3):
    - c₂ = (4+√2)ζ(2)/4;
    - c₃ = (4+√2)ζ(3)/4;
    - c₄ = 5ζ(2)²/8 + 13ζ(4)/8 + √2(25ζ(4) − 7ζ(2)²)/128;
    - c₅ = 5ζ(2)ζ(3)/4 + 9ζ(5)/4 + √2(9ζ(5) − 7ζ(2)ζ(3))/64;
    - c₆ is weight-6 homogeneous.
  - Numerics (`r5_scope_output.txt`, dps 60–70):
    - the exact v-form equals θ-quadrature at τ = 0.5, 0.2 and 0.07, to ≤ 1.3·10⁻⁶⁰;
    - c₀ = 4 + 4√2, the perimeter of A·B₁;
    - Neville extrapolation of θ-quadrature alone (τ = 1/(20·2^i), i < 11) confirms c₄, c₅ and c₆ to 1·10⁻²², 2·10⁻²¹ and 3·10⁻²¹.
    - So this linear image is pure through τ⁶, as Theorem B predicts.
  - *Quadrature caveat for anyone repeating this:* h is only C^{⌊1/τ⌋} where x = 0 or y = 0, so those angles must be quadrature breakpoints. Without them I first saw a spurious 5·10⁻²¹ at τ = 0.2.
- **FIX.** "Coefficients in ℚ(entries)·ζ" is false as stated: this integer matrix produces √2. The per-kink prefactor is Q(1,±1)^{−3/2}/det A. The correct statement is: coefficients in ℚ(entries of A, √Q(1,1), √Q(1,−1)), each times a weight-n ζ-polynomial. The diag(1,2) example's √5 is the same phenomenon.
- Unproved but harmless: "the evenness is about v₀ = ½log(a/b)".

**(ii) Designed tents: HOLDS WITH FIX.**
- Checked:
  - (e^Λ)′ = (T/H²)M, so M du = (H²/T)e^v dv;
  - h′(0±) = ±2^τ(T − HΛ′(0)) with Λ′(0) = T/H;
  - the bracket generalises to (Me^{−Λ}/c)[−3M′ + M²e^{−Λ}(1+2K)/c], with c = H²/T;
  - M² − ce^ΛM′ vanishes at 0, and its derivative M(M′ + ce^Λ) = M(c+T) is > 0, since M′ + ce^Λ ≡ c + T identically;
  - the bracket at 0 equals T(1/τ − 1).
- **Gap.** The positivity of the h′-jump at the split angle of the rectangle is only numerical in the text. It has a one-line proof, which should be inserted:
  - the split angle exists and is unique, because Λ₁(θ) increases from 0 and Λ₂(θ) decreases to 0 on [0, π/2] (both Λ′ = Me^{−Λ}/c with M > 0), so the intermediate value theorem applies;
  - at the split, m and f agree from both sides;
  - so jump = m(f₂′ − f₁′) = 2mfσ(Λ₁′ − Λ₂′) > 0, because Λ₁′ > 0 > Λ₂′ (derivatives in θ).
  - My numerics for A = 1, B = 2: split θ* = 0.2660800472, jump 1.875, 0.1946, 2.19·10⁻⁴ at τ = 0.9, 0.3, 0.1.

**(iii) Polygon soft-max: FAILS as stated.**
- The local computation is right: φ = 2H/T + 2tan u with v = T sin u, and ψ₁ = 2(1−H)/T per kink.
- But ψ₁ enters F linearly, so the γ₁₂ coefficient of c₄ is **Σ_j 2(1−H_j)/T_j**, summed over the kinks, and this can vanish with H_j ≠ 1. Then the whole k = 1 odd tower Σ_jψ_{1,j}·τ²∂G(τ,τ) cancels at every order.
- **Counterexample.** The rectangle [−6/5, 6/5] × [−3/5, 3/5]:
  - (1 − 6/5)/(3/5)·2 + (1 − 3/5)/(6/5)·2 = −2/3 + 2/3 = 0;
  - the ψ-coefficients of the two kink types are ψ₁ = −2/3 and +2/3;
  - prediction: c₂ = (5/2)ζ(2), c₃ = (5/2)ζ(3), **c₄ = (5/16)ζ(4)**, which is pure.
  - Direct θ-quadrature of h = m[Σ_v e^{(v·n−m)/τ}]^τ (dps 60, τ = 1/(100·2^i), i < 9), with Neville extrapolation, gives c₄ = 0.3382260105347306848487907 against (5/16)ζ(4) = 0.3382260105347306848487512, a difference of 4·10⁻²³. A single kink pair's γ₁₂ term would be 0.029, so the test is sensitive. c₅ also matches its prediction to 8·10⁻²³.
- **Correct statement:** "pure through τ³ always; the τ⁴ coefficient contains γ₁₂ with coefficient Σ_j 2(1−H_j)/T_j."
- "τ⁶ failure unavoidable" is right for H_j ≡ 1, where every kink has ψ₃ > 0. In general its γ₃₂ coefficient is Σ_jψ_{3,j}, and the claim is conditional on γ₃₂ ∉ ℚ[ζ].
- The general-polygon window reduction remains a sketch, as the author says.
- One caution: this "soft-max" m(Σe^{(v·n−m)/τ})^τ is not translation-covariant, so the H_j depend on the origin. This is consistent with the formula above.

**Other unproved assertions in (g) and in the Outcome.**
- "It gives an independent complete proof … of thm:pw" is right for the planar squircle perimeter only. The n-dimensional extension is explicitly not examined.
- The Deligne / period-conjecture route to non-purity is cited, not carried out. It is labelled CONDITIONAL, which is fine.
- The reducibility of odd-weight alternating linear Euler sums is cited from memory. It is NUMERICAL here (PSLQ), which is fine.

---

## Statements in the status table I would downgrade

1. **(g) "polygon soft-max: τ⁴ impure unless all H_j = 1".** Change it to FALSE AS STATED. The γ₁₂ coefficient is Σ_j 2(1−H_j)/T_j; the rectangle 6/5 × 3/5 has a pure c₄ = (5/16)ζ(4), checked to 23 digits. "τ⁶ always" becomes CONDITIONAL, and only for H_j ≡ 1.
2. **(g) "linear images of B_p … coefficients in ℚ(entries)".** Keep PROVED, but correct the field to ℚ(entries, √Q(1,±1)): A = [[1,1],[0,2]] gives √2.
3. **(g) "designed tents … convex".** It is PROVED only after the one-line jump-at-split argument is inserted. As written, that step is NUMERICAL.
4. **(e) "γ_{a,2} = Euler sum: PROVED".** Keep PROVED, with the Tonelli justification replacing the false absolute-convergence claim.
5. **(e) body text "PROVED (the sums are not formally pure)".** Downgrade to a REMARK (heuristic).
6. **(a) hypothesis "automatic for V < ∞".** Add the boundedness or L¹ clause. No status change.

Everything else in the table I would keep as it stands:
- (a), (b), (c) all rows;
- (d) under SM6's reduction;
- (e) Proposition E;
- (f) closed form, series and convexity for every τ ∈ (0,1].
