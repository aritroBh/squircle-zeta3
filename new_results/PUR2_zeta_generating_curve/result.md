# PUR2: a convex curve whose perimeter generates every ζ(k) with coefficient 2

**Outcome.** The planner's claims (A)–(D) hold. They need these corrections:

- (A) holds for 0 < τ < 2/(1+|c|) as stated, and for every τ > 0 by analytic continuation (poles cancel).
- (D): ψ''''(0) does not enter τ⁴ or τ⁵. It first enters at τ⁶.
- (C): the formula holds for every κ. Admissibility (φ > 0 on [0,V]) is the condition κ < κ_max = 1.5078108109….

(E) is proved in a stronger form than asked:

- For every real c, K_{τ,c} is a convex body for every τ ∈ (0,1]. This is a three-line inequality, with no grid.
- The same holds for the τ⁴-free curve κ² = 9/7. That proof uses an interval-arithmetic certificate for one numerical constant.
- For no member of the class, and no τ > 1, is K_{τ,c} convex.

One extra result came out (D′): within the class, the perimeter series equals 8 + 2Σζ(k)τ^k to order N if and only if ψ agrees with cosh v to order N − 2 at 0.

Scripts (mpmath 1.3, sympy; project venv):

| script | output | what it does |
|---|---|---|
| `pur2_symbolic.py` | `pur2_symbolic_output.txt` | (A) coefficients in ab, (B) series, (C) τ⁴ zero, (D) universal coefficients, squircle cross-check |
| `pur2_numerics.py` (+ `pur2_numerics_lib.py`, the same class without the driver) | `pur2_numerics_output.txt` | (F) u-quadrature of h vs closed form − tail; h+h″ formula; grid minima; jumps |
| `pur2_convexity.py` | `pur2_convexity_output.txt` | κ_max, interval certificate at κ² = 9/7, κ_s, κ_D, real-c D minima |
| `pur2_kappa_scan.py` | `pur2_kappa_scan_output.txt` | κ ∈ (κ_D, κ_max): proved τ₀ and numerical τ* |
| `pur2_checks.py` | `pur2_checks_output.txt` | Taylor coefficients of the closed forms, the elementary (B) integral, the Beta lemma at complex α |

---

## 0. Setting (from PUR1 §0, §6)

**The class.** Fix φ = e^vψ with:
- ψ smooth on [0, V];
- ψ(0) = 1;
- the Taylor series of ψ at 0 even;
- φ > 0 on [0, V].

Let Φ(v) = ∫₀^vφ, and let V be defined by Φ(V) = 1. On 0 ≤ u ≤ π/4 put:
- M = cos u + sin u and M′ = cos u − sin u;
- p(u) = 1 + sin u − cos u = ∫₀^u M, which runs over [0, 1];
- Λ(u) = Φ^{−1}(p(u)).

Λ is smooth and increasing, by the inverse function theorem, since Φ′ = φ > 0. Also Λ(0) = 0, Λ(π/4) = V and Λ′ = M/φ(Λ).

**The body.** In each window |u| ≤ π/4 around jπ/2, set

  h = m·f,  m = cos u + |sin u|,  f = (1 + e^{−2Λ(|u|)/τ})^τ.

h is symmetric under u ↦ −u and under reflection in π/4. Substituting v = Λ(u) gives M du = φ dv, so

  ∫₀^{2π} h = 8∫₀^V φ(v)(1+e^{−2v/τ})^τ dv = 8 + 8F,  F = ∫₀^V φ B_τ(v) dv,  B_τ(v) = (1+e^{−2v/τ})^τ − 1.

**The family.** φ_c = e^v cosh(cv) = ½(e^{2av} + e^{2bv}), with a = (1+c)/2, b = (1−c)/2, a + b = 1 and ab = (1−c²)/4.
- c real: ψ = cosh(cv).
- c = iκ: ψ = cos(κv), ab = (1+κ²)/4.
- c = 0 is PUR1's designed square.

---

## 1. A Beta lemma (self-contained; it is PUR1's (R) at λ = 2α)

**Lemma 1.** Let α, β ∈ ℂ with Re α < 1 and Re β < 1, and put τ := α + β. Then

  ∫₀^∞ (e^{2αs} + e^{2βs})[(1+e^{−2s})^τ − 1] ds = ½[Γ(−α)Γ(−β)/Γ(−τ) + 1/α + 1/β].

The right side is continued analytically across α = 0 and across β = 0.

*Proof.*
1. **Left side is analytic.** For |x| ≤ 1, |(1+x)^τ − 1| ≤ |τ|2^{|τ|}|x|. So the left side converges locally uniformly on the convex, hence connected, domain D = {Re α < 1, Re(τ−α) < 1} ⊂ ℂ², and is analytic there.
2. **Right side is analytic on D.**
   - 1/Γ(−τ) is entire.
   - Inside D, the only poles of Γ(−α) and Γ(−β) are α = 0 and β = 0.
   - At α = 0, Γ(−α) = −1/α + O(1) and Γ(−β) → Γ(−τ). So the α-pole cancels against 1/α. The same holds at β = 0.
3. **The identity on an open subset of D.** On the set U = {Re τ < Re α < 0} ⊂ D (nonempty and open), substitute x = e^{−2s}. Each term splits as
   ½∫₀¹ x^{−α−1}(1+x)^τ dx + 1/(2α).
   Map x ↦ 1/x in the β-term and use β − τ = −α. This turns ∫₀¹ x^{−β−1}(1+x)^τ dx into ∫₁^∞ x^{−α−1}(1+x)^τ dx. Hence the left side equals
   ½∫₀^∞ x^{−α−1}(1+x)^{τ} dx + 1/(2α) + 1/(2β) = ½[B(−α, α−τ) + 1/α + 1/β].
   This uses the Beta integral ∫₀^∞x^{p−1}(1+x)^{−p−q}dx = B(p,q) with Re p, Re q > 0. Here p = −α and q = α − τ = −β.
4. **Conclusion.** By the identity theorem the two sides agree on D. ∎

*Check:* α = 0.3 + 0.7i, τ = 0.6 agree to 5·10⁻⁴¹ (quadrature-limited; `pur2_checks_output.txt`).

**Rewriting.** Take α = aτ and β = bτ, with a + b = 1. Use Γ(−z) = −Γ(1−z)/z and 1/Γ(−τ) = −τ/Γ(1−τ). The bracket becomes (1 − R)/(abτ), where

  R(τ) := Γ(1−aτ)Γ(1−bτ)/Γ(1−τ).

---

## 2. (A) The family: PROVED

**Theorem A.** Let c be real with a = (1+|c|)/2, or c = iκ. Suppose 0 < τ < 1/Re a. That is, τ < 2/(1+|c|) for real c, and τ < 2 for imaginary c. Then

  **perimeter(K_{τ,c}) = 8 + (2/(ab))[1 − Γ(1−aτ)Γ(1−bτ)/Γ(1−τ)] − E_c(τ),**

  **E_c(τ) = 2τ Σ_{n≥1} C(τ,n) e^{−2nV/τ} [e^{2aV}/(n−aτ) + e^{2bV}/(n−bτ)].**

Here C(τ,n) is the binomial coefficient. The leading term of the tail is E_c = 2τ²e^{−2V/τ}[e^{2aV}/(1−aτ) + e^{2bV}/(1−bτ)](1 + O(e^{−2V/τ})).

For c = iκ the right side is real, since a and b are conjugate. The statement holds for every κ ≥ 0. The *curve* exists only for κ < κ_max (§4).

At c = 0 the theorem is PUR1 §6: 16 − 8Γ(1−τ/2)²/Γ(1−τ) − 16τΣC(τ,n)4^{−n/τ}/(2n−τ).

*Proof.*
1. **Split off the tail.** F = F_∞ − T, where F_∞ = ∫₀^∞φB_τ and T = ∫_V^∞φB_τ. Both converge: |φ| ≤ e^{2(Re a)v} and 0 ≤ B_τ ≤ τ2^τe^{−2v/τ}.
2. **The main part.** Substitute v = τs:
   F_∞ = (τ/2)∫₀^∞(e^{2aτs} + e^{2bτs})[(1+e^{−2s})^τ − 1]ds.
   By Lemma 1 with α = aτ and β = bτ (Re α, Re β < 1 by hypothesis), F_∞ = (τ/2)·½(1−R)/(abτ) = (1−R)/(4ab).
3. **The tail.** For v ≥ V we have x = e^{−2v/τ} ≤ e^{−2V/τ} < 1, so B_τ = Σ_{n≥1}C(τ,n)x^n converges absolutely, with |C(τ,n)| bounded polynomially in n.
   - Σ_n|C(τ,n)|∫_V^∞|φ|e^{−2nv/τ}dv < ∞, because n/τ > Re a. So termwise integration is legitimate (Fubini).
   - Each term is ½C(τ,n)∫_V^∞(e^{2av}+e^{2bv})e^{−2nv/τ}dv = ½C(τ,n)τ e^{−2nV/τ}[e^{2aV}/(2n−2aτ) + e^{2bV}/(2n−2bτ)].
   - Hence 8T = E_c.
4. **Assemble.** Perimeter = 8 + 8F_∞ − 8T. ∎

**Remark A′ (all τ > 0; PROVED).**
- The left side, 8∫₀^Vφ exp(τ log(1+e^{−2v/τ}))dv, is real-analytic in τ > 0.
- On the right, the Gamma term and E_c are meromorphic in τ near (0, ∞). E_c converges there, since e^{−2nV/τ} decays geometrically in n. Their poles are at aτ ∈ ℕ or bτ ∈ ℕ.
- The identity on (0, 1/Re a) therefore continues to all τ > 0 off those points, and the poles cancel.

Example, c = 1 at τ = 1. The direct value is 8 + 8∫₀^V φe^{−2v}dv = 8 + 4V + 2(1 − e^{−2V}) = 12.07161687483714381923862516…, equal to the u-quadrature to 60 digits. At τ = 1 − 10⁻²⁵ the closed form minus tail differs from it by 7.5·10⁻²⁵.

**Coefficients.** For |τ| < min(1, 1/|a|, 1/|b|) use log Γ(1−z) = γz + Σ_{k≥2}ζ(k)z^k/k. Euler's γ cancels because a + b − 1 = 0:

  **log R = Σ_{k≥2} (a^k + b^k − 1) ζ(k) τ^k / k.**

**Lemma (polynomiality).** a^k + b^k − 1 is a polynomial in e := ab, divisible by e.
- By Newton's identity p_k = p_{k−1} − e·p_{k−2}, with p₀ = 2 and p₁ = 1, so p_k ∈ ℤ[e]. Also p_k(0) = 1, since a + b = 1 and e = 0 force {a, b} = {1, 0}.
- Explicitly:
  - k = 2: −2e
  - k = 3: −3e
  - k = 4: −4e + 2e²
  - k = 5: −5e + 5e²
  - k = 6: −6e + 9e² − 2e³
  - k = 7: −7e(1−e)²

Hence log R ∈ e·ℚ[e][ζ(2), ζ(3), …][[τ]], so 1 − R is divisible by e. Every perimeter coefficient is therefore a polynomial in ab and the ζ(k), homogeneous of weight n at τⁿ.

Through τ⁷ (sympy, `pur2_symbolic_output.txt`, before ζ(2k) relations):

| n | coefficient of τⁿ in 8 + (2/e)(1−R) |
|---|---|
| 2 | 2ζ(2) |
| 3 | 2ζ(3) |
| 4 | 2ζ(4) − e(ζ(4) + ζ(2)²) = **(2 − 7e/2) ζ(4)** |
| 5 | 2ζ(5) − 2e(ζ(5) + ζ(2)ζ(3)) |
| 6 | (2 − 13e/2 + 31e²/8) ζ(6) − e ζ(3)² |
| 7 | 2ζ(7) − e(2ζ(2)ζ(5) + 2ζ(3)ζ(4) + 4ζ(7)) + e²(ζ(2)²ζ(3) + 2ζ(2)ζ(5) + ζ(3)ζ(4) + 2ζ(7)) |

At e = 1/4 (c = 0) these reproduce PUR1 §6. Taylor coefficients of the Gamma form at c = ½ and at κ² = 9/7 match this table to 3·10⁻⁵¹ and 6·10⁻⁵¹ (`pur2_checks_output.txt`).

---

## 3. (B) The limit c = 1: PROVED directly

**Data.**
- φ = (e^{2v} + 1)/2 and Φ = (e^{2v} − 1)/4 + v/2.
- Φ(Λ) = p(u) gives e^{2Λ} + 2Λ = 1 + 4p. Hence

  **Λ(u) = ½[q − W(e^q)],  q = 5 + 4 sin u − 4 cos u,  V = ½[5 − W(e⁵)] = 0.65327932051967509783813032715…**

  Here W is the principal Lambert W: W(e^q) solves E + log E = q, and E = e^{2Λ}.

**Theorem B.** For 0 < τ < 1,

  **perimeter(K_{τ,1}) = 8 − 2τ(γ + ψ(1−τ)) − 2τ Σ_{n≥1} C(τ,n) e^{−2nV/τ}[W(e⁵)/(n−τ) + 1/n]**
  **= 8 + 2 Σ_{k≥2} ζ(k) τ^k − (same tail),**

where e^{2V} = W(e⁵) = 3.69344135896…. At τ = 1 the perimeter is 8 + 4V + 2(1 − e^{−2V}), which is the limit (Remark A′). The leading tail is 2τ²e^{−2V/τ}[W(e⁵)/(1−τ) + 1](1 + O(e^{−2V/τ})).

*Direct proof (no limit c → 1).*
1. **Tail.** It is §2 step 3 with a = 1 and b = 0.
2. **Reduction.** For the main part, substitute v = τs and x = e^{−2s}:
   8F_∞ = 2τ·I(τ),  I(τ) := ∫₀¹(x^{−τ−1} + x^{−1})[(1+x)^τ − 1]dx.
3. **Map x ↦ 1/x in the first piece.** This gives ∫₁^∞[(1+y)^τ − y^τ]/y dy. It converges because (1+y)^τ − y^τ = O(y^{τ−1}) and τ < 1.
   For Y > 1, ∫₁^Y(y^τ − 1)/y dy = (Y^τ − 1)/τ − log Y. So

   I = lim_{Y→∞}[J(Y) − (Y^τ−1)/τ + log Y],  J(Y) = ∫₀^Y((1+x)^τ−1)/x dx.

4. **Evaluate J(Y).** Put 1 + x = 1/(1−t), with T = Y/(1+Y). Then dx/x = dt/(t(1−t)), and 1/(t(1−t)) = 1/t + 1/(1−t). So

   J(Y) = ∫₀^T((1−t)^{−τ} − 1)/t dt + ((1+Y)^τ − 1)/τ − log(1+Y).

   Hence I = lim[∫₀^T((1−t)^{−τ}−1)/t dt + ((1+Y)^τ − Y^τ)/τ − log(1 + 1/Y)] = ∫₀¹((1−t)^{−τ} − 1)/t dt.
5. **Series.** Map t ↦ 1 − t and expand 1/(1−t) = Σtⁿ. Monotone convergence applies, since t^{n−τ} − tⁿ ≥ 0:

   I = Σ_{n≥0}[1/(n+1−τ) − 1/(n+1)] = −γ − ψ(1−τ) = Σ_{k≥1}ζ(k+1)τ^k  (|τ| < 1).

   The last two equalities are the standard series for ψ and its expansion in τ, with absolute convergence. ∎

*Checks.*
- I(τ) + γ + ψ(1−τ) = −3·10⁻⁵² at τ = 0.3 and 8·10⁻⁵¹ at τ = 0.77.
- sympy confirms −2τ(γ + ψ(1−τ)) = 2Σζ(k)τ^k through τ⁷.
- The general formula at e = 0 gives 2ζ(n) at every n ≤ 8.

---

## 4. (C) Imaginary c = iκ and the τ⁴-free curve

**Admissibility: PROVED (κ_max NUMERICAL to 20 digits).**
- φ = e^v cos κv and Φ(v) = [e^v(cos κv + κ sin κv) − 1]/(1+κ²).
- Φ increases on [0, π/(2κ)] and Φ(π/(2κ)) = (κe^{π/(2κ)} − 1)/(1+κ²).
- **The curve exists, with φ > 0 on [0,V], iff Φ(π/(2κ)) > 1, i.e. κe^{π/(2κ)} > 2 + κ², i.e. 0 ≤ κ < κ_max = 1.50781081090595335…** (κ_max² = 2.2734934414848686…).

*Why the condition is an interval.* q(κ) = log κ + π/(2κ) − log(2+κ²) is strictly decreasing on (0, ∞). Its derivative has the sign of −2κ³ − πκ² + 4κ − 2π, which is < 0 because 4κ − πκ² ≤ 4/π < 2π.

*At κ = κ_max* φ(V) = 0, so Λ′(π/4) = ∞. *Beyond κ_max* Φ turns down before reaching 1, and Λ does not exist.

**Formula: PROVED.** Theorem A holds for all 0 < τ < 2, since Re a = ½. ab = (1+κ²)/4 ranges over [¼, 0.8184) on the admissible range.

**The τ⁴ zero: PROVED.** The τ⁴ coefficient (2 − 7ab/2)ζ(4) vanishes iff ab = 4/7, i.e. κ² = 9/7 (κ = 1.1338934…). This is admissible, since κ²/κ_max² = 0.566.

The τ⁴-free curve has

  **perimeter = 8 + 2ζ(2)τ² + 2ζ(3)τ³ + 0·τ⁴ + (6/7 ζ(5) − 8/7 ζ(2)ζ(3))τ⁵ − (22/49 ζ(6) + 4/7 ζ(3)²)τ⁶ + O(τ⁷) − E(τ),**

with V = 0.776780305889861281…, φ(V) = 1.38414….

- The τ⁵ coefficient is −1.37098118164544991…; the τ⁶ coefficient is −1.28244672897550285….
- The Taylor coefficient of the Gamma form at τ⁴ is 0 to 50 digits.

**Range of the τ⁴ coefficient** over the admissible family:
- real c gives (9 + 7c²)ζ(4)/8 ≥ (9/8)ζ(4);
- imaginary c gives (9 − 7κ²)ζ(4)/8, down to −0.86ζ(4) as κ → κ_max.

So the sign of c₄ can be chosen.

---

## 5. (D) Universality in the class: PROVED

**Theorem D.** For any φ = e^vψ in the class of §0 (ψ even, smooth, ψ(0) = 1), as τ → 0,

  perimeter = 8 + Σ_{n≥2} c_n τⁿ + O(τ^{N+1}) for every N,

where

  c₂ = 2ζ(2) = π²/3,  c₃ = 2ζ(3),
  **c₄ = (9 + 7ψ″(0)) ζ(4)/8**  [= (7/4)ζ(4) − ζ(2)²/4 + ψ″(0)(ζ(2)² + ζ(4))/4],
  **c₅ = (3 + ψ″(0)) ζ(5)/2 + (ψ″(0) − 1) ζ(2)ζ(3)/2,**
  c₆ = (79/128 + (73/64)ψ″(0) + (31/128)ψ⁗(0)) ζ(6) + (ψ″(0) − 1) ζ(3)²/4,
  c₇ = (9/8 + ¾ψ″ + ⅛ψ⁗)ζ(7) + (−3/8 + ¼ψ″ + ⅛ψ⁗)ζ(2)ζ(5) + (−9/32 + ψ″/16 + 7ψ⁗/32)ζ(3)ζ(4),

with every derivative taken at 0.

**Correction to the claim.** ψ⁗(0) does not enter c₄ or c₅. In general c_n depends only on ψ^{(2k)}(0) with 2k ≤ n − 2, so ψ^{(2k)}(0) first appears at τ^{2k+2}.

*Proof.*
1. **Hypotheses.** φ is bounded and smooth on [0, V], as PUR1 (a) requires after the referee's fix. So PUR1 Prop. A and Theorem B apply, and perimeter = 8 + 8F with an O(τ^{N+1}) remainder.
2. **Expansion.** Since ψ is even, (B3) gives F ~ τ²Σ_k ψ^{(2k)}(0)τ^{2k}[x^{2k}]𝓡(x,τ). Here 𝓡 = (1−Q)/(τ²−x²) and log Q = Σζ(k)/k[((τ+x)/2)^k + ((τ−x)/2)^k − τ^k].
3. **Degree count.** 𝓡 is a power series whose degree-d part has weight d + 2 (PUR1 §0 and the referee's divisibility argument). So the τⁿ coefficient uses only [x^{2k}τ^{n−2−2k}]𝓡.
4. **Computation.** The table above is a sympy computation: the division by τ² − x² is exact degree by degree (`pur2_symbolic.py`), and the result is reduced by ζ(2)² = (5/2)ζ(4), ζ(2)ζ(4) = (7/4)ζ(6) and ζ(2)³ = (35/8)ζ(6). ∎

*Cross-checks.*
- ψ = cosh cv (ψ^{(2k)}(0) = c^{2k}) reproduces §2's table exactly, through τ⁸.
- The PUR1 squircle: ψ = (√2/2)/cosh^{3/2}(2v) has ψ″ = −6, ψ⁗ = 156 and ψ⁽⁶⁾ = −9096 after normalising. It reproduces PUR1's c₄ … c₇ (difference 0 / 10⁻¹⁷⁰).

**Corollary D′ (rigidity; PROVED).**
- The coefficient of ψ^{(2k)}(0) in c_{2k+2} is 8[x^{2k}]𝓡(x,0) = 8[x^{2k+2}]Q(x,0). Here Q(x,0) = Γ(1−x/2)Γ(1+x/2) = (πx/2)/sin(πx/2) = Σ_n 2η(2n)(x/2)^{2n}, so the coefficient is **4^{1−k}η(2k+2) > 0** (k = 1: 7ζ(4)/8; k = 2: 31ζ(6)/128, matching the table).
- Hence, by induction on k: **c_n = 2ζ(n) for all n ≤ N if and only if ψ^{(2k)}(0) = 1 for all 2k ≤ N − 2.** That is, ψ has the Taylor jet of cosh v.
- So among the class, c = 1 is, up to flat perturbations, the only curve whose perimeter series is 8 + 2Σζ(k)τ^k.

---

## 6. (E) Convexity

**The formula (PROVED).** On 0 < u < π/4 put y = 2Λ/τ, σ = 1/(1+e^y) ∈ (0, ½] and K = σ + (1−σ)/τ. From f′ = −2fσΛ′, f″ = 4fσKΛ′² − 2fσΛ″, M″ = −M, Λ′ = M/φ and Λ″ = M′/φ − M²φ′/φ³,

  **h + h″ = 2fσ(M/φ)·B,  B = −3M′ + (M²/φ)(2K + φ′/φ),**

with φ evaluated at Λ(u). At c = 0 this is PUR1's formula. It matches mp.diff to relative error ≤ 3·10⁻⁶⁰ for c = 1, ½, 0 and κ² = 9/7.

Split off the τ-dependence: 2K = 2 + 2(1−σ)(1/τ − 1). Then

  **B = D + (M²/φ)·2(1−σ)(1/τ − 1),  D := −3M′ + (M²/φ)(2 + φ′/φ),**

where D does not depend on τ.

**Key identity (PROVED).** Use p = Φ(Λ) ∈ [0,1] as the variable. Then w := M′ = 1 − p and M² = 2 − w². Put χ(p) := ψ′/ψ at v = v(p). Since dφ/dp = φ′/φ = 1 + χ, we have φ = 1 + p + ∫₀^pχ, and

  **Dφ = 6p + (2 − w²)χ(p) − 3w∫₀^p χ(p′)dp′.**

(The constant part: 3(2 − w²) − 3w(2 − w) = 6p.) This was checked numerically at u = 0.5 for c = ½ and κ² = 9/7, to 20 digits.

**Boundary behaviour (PROVED, whole class).**
- **u = 0.** h′(0⁺) = f(0)[M′(0) − 2σ(0)Λ′(0)] = 2^τ(1 − 1/φ(0)) = 0. h is even, so h ∈ C² across the edge normals. Also B(0) = 1/τ − 1.
- **u = π/4.** M′ = 0, so h′(π/4⁻) = √2f′ = −4fσ/φ(V) < 0. By the reflection in π/4,

  **jump of h′ at π/4 = −2h′(π/4⁻) = 8 f σ / φ(V) = 8(1+e^{−2V/τ})^τ/[(1+e^{2V/τ})φ(V)] > 0.**

  At c = 0 this is the referee's 4(1+4^{−1/τ})^τ/(1+4^{1/τ}). It is a positive atom of h + h″, that is, a short edge of length 8fσ/φ(V), not a concave corner. It agrees with −2h′(π/4⁻) from mp.diff to 12 printed digits.

**Support-function step.** h is continuous and piecewise C², is C¹ at the edge normals, and has positive h′-jumps at the diagonals. Then h + h″ ≥ 0 as a measure, provided it is ≥ 0 on each open piece, and h is the support function of a convex body (Schneider, *Convex Bodies*, §1.7 and §2.5, cited from memory). Its perimeter is ∫h by Cauchy's formula.

### E1. Real c, and every ψ with log ψ convex: convex for every τ ∈ (0,1] (PROVED)

**Theorem E1.** If log ψ is convex on [0, V], then K_{τ} is a convex body for every 0 < τ ≤ 1, and h + h″ > 0 on (0, π/4). This holds for ψ = cosh(cv) with any real c, so it covers c = 0, ½, 1 and all c > 1.

*Proof.*
1. ψ is even, so ψ′(0) = 0. log ψ convex then makes χ = ψ′/ψ nondecreasing and ≥ 0 in v, and hence in p.
2. So ∫₀^pχ ≤ pχ(p). Since w ≥ 0,

   Dφ ≥ 6p + χ(p)[(2 − w²) − 3wp] = 6p + χ(p)(1 − p + 2p²) ≥ 6p > 0  on (0, 1].

   (1 − p + 2p² has negative discriminant.)
3. For τ ≤ 1 the second term of B is ≥ 0, so B ≥ D > 0 on (0, π/4].
4. With the boundary facts above, h + h″ is a nonnegative measure. ∎

(This replaces the case analysis I first had, which needed |c| < 6.) Grid minima of Dφ/M² are positive for c = 0, ½, 1, 2, 3, 6, 10, 30 and 100.

### E2. Imaginary c = iκ (PROVED where stated)

Here χ = −κ tan(κv) ≤ 0 and χ is decreasing.

**(i) Convex for all τ ∈ (0,1] if κ ≤ κ_s, with κ_s = 1.18271 (κ_s² = 1.39880).**
1. The integral term of Dφ is ≥ 0, so Dφ ≥ 6p − (1 + 2p − p²)κ tan(κv(p)) ≥ 6p − 2κ tan(κv(p)).
2. Let k(p) = κ tan(κv(p)). Then k(0) = 0 and k′(p) = κ²sec²(κv)/φ(v).
3. So k′ ≤ 3, and hence k ≤ 3p and D ≥ 0, whenever κ² ≤ 3e^v cos³(κv) on [0, V].
4. e^v cos³κv is unimodal on [0, π/(2κ)], because its derivative has the sign of 1 − 3κ tan κv. So it suffices that **κ² ≤ 3·min(1, e^V cos³(κV)).**

**At κ² = 9/7 this is certified by interval arithmetic** (`pur2_convexity_output.txt`):
- Φ(V ∓ 10⁻²⁰) brackets 1;
- on that V-interval, 3e^V cos³(κV) ∈ [1.6825197373002186654…, 1.6825197373002186656…] > 9/7.

So **the τ⁴-free curve is a convex body for every τ ∈ (0,1]: PROVED.** κ_s is the root of the sufficient condition (NUMERICAL).

**(ii) Every admissible κ, for small τ.** Use 2(1−σ) ≥ 1 and −Dφ/M² ≤ −χ(p) − 6p/(2−w²) ≤ κ tan(κV). Then B ≥ 0 whenever

  **τ ≤ τ₀(κ) := 1/(1 + κ tan(κV))  (PROVED; V numerical).**

More generally, for any ψ in the class, the body is convex for τ ≤ 1/(1 + δ), with δ = max(0, sup(−Dφ/M²)) < ∞. δ is finite by continuity of D on [0, π/4] and φ > 0. So **every member of the class is convex for all sufficiently small τ, with an explicit τ₀.**

**(iii) Numerical sharp picture** (grid of 300 points, `pur2_convexity_output.txt`, `pur2_kappa_scan_output.txt`):
- D ≥ 0, and hence convexity for all τ ≤ 1, holds up to κ_D = 1.35542 (κ_D² = 1.83717).
- For κ_D < κ < κ_max, D < 0 on an interval ending at π/4, and convexity fails for τ near 1 (at τ = 1, B = D). It holds up to a numerical τ*(κ):

| κ | κ² | φ(V) | min Dφ/M² | τ₀ proved | τ* numerical |
|---|---|---|---|---|---|
| 1.3 | 1.69 | 1.094 | +0.011 | 0.296 | 1 |
| 1.36 | 1.85 | 0.947 | −0.066 | 0.246 | 0.963 |
| 1.4 | 1.96 | 0.825 | −0.80 | 0.208 | 0.698 |
| 1.45 | 2.10 | 0.624 | −2.57 | 0.152 | 0.434 |
| 1.48 | 2.19 | 0.445 | −5.41 | 0.106 | 0.270 |
| 1.5 | 2.25 | 0.243 | −13.5 | 0.057 | 0.129 |
| 1.507 | 2.271 | 0.080 | −49.1 | 0.019 | 0.039 |

τ* → 0 as κ → κ_max, where φ(V) → 0 and φ′/φ → −∞ at π/4. Wherever convexity fails, ∫h is only the signed length of the hedgehog, not a perimeter.

### E3. τ > 1: never convex (PROVED)

For every member of the class, B(0) = 1/τ − 1 < 0, and the prefactor at 0 is 2fσM/φ = 2^τ. So (h + h″)(0) = 2^τ(1/τ − 1) < 0. h is C² there, so h + h″ < 0 on a neighbourhood of each edge normal, and h is not a support function.

Numerically, (h + h″)(10⁻⁶) = −0.38289 at τ = 1.2 for all four curves, against 2^{1.2}(1/1.2 − 1) = −0.38289. The formula of Theorem A still holds for τ > 1 (Remark A′), but only as a signed length.

**Summary of (E).**
- Real c (including c = 0, ½, 1): convex for all 0 < τ ≤ 1. PROVED.
- κ ≤ κ_s, including κ² = 9/7: convex for all 0 < τ ≤ 1. PROVED, with an interval certificate.
- κ < κ_max: convex for τ ≤ 1/(1 + κ tan κV). PROVED.
- Sharp κ-threshold κ_D and the thresholds τ*(κ): NUMERICAL.
- τ > 1: never convex. PROVED.

---

## 7. (F) Numerics (own code, dps 60)

Perimeter = 8∫₀^{π/4} h du. This is direct tanh-sinh quadrature in u, with Λ computed at each node:
- c = 1: by Lambert W;
- c = ½: by Cardano, since z = e^{v/2} solves z³ + 3z = 4 + 3p;
- κ² = 9/7: by Newton.

The geometric breakpoints are π/4·2^{−k}. Each result is compared with Theorem A/B: closed form minus E.

| curve | τ | ∫₀^{2π} h (45 digits) | closed − E − ∫h | E |
|---|---|---|---|---|
| c = 1 | 0.5 | 9.08083797485323736848623271077471697789365114 | 1.2·10⁻⁶⁰ | 0.30546 |
| c = 1 | 0.2 | 8.15446356642325109885142123516278617123205254 | 0 | 6.5·10⁻⁴ |
| c = 1 | 0.1 | 8.03554204092969452808047118366066193485885701 | 0 | 2.2·10⁻⁷ |
| c = 1 | 0.05 | 8.00853939513026973543145328292188661040949078 | −1.2·10⁻⁶⁰ | 1.1·10⁻¹³ |
| c = 1 | 1 | 12.0716168748371438192386251629641173275005311 | 0 (vs 8+4V+2(1−e^{−2V})) | — |
| c = ½ | 0.5 | 9.06303713755488051896767763605563383980218886 | 0 | 0.19670 |
| c = ½ | 0.2 | 8.15309087667533640027170001307298351057175022 | 1.2·10⁻⁶⁰ | 4.2·10⁻⁴ |
| c = ½ | 0.1 | 8.0354582759934502776897133738987979089007784 | 0 | 1.1·10⁻⁷ |
| c = ½ | 0.05 | 8.00853457997715587461304994280427265214541039 | 0 | 3.1·10⁻¹⁴ |
| c = ½ | 1 | 12.0190894355127902046797118180477444956475808 | −2.5·10⁻⁶⁰ | — |
| κ² = 9/7 | 0.5 | 9.02004451543642351432486129685425748377291271 | 0 | 0.038848 |
| κ² = 9/7 | 0.2 | 8.15021505033253291083355832602676979972141777 | −1.2·10⁻⁶⁰ | 8.7·10⁻⁵ |
| κ² = 9/7 | 0.1 | 8.0352877499955753527152882476406415624873461 | 0 | 9.6·10⁻⁹ |
| κ² = 9/7 | 0.05 | 8.0085247357384665959865369587426053437869009 | 2.5·10⁻⁶⁰ | 4.4·10⁻¹⁶ |
| κ² = 9/7 | 1 | 11.8830165085179060081404356149009550656536909 | 1.2·10⁻⁶⁰ | — |

So the agreement is ≥ 59 digits in all 15 cases. The tails are non-negligible at τ ≥ 0.2, so the test sees them.

**min of h + h″ on a 400-point grid of (0, π/4]:**

| τ | c = 1 | c = ½ | κ² = 9/7 | c = 0 (PUR1, for comparison) |
|---|---|---|---|---|
| 1 | 0.0274 (at u = π/1600, → 0 as u → 0) | 0.0244 | 0.0184 | 0.0234 |
| 0.5 | 0.395 (at π/4) | 0.412 | 0.457 | 0.419 |
| 0.2 | 0.0173 | 0.0157 | 0.0120 | 0.0152 |
| 0.1 | 4.69·10⁻⁵ | 3.26·10⁻⁵ | 1.04·10⁻⁵ | 2.83·10⁻⁵ |
| 0.05 | 1.91·10⁻¹⁰ | 7.57·10⁻¹¹ | 3.75·10⁻¹² | 5.27·10⁻¹¹ |

For τ ≤ 0.5 the c = 0 column agrees with PUR1's referee table to all printed digits. At τ = 1 the referee's grid included u = 10⁻⁶⁰, where h + h″ → 0; mine starts at π/1600.

**Further measured values.**
- Jumps of h′ at π/4, at τ = 0.5, 0.2, 0.1: 0.2412, 4.95·10⁻³, 7.22·10⁻⁶ (c = 1); 0.2529, 2.44·10⁻³, 1.03·10⁻⁶ (κ² = 9/7).
- h′(0⁺) ≈ 6·10⁻⁶⁵ (rounding).

---

## 8. Literature (modest; references from memory, not reopened)

Nothing here is new as analysis.
- Σ_{k≥2}ζ(k)τ^k = −τ(γ + ψ(1−τ)) is the Taylor series of the digamma function (DLMF §5.7).
- I(τ) = ∫₀¹(t^{−τ} − 1)/(1−t)dt is Gauss's integral for ψ (DLMF §5.9).
- log Γ(1−z) = γz + Σζ(k)z^k/k is classical (DLMF §5.7).
- R = Γ(1−aτ)Γ(1−bτ)/Γ(1−τ) = (1−τ)B(1−aτ, 1−bτ) is a Beta value, and Lemma 1 is a regularised Beta integral (DLMF §5.12).

The point of PUR2 is only the geometric realisation: an explicit convex curve, convexity proved, whose perimeter equals these classical generating functions up to an explicit exponentially small series. I make no novelty claim for that either; I did not search the literature.

## 9. Status

| item | status |
|---|---|
| Lemma 1 (Beta lemma, complex α, β, Re < 1) | **PROVED** (+ check 5·10⁻⁴¹) |
| (A) perimeter = 8 + (2/ab)(1−R) − E_c, 0 < τ < 2/(1+\|c\|) (real c), 0 < τ < 2 (c = iκ) | **PROVED** + NUMERICAL (≥ 59 digits) |
| (A′) same identity for all τ > 0 by continuation, poles cancel | **PROVED** (+ τ = 1 checks) |
| (A) log R = Σ(a^k+b^k−1)ζ(k)τ^k/k; coefficients ∈ ℚ[ab][ζ], weight n | **PROVED**; table through τ⁷ symbolic |
| (B) c = 1: Λ via Lambert W; perimeter = 8 − 2τ(γ+ψ(1−τ)) − E, 0 < τ < 1; = 8 + 4V + 2(1−e^{−2V}) at τ = 1 | **PROVED directly** + NUMERICAL |
| (C) admissible iff κ < κ_max = 1.507810810906… | **PROVED** (κ_max numerical root) |
| (C) τ⁴ coefficient (2 − 7ab/2)ζ(4), zero at κ² = 9/7; τ⁵, τ⁶ values | **PROVED** |
| (D) c₂ = 2ζ(2), c₃ = 2ζ(3) universal; c₄, c₅ depend only on ψ″(0); ψ⁗(0) first at τ⁶ | **PROVED** (claim corrected) |
| (D′) series = 8 + 2Σζ(k)τ^k to order N ⟺ ψ has cosh's jet to order N−2 | **PROVED** |
| (E) h ∈ C² at edge normals; positive h′-jump 8fσ/φ(V) at diagonals (whole class) | **PROVED** |
| (E1) log ψ convex (all real c): convex for every τ ∈ (0,1] | **PROVED** |
| (E2) κ² = 9/7: convex for every τ ∈ (0,1] | **PROVED** (interval certificate for one constant) |
| (E2) any admissible κ: convex for τ ≤ 1/(1+κ tan κV); any ψ: some explicit τ₀ > 0 | **PROVED** |
| (E2) sharp κ_D = 1.3554, τ*(κ) table | **NUMERICAL** |
| (E3) τ > 1: never convex, whole class | **PROVED** |
| (F) 15 quadrature checks, grid minima | **NUMERICAL** (≥ 59 digits) |
| novelty | **not claimed** |

## 10. Five-sentence summary

1. For φ_c = e^v cosh(cv), the convex smoothing K_{τ,c} of the square, with h = m(1+e^{−2Λ_c/τ})^τ, has perimeter 8 + (2/ab)[1 − Γ(1−aτ)Γ(1−bτ)/Γ(1−τ)] minus an explicit series in e^{−2nV/τ}, where a = (1+c)/2 and b = (1−c)/2. The proof is a regularised Beta integral, and every τⁿ coefficient is a weight-n polynomial in zeta values and ab.
2. At c = 1, with Λ = ½[q − W(e^q)] and q = 5 + 4 sin u − 4 cos u, the perimeter is 8 − 2τ(γ + ψ(1−τ)) − tail = 8 + 2Σ_{k≥2}ζ(k)τ^k − tail. This is proved by an elementary computation that reduces to Gauss's integral for ψ.
3. Within the whole class φ = e^vψ (ψ even, ψ(0) = 1), the τ² and τ³ coefficients are always 2ζ(2) and 2ζ(3). The τ⁴ and τ⁵ coefficients depend only on ψ″(0), and ψ^{(2k)}(0) enters first at τ^{2k+2} with weight 4^{1−k}η(2k+2) > 0, so the series 8 + 2Σζ(k)τ^k characterises cosh's Taylor jet.
4. Taking c = iκ with κ² = 9/7 kills the τ⁴ term exactly. This curve is admissible (κ < κ_max ≈ 1.5078), and it is convex for every τ ∈ (0,1] by an inequality plus a one-constant interval certificate.
5. Every curve with log ψ convex (all real c) is a convex body for every τ ∈ (0,1]; every member is convex for small τ with an explicit τ₀; none is convex for τ > 1. Fifteen direct u-quadratures match the closed forms to ≥ 59 digits.

## Planner notes after the referee (2026-11-22)

Referee: `../PUR2R_referee_zeta_curve/REFEREE.md` (independent, own code). Mathematics correct; no claim false. Repairs to read into the text above:
1. Lemma 1, step 2: the point α = β = 0 must be treated. Corrected: "The right side equals ½(τ/(αβ))(1 − R), with R = Γ(1−α)Γ(1−β)/Γ(1−α−β). Since R(0,β) = R(α,0) = 1, 1 − R = αβ·g with g analytic, so the right side is ½τg, analytic on all of D."
2. At c = ±1, ab = 0 and (2/ab)(1 − R) means its limit, analytic in c; Theorem B gives it.
3. Convexity criterion: use the referee's self-contained curve argument (it also shows ∫h is the perimeter, short edges included) instead of the citation from memory. Two windows meet at each diagonal, not four.
4. The "min Dφ/M²" column of the κ table and the real-c grid minima in §6 are grid artefacts (D(0) = 0, so the true minimum is 0); only the sign near π/4 means anything.
Numerical only: κ_s ≈ 1.18271, κ_D ≈ 1.35542, the τ*(κ) table. κ_max certified in (1.50781081090595335, 1.50781081090595336) by the referee (arb).
