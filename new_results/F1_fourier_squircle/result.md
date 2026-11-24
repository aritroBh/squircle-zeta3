# F1: Fourier meets the squircle

Folder: `research/10_frontier/F1_fourier_squircle/`. Every number below is printed by a script in this folder; the
`.out` files hold the runs. Python: the project `.venv` (mpmath, sympy). Started fresh (no partial results from
the earlier killed run were found, either here or in `<scratch>/`).

Notation: squircle |x|^p+|y|^p = 1, ε = 1/p, dual exponent q = p/(p−1), **τ = 1/q = 1 − 1/p** (τ = 0 is the diamond
p = 1, τ = ½ the circle, τ = 1 the square), λ = log 2, L = log(1+√2), G = Catalan. Every function below has the symmetry
of the square, so f = a₀/2 + Σ_{k≥1} a_k cos 4kθ with a_k = (4/π)∫₀^{π/2} f cos 4kθ dθ.

---

## 0. Summary (one paragraph)

The right Fourier object for the squircle is the **support function** h_p(θ) = (|cos θ|^q+|sin θ|^q)^{1/q}
(= 1/ρ_q, the reciprocal radial function of the dual squircle), because Cauchy's formula makes the perimeter a single
Fourier coefficient, **P(p) = π·a₀(h_p)**. With it:

1. **(Diamond end, main theorem, PROVED.)** For every k ≥ 0,
   **a_k(p) = (−1)^k (√2/π) [ −4/(16k²−1) + ζ(2)τ² + ζ(3)τ³ − (33/16 + 7k²)ζ(4)τ⁴ + O_k(τ⁵) ].**
   Every Fourier coefficient of the squircle's support function picks up the same ζ(2) at second order and the same
   ζ(3) at third order. k = 0 gives **P/√2 = 4 + ζ(2)τ² + ζ(3)τ³ − (33/16)ζ(4)τ⁴ − (¾ζ(5)+⁷⁄₄ζ(2)ζ(3))τ⁵ + …**,
   which reproduces P″(1) = √2π²/3 and P‴(1) = 6√2(ζ(3) − π²/3). Equivalently, the rounding of the diamond deposits
   Dirac masses (√2/4)(ζ(2)τ² + ζ(3)τ³) at the four edge normals.
2. **(Pure weight, PROVED.)** In τ = 1 − 1/p the n-th coefficient has pure weight n (alternating MZVs).
   D11's d-expansion mixes weights only because d = τ/(1−τ). Checked on D11's exact table for n ≤ 10.
3. **(Circle end, PROVED.)** The velocity of the homotopy at the circle,
   u = ∂_p log ρ_p|₂ = ∂_p h_p|₂ = −½(cos²θ log|cos θ| + sin²θ log|sin θ|), has Fourier coefficients
   **−1/(4k(4k²−1))**. That is an *odd* rational function of k, which no polynomial cosine series can have.
   Pairing u with Aritro's own polynomials gives every odd zeta value:
   **ζ(3) = (4/π)∫₀^{π/2} u(θ)[(4θ−π)² + 8 − π²/3] dθ**, and ζ(2j+1) = 8⟨u, β_{2j} − 4β_{2j−2}⟩.
4. **(Circle end, Parseval, PROVED + NUMERICAL.)** Hurwitz's isoperimetric Parseval sum Σ(16k²−1)a_k(h_p)² behaves
   differently at second and third order:
   - its ε² term (first-order data only) is 7π²/48 − 1. That is even weight, as the parity lock predicts.
   - its ε³ term is **Σ(16k²−1)α_kβ_k = −(15/32)ζ(3) + 43π²/384 − (7/48)π²λ + λ − 5/16**. Here the second-order
     coefficients β_k carry harmonic numbers. It is exact by Hurwitz + P, A, and checked to 10⁻²¹.

Why this works where the Fourier method failed: ∂_p|x|^p = |x|^p log|x|. The non-linear homotopy manufactures the
logarithm, and log|cos θ| = −log 2 − Σ(−1)^n cos(2nθ)/n is the first Fourier series whose cosine coefficients carry an
odd power of 1/n. Polynomials, and any linear homotopy of them, never do (Section 6).

---

## 1. The representation (Theorem 1, PROVED)

For p ≥ 1 the support function of B_p = {|x|^p+|y|^p ≤ 1} is the dual norm h_p(θ) = ‖(cos θ, sin θ)‖_q (Hölder with its
equality case), and Cauchy's formula for planar convex bodies gives P = ∫₀^{2π} h dθ. Hence

    P(p) = ∫₀^{2π} (|cos θ|^q + |sin θ|^q)^{1/q} dθ = π a₀(h_p)
         = 8 ∫₀¹ (1+u^q)^{1/q} (1+u²)^{-3/2} du                       (u = tan θ on [0, π/4])
         = 4√2 + 8τ ∫₀^∞ [(1+e^{-s})^τ − 1] e^{-τs} (1+e^{-2τs})^{-3/2} ds     (u = e^{-τs}).

*Evidence:* `f1_01_duality_check.py`. All four forms agree with the note's arclength formula (P/4 = 2 − εF(ε)) to
**1e-50** at p = 1.05, 1.25, 1.5, 2, 3, 7, 30.

*By-products:*
- **M₀ in three lines.** dP/dq at q = 1 equals 4∫₀^{π/2}[c log c + s log s − (c+s)log(c+s)] = 4[2(λ−1) − (2√2L − 2)]
  = 8(λ − √2L) = 8M₀ (checked to 1e-50).
- P′(2) = π(log 2 − ½), checked to 1e-51.
- **Polar duality is ε ↦ −ε.** With ε = 2/p − 1 we have q = 2/(1−ε), so h_{p(ε)} = 1/ρ_{p(−ε)}.

*Prior art:* the ingredients are classical: Cauchy's formula / mean width, the support function of a norm ball being
the dual norm, and Hurwitz 1901. I grepped the full texts of Montilla–Rodríguez–Linares (2607.04592) and
Rodríguez–Montilla (2607.09048): neither uses the support function, Cauchy's formula or Fourier series; both use
hypergeometric series of the arclength. Three web searches found no paper that applies the dual-norm/Cauchy route to
the squircle perimeter or its expansions. This was a limited search.

---

## 2. Closed forms of Fourier coefficients (question (a))

All PROVED by elementary integration unless marked. Every line is checked against quadrature to **≤ 4e-71**
(`f1_02_endpoint_coeffs.py`, k = 0..6).

| function | a_k (k ≥ 1) | a₀ |
|---|---|---|
| h₁ = max(\|c\|,\|s\|) (diamond support) | −(−1)^k 4√2/(π(16k²−1)) | 4√2/π |
| h_∞ = \|c\|+\|s\| (square support) | −8/(π(16k²−1)) | 8/π |
| ρ₁ = 1/(\|c\|+\|s\|) (diamond radial) | (8/π)(−1)^k Σ_{n odd > 4k} χ₈(n)/n  (tails of L(χ₈,1) = L/√2) | 4√2L/π |
| log ρ₁ | (−1)^k (2/(πk)) Σ_{j≥k} (−1)^j/(2j+1)  (tails of Leibniz) | log 2 − 4G/π |
| u = ∂_p log ρ_p\|₂ = ∂_p ρ_p\|₂ = ∂_p h_p\|₂ | **−1/(4k(4k²−1))** | log 2 − ½ |
| ρ_p^{−p} = \|c\|^p+\|s\|^p (all p) | 2γ_k(p), γ_m(p) = 2^{1−p}Γ(p+1)/(Γ(1+p/2+2m)Γ(1+p/2−2m)) | 2γ₀(p) |

*Proofs (sketch):*
- **h₁:** ∫_{−π/4}^{π/4} cos θ cos 4kθ = −(−1)^k√2/(16k²−1).
- **ρ₁:** J_m = ∫_{−π/4}^{π/4} cos 2mφ/cos φ satisfies J_m = 4 sin((2m−1)π/4)/(2m−1) − J_{m−1}, with J₀ = 2L.
- **log ρ₁:** expand log cos φ, then sum ∫_{−π/4}^{π/4} cos 2mφ cos 4kφ against χ₄. The mean comes from
  ∫₀^{π/4} log cos = G/2 − (π/4)λ. (The PSLQ line in the `.out` returns [−1,0,1,−4,0], i.e. a₀ = log 2 − 4G/π; the
  "check" line printed beside it has a sign typo.)
- **u:** multiply log|cos θ| = −λ − Σ(−1)^m cos 2mθ/m by cos²θ = (1+cos 2θ)/2.
- **γ_m:** Gradshteyn–Ryzhik 3.631.9, the generalized binomial C(p, p/2+2m)/2^p. Checked at p = 0.7, 1.3, 3, 5.5.

**p-derivatives of γ_m** (PROVED, checked to 1e-72):
- At p = 2: ∂_pγ_m(2) = 1/(4m(4m²−1)) for m ≥ 1 (from the pole of 1/Γ(1+p/2−2m)), and ∂_pγ₀(2) = ½ − log 2.
  Rational, apart from the constant term.
- At p = 1: ∂_pγ_m(1) = −(4/(π(16m²−1)))·[log 2 + 1 − 2Σ_{j=1}^{2m−1} 1/(2j−1) − 8m/(16m²−1)], and
  ∂_pγ₀(1) = (4/π)(log 2 − 1). Weight ≤ 1.

**Derivative at p = 1 of log ρ_p** (v = ∂_p log ρ_p|₁):
- In u = tan θ the log(1+u²) terms cancel: v = log(1+u) − u log u/(1+u), and v′(θ) = −log tan θ/(1 + sin 2θ). The
  letters are {0, −1, ±i}, so the coefficients have weight ≤ 2 in {1, λ, π, π², πλ, λ², G} (PROVED).
- Integrating by parts once, then w = tan(θ − π/4), gives
  a_k(v) = ((−1)^k/(2πk)) ∫_{−1}^{1} log((1+w)/(1−w)) Im(1+iw)^{4k}/(1+w²)^{2k} dw. So for k ≥ 1,
  a_k(v) ∈ (−1)^k π/4 + ℚ + ℚ/π (structure PROVED).
- The values are NUMERICAL (PSLQ, 70 digits, basis of 7; residual ≤ 2.7e-68; max coefficient 16, 16, 48, 720,
  5040; the growth is the lcm of the denominators):

| k | a_k(v) |
|---|---|
| 0 | 4G/π + log 2 − π/4 |
| 1 | π/4 − 1 |
| 2 | −π/4 + 4/3 − 2/π |
| 3 | π/4 − 7/5 + 16/(9π) |
| 4 | −π/4 + 152/105 − 98/(45π) |

**First derivatives at p = 1 of the support coefficients vanish** (one-sided, p → 1+). By Theorem 3,
a_k(p) − a_k(1) = O(τ²). The second and third derivatives are where ζ(2) and ζ(3) live (Section 3).

---

## 3. Diamond end: every Fourier coefficient carries ζ(2)τ² + ζ(3)τ³ (Theorem 3, PROVED)

**Exact representation.** With A_k = a_k(h₁) = −(−1)^k4√2/(π(16k²−1)):

    a_k(τ) − A_k = (8/π) τ ∫₀^∞ [(1+e^{-s})^τ − 1] G_k(τs) ds,
    G_k(x) = cos(4kθ) cos²θ sin θ |_{θ = arctan e^{-x}} = (−1)^k cos(2k·gd x) · 2^{-3/2} e^{x/2} sech^{3/2} x .

**Asymptotic expansion.** a_k − A_k ~ Σ_{n≥2} c_n(k) τ^n with

    c_n(k) = (8/π) Σ_{j=1}^{n−1} (1/j!) I_{n−1−j, j} g_{k, n−1−j},
    I_{m,j} = ∫₀^∞ s^m log^j(1+e^{-s}) ds,   g_{k,m} = [x^m] G_k(x).

*Proof.* Use 0 ≤ ℓ(s) := log(1+e^{−s}) ≤ e^{−s}.
- Remainder in τ: (1+e^{−s})^τ − 1 = Σ_{j≤J} τ^j ℓ^j/j! + O(τ^{J+1}ℓ^{J+1}).
- Remainder in x: G_k is C^∞ on [0, ∞) with bounded derivatives, so G_k(x) = Σ_{m<M} g_{k,m}x^m + O(x^M).
- Both remainders integrate to O(τ^{min(J, M)+2}) (Watson-lemma type).

**The first coefficients:**
- Ingredients:
  - I_{0,1} = η(2) = π²/12;
  - I_{1,1} = η(3) = ¾ζ(3);
  - I_{0,2} = ∫₀¹ log²(1+x)dx/x = ζ(3)/4;
  - g_{k,0} = (−1)^k√2/4 and g_{k,1} = (−1)^k√2/8, independent of k because sin 4kθ = 0 at θ = π/4.
- Hence:
  - **c₂(k) = (8/π)g_{k,0}I_{0,1} = (−1)^k(√2/π)ζ(2)**;
  - **c₃(k) = (8/π)[g_{k,0}I_{0,2}/2 + g_{k,1}I_{1,1}] = (8/π)(−1)^k√2(1/32 + 3/32)ζ(3) = (−1)^k(√2/π)ζ(3)**. ∎
- **τ⁴ (PROVED):** g_{k,2} = −(5√2/32 + (√2/2)k²)(−1)^k and I_{2,1} = 2η(4) = (7/4)ζ(4). Use the classical
  - I_{0,3} = ∫₀¹ log³(1+x)/x = 6ζ(4) + (3/2)ζ(2)λ² − (21/4)ζ(3)λ − λ⁴/4 − 6Li₄(½);
  - I_{1,2} = −(15/4)ζ(4) + 4Li₄(½) + λ⁴/6 − ζ(2)λ² + (7/2)ζ(3)λ (PSLQ in `.out` agrees).

  The Li₄(½), λ⁴, ζ(2)λ² and ζ(3)λ parts cancel exactly between g_{k,0}I_{0,3}/6 and g_{k,1}I_{1,2}/2, which leaves
  **c₄(k) = −(−1)^k(√2/π)(33/16 + 7k²)ζ(4)**.
- **τ⁵ (NUMERICAL):** c₅(k) = −(−1)^k(√2/(4π))[(3+16k²)ζ(5) + (7+16k²)ζ(2)ζ(3)]. PSLQ, k = 0..3, basis
  {ζ(5), ζ(2)ζ(3)} at 60 digits: residual ≤ 2.5e-60, max coefficient 151. The k²-part is equivalent to the classical
  ∫₀¹ log²x log²(1+x)/x dx = 2ζ(2)ζ(3) − (29/8)ζ(5).

**Pure weight (PROVED).** g_{k,m} ∈ (−1)^k √2·ℚ[k²], and I_{m,j} is an alternating MZV of weight m+j+1. So
(−1)^kπc_n(k)/√2 is a ℚ[k²]-combination of weight-n alternating MZVs, of degree ≤ ⌊(n−2)/2⌋ in k². The individual
I_{m,j} contain Li₄(½), but the combinations seen so far reduce to single-zeta products:
- PROVED for n ≤ 4;
- NUMERICAL for n = 5 (k ≤ 3);
- for k = 0, all n ≤ 10 from D11's exact table (next item).

**Check against D11 (`f1_03_diamond_fourier.py`).** P/(4√2) = (1+η)Σf_nη^n with η = −τ, so its τ-coefficient is
e_n = (−1)^n(f_n + f_{n−1}).
- **Every e_n, n = 1..10, is homogeneous of weight n.** D11's f_n themselves mix weights {0, 2, 3, …, n}.
- πc_n(0) = 4√2e_n to ≤ 4e-59 for n = 2..8. This is an independent re-derivation of D11's a₂…a₈ by a different
  representation.
- In τ: e₂ = ζ(2)/4, e₃ = ζ(3)/4, e₄ = −33ζ(4)/64, e₅ = −(3ζ(5) + 7ζ(2)ζ(3))/16.

**Remainder check (NUMERICAL).** (a_k − A_k − Σ_{n≤6} c_nτ^n)/τ⁷ → c₇(k) as τ = .02, .01, .005:
- k = 0: 18.17, 19.23, 19.74 → c₇ = 20.23;
- k = 1: −136.0, −152.6, −160.7 → c₇ = −168.4;
- k = 2: 709.6, 880.1, 963.4 → c₇ = 1044.5.

The expansion is non-uniform in k (the k = 5 column converges slowly), as expected.

**Corollary 3a (Dirac comb, PROVED with the same argument).** For φ ∈ C³ on the circle, put θ_j = π/4 + jπ/2. Then

    ∫₀^{2π} (h_p − h₁) φ dθ = (√2/4)(ζ(2)τ² + ζ(3)τ³) Σ_{j=0}^{3} φ(θ_j) + O(τ⁴).

The φ′ terms cancel because h_p is symmetric about each θ_j. Numerical check with the non-symmetric
φ = exp(sin θ + ½cos 3θ): residual/τ⁴ = −4.12, −3.99, −3.92, −3.89 at τ = .04, …, .005, so it is bounded
(`f1_06_comb_and_figure.py`). For φ = cos 4kθ this is Theorem 3.

**Proposition 3b (uniform boundary-layer form, NUMERICAL plus derivation).**

    a_k = A_k + (−1)^k (√2/π)[1/(4k²) − πτ/(2k sinh(2πkτ))] + O(τ³),

uniformly in k.
- Derivation: replace G_k(τs) by g_{k,0}cos(2kτs), and use ∫₀^∞ ℓ(s)cos(bs)ds = Σ(−1)^{n+1}/(n²+b²)
  = 1/(2b²) − π/(2b sinh πb).
- Evidence: sup over k ≤ 200 of the error is 0.518τ³, 0.530τ³, 0.536τ³ at τ = .02, .01, .005.
- For kτ ≫ 1 the coefficients become exponentially small (k = 200, τ = .01: a_k = −2.4e-10 against A_k = −2.8e-6),
  which is how the kinks of h₁ get smoothed.
- **Consistency check (PROVED given 3b).** Insert 3b into Hurwitz's identity. Use the η-quotient asymptotics
  Σ_k 1/(k sinh 2πkτ) = π/(12τ) − log 2 + πτ/6 + O(e^{−c/τ}) (checked to 1e-27) and
  ∫₀^∞(1/y − csch y)²dy = 2 log 2 − 1 (proof: sine-transform Parseval with
  ∫₀^∞(1/y − csch y)sin ωy dy = (π/2)(1 − tanh(πω/2))). The τ-coefficient of Σ(16k²−1)(a_k² − A_k²) is then
  −(8/π)log 2 + (4/π)(2 log 2 − 1) = **−4/π**. This is exactly the value forced by P = √2(4 + O(τ²)) and
  A = 2(1+τ) + O(τ²). The log 2 terms cancel, and termwise differentiation would have given 0.

**Figure:** `f1_boundary_layer.png`. Left: (h_p − h₁)/τ² concentrates at θ = π/4. Right:
(−1)^kπ(a_k − A_k)/(√2τ²) starts at ζ(2) + ζ(3)τ for every small k and follows the 1/sinh profile (dotted).

---

## 4. Circle end p = 2 (questions (a), (b))

**Cumulant form (PROVED).** Use ε = 2/p − 1. Put C = cos², a = log|cos|, b = log|sin|, and κ₁ = Ca + (1−C)b,
κ₂ = C(1−C)(a−b)², κ₃ = C(1−C)(1−2C)(a−b)³ (the cumulants of the two-point law). Then

    log h_p = κ₁ε + κ₂ε² + (κ₂ + ⅔κ₃)ε³ + O(ε⁴)        (checked numerically at a point to 20 digits).

**P₀…P₃ by Beta derivatives (PROVED).** Integrating monomials by Beta derivatives gives
P/π = 2 + (1 − 2λ)ε + (λ² − λ − ¾ + 11π²/96)ε² +
(−(23/32)ζ(3) − (11/96)π²λ − λ³/3 + λ²/2 + ¾λ + 1/16 + 37π²/384)ε³. This is **exactly** the result of
`05_final_verify/circle_end` (symbolic difference 0 for P₀…P₃), which used the arclength. So it is a second,
independent exact derivation (`f1_04_circle_end.py`).

**First-order data are locked to even weight.**
- α_k := [κ₁]_k = 1/(2k(4k²−1)) (= −2× the coefficients of u; PROVED).
- Hurwitz's identity (Hurwitz 1901, checked here for the squircle at p = 1.2 to 1e-29 and at p = 1.5 to 1e-12 with
  K = 120) says
  P² − 4πA = 2π² Σ_{k≥1}(16k²−1) a_k(h_p)².
- At order ε²: **Σ(16k²−1)α_k² = 7π²/48 − 1** (sympy summation, exact). It equals D₂/(2π²) computed from P and
  A = 4Γ(1+1/p)²/Γ(1+2/p) (difference 0). The log 2 terms in P² and 4πA cancel, as they must: an even rational summand
  can only give ℚ + ℚπ² + … (residue theorem).

**Second-order data carry harmonic numbers, and cubic Parseval gives ζ(3).**
- Write X = log|2 sin 2θ| = −Σ cos 4nθ/n and Y = log|tan θ| = −2Σ_{m odd} cos 2mθ/m. Then
  κ₂ + κ₁²/2 = ⅛(X−2λ)² − ¼(X−2λ)cos 2θ·Y + ⅛(1+sin²2θ)Y².
- [X²]_k = (2H_k − 1/k)/k, [Y²]_k = 4O_{2k}/k with O_{2k} = Σ_{i≤k} 1/(2i−1), [cos 2θ·Y]_k = −4k/(4k²−1).
- So the second-order coefficient is

      β_k = (2H_k − 1/k)/(8k) + λ/(2k) − ¼C_k − 2λk/(4k²−1) + ¾O_{2k}/k − (Y_{k−1} + Y_{k+1})/32,

  with Y₀ = π²/2, Y_j = 4O_{2j}/j, and C_k = [X cos 2θ Y]_k an explicit digamma/harmonic expression (in the script).
  It matches quadrature to ≤ 6e-52 for k = 1..6.
- At order ε³, Hurwitz gives the **explicit Parseval sum of squircle coefficients**

      Σ_{k≥1} (16k²−1) α_k β_k = D₃/(4π²) = −(15/32)ζ(3) + 43π²/384 − (7/48)π² log 2 + log 2 − 5/16
                                 = −0.0752856509449585376492485782…

- Status:
  - The right side is exact (sympy, from P₃ and A). The identity follows from Hurwitz by termwise differentiation:
    PROVED modulo the routine dominated-convergence step (|∂^j_ε a_k| = O(k⁻³ log² k) near ε = 0; not written out).
  - Numerically, Σ_{k≤16000} plus a 7-term tail fit gives a difference of **1.7e-26**; the tail fit is stable to
    3e-22.
  - ζ(3) enters only through Euler sums: Σ H_k/k² = 2ζ(3) and Σ O_{2k}/k² = (7/4)ζ(3).

**Every odd ζ from the first-order velocity and a polynomial (PROVED).** Let ⟨f,g⟩ = (2/π)∫₀^{π/2}fg and
β_{2j}(θ) = Σ_{n≥1} cos 4nθ/n^{2j} = (−1)^{j+1}(2π)^{2j}B_{2j}(2θ/π)/(2(2j)!). This is a Bernoulli polynomial, and
β₂ = ((4θ−π)² − π²/3)/4 is Aritro's x² with x = 4θ − π. Put Z_j = Σ 1/(n^{2j+1}(4n²−1)). From
1/(n^{2j+1}(4n²−1)) = −1/n^{2j+1} + 4/(n^{2j−1}(4n²−1)) we get ζ(2j+1) = 4Z_{j−1} − Z_j, with Z₀ = 2λ − 1 = 4⟨u,1⟩
and Z_j = −8⟨u, β_{2j}⟩. Hence

    ζ(3)    = 16⟨u,1⟩ + 8⟨u,β₂⟩ = (4/π) ∫₀^{π/2} u(θ) [ (4θ−π)² + 8 − π²/3 ] dθ,
    ζ(2j+1) = 8 ⟨u, β_{2j} − 4β_{2j−2}⟩            (j ≥ 2).

Checked to ≤ 1.6e-61 for ζ(3), ζ(5), ζ(7), ζ(9), ζ(11) (`f1_05_odd_zeta_parseval.py`). The general formula for
⟨u, x^{2j}⟩ (x = 4θ − π) is checked for j ≤ 4. *Novelty:* this is Euler-1772-type log-cos integration (see
`01_survey/T3_fourier_elementary` §1.1–1.2). What is new is only the reading of the integrand as the squircle's
velocity at the circle.

---

## 5. Answers to the specific questions

- **(a)** See Section 2.
  - At p = 2 every p-derivative coefficient is rational (u: −1/(4k(4k²−1)); γ_m: 1/(4m(4m²−1))), apart from a log 2
    in the constant term.
  - At p = 1 the first derivatives are weight ≤ 2 (log ρ: π, ℚ, ℚ/π, plus G and λ in a₀) or weight ≤ 1 (γ_m).
    The support coefficients have zero first derivative.
  - No zeta value appears at first order at either end. That is the parity lock in action.
- **(b)** Yes, in three precise ways:
  - (i) the diamond theorem (each support Fourier coefficient carries ζ(3)τ³; k = 0 is the perimeter);
  - (ii) the cubic Hurwitz–Parseval sum at the circle, −(15/32)ζ(3) + …;
  - (iii) first-order velocity × polynomial, which gives ζ(3) exactly, and every ζ(2j+1).
  
  The homotopy supplies the logarithm, i.e. the "wrong-parity" Clausen-type coefficients that polynomial Fourier
  series miss. At the circle it appears as u and as the H_k in β_k. At the diamond it appears as the boundary-layer
  profile log(1+e^{−s}) = Σ(−1)^{n+1}e^{−ns}/n, the cosine series of log(2cos(x/2)) continued to x = is.
- **(c) d²/dp², d³/dp³ at p = 1 through Fourier coefficients.**
  - Yes. P = πa₀(h_p) and τ = d − d² + d³ − …, so Theorem 3 gives P/√2 = 4 + (π²/6)d² + (ζ(3) − π²/3)d³.
  - Hence P″(1) = 2√2ζ(2) = √2π²/3 and P‴(1) = 6√2(ζ(3) − π²/3), matching the ledger (R2/D11).
  - The same holds for every coefficient: d²a_k/dp² = (−1)^k√2π/3 and d³a_k/dp³ = 6(−1)^k(√2/π)(ζ(3) − π²/3).
  - The −π²/3 is not a separate constant: it is the ζ(2)τ² term re-expanded in d. In τ the series is pure weight.
  - Caveat: the radial-function Parseval (area) cannot be differentiated termwise at p = 1, because the boundary layer
    makes the expansion non-uniform in k; Proposition 3b handles that correctly.

---

## 6. Why it works where the Fourier approach failed (in the language of the paper)

1. **The lock.**
   - Integrating a polynomial against cos nx by parts gives each power 1/n^j attached alternately to sin nx and
     cos nx at the endpoints.
   - For x², x⁴ only even powers survive in the cosine series, so Parseval or x = π gives ζ(2), ζ(4).
   - For x³ the odd powers 12/n³ come attached to sin nx and die at x = π (his π³ = 0). Making the function even
     (|x|³) moves to a cosine series and again leaves only even powers.
   - In general a piecewise-polynomial coefficient is a finite sum Σ_j e^{inβ_j}P_j(1/n). Odd powers in a cosine
     series appear only with a sin(nβ) phase, which gives β(3) = π³/32, never ζ(3).
2. **A linear homotopy stays locked.** Mixing the diamond and the square linearly (Minkowski sum
   (1−t)◇ + t□, an octagon) mixes their support functions linearly:
   - a_k(t) = (1−t)A_k + tB_k ∈ ℚ(√2)/π;
   - P(t) = 4√2(1−t) + 8t;
   - A(t) = 2 + 4t − 2t² (mixed area V(◇,□) = 4).
   
   No zeta value can ever appear. This is exactly his x² → x⁴ observation.
3. **The squircle unlocks it.**
   - The exponent is the homotopy parameter, and ∂_p|x|^p = |x|^p log|x|. The velocity u at the circle contains
     log|cos θ|, whose cosine coefficients ±1/n are a bare odd power.
   - u's coefficients −1/(4n(4n²−1)) = −Σ_{j≥0} 4^{−j−2} n^{−2j−3} form an *odd* function of n with
     infinitely many odd powers. No piecewise polynomial has this.
   - Parseval of u against x² gives Σ 1/(n³(4n²−1)) = −ζ(3) + 4(2 log 2 − 1).
4. **At the diamond.**
   - The rounding lives in a layer of width τ around each edge normal, where (1+tan^qθ)^{1/q} = (1+e^{−s})^τ
     = 1 + τ log(1+e^{−s}) + …. The 1/q-th root brings out the logarithm.
   - Integrating its series Σ(−1)^{n+1}e^{−ns}/n against 1 and against s gives Σ(−1)^{n+1}/n² = ζ(2)/2 and
     Σ(−1)^{n+1}/n³ = ¾ζ(3): ∫₀^∞ s e^{−ns} ds = 1/n² plays the role of his integration by parts.
   - The layer is a point mass in θ, so it shows up in *every* Fourier coefficient with the same ζ(2) and ζ(3).

This agrees with the survey framing in `01_survey/T3_fourier_elementary` §5 (the "phase–order lock", log-sine
functions). T3/T4b had already identified the *area* as the log-sine moment generating function (Γ-function). F1 is
the perimeter / support-function / Fourier-coefficient side, and its constants are not Γ-function constants.

---

## 7. Evidence table

| claim | status | evidence | script |
|---|---|---|---|
| P = ∫h_p = πa₀(h_p), three integral forms | PROVED | 1e-50 vs arclength, 7 values of p | f1_01 |
| M₀ from dP/dq at q = 1; P′(2) = π(λ − ½) | PROVED | 1e-50, 1e-51 | f1_01 |
| endpoint coefficients h₁, h_∞, ρ₁, log ρ₁, u, γ_m | PROVED | ≤ 4e-71, k ≤ 6 | f1_02 |
| ∂_pγ_m at p = 1, 2 | PROVED | ≤ 1e-72 | f1_02 |
| a_k(v), k = 0..4 (log ρ derivative at p = 1) | structure PROVED, values NUMERICAL | PSLQ, 70 digits, residual ≤ 2.7e-68, max coefficient ≤ 5040, basis of 7 | f1_02 |
| c₂ = (−1)^k√2ζ(2)/π, c₃ = (−1)^k√2ζ(3)/π, all k | PROVED | 60 digits, k ≤ 5 | f1_03 |
| c₄ = −(−1)^k(√2/π)(33/16 + 7k²)ζ(4) | PROVED | PSLQ [16, 33+112k²], residual ≤ 2e-59 | f1_03 |
| c₅(k) formula | NUMERICAL | PSLQ, k ≤ 3, residual ≤ 2.5e-60, max coefficient 151 | f1_03 |
| pure weight in τ; D11 e_n homogeneous, n ≤ 10 | PROVED (structure), checked | πc_n(0) = 4√2e_n to 4e-59, n ≤ 8 | f1_03 |
| Dirac-comb corollary | PROVED (sketch) | residual/τ⁴ bounded (−3.9) | f1_06 |
| uniform 1/sinh form | NUMERICAL | sup error ≈ 0.53τ³ | f1_03 |
| ∫(1/y − csch y)² = 2 log 2 − 1; deficit slope −4/π | PROVED | 25 digits; exact match | f1_05 |
| P₀…P₃ at the circle via the support function | PROVED | symbolic difference 0 vs circle_end | f1_04 |
| Σ(16k²−1)α_k² = 7π²/48 − 1 | PROVED | sympy exact | f1_04 |
| cubic Parseval Σ(16k²−1)α_kβ_k = −(15/32)ζ(3) + … | PROVED mod dominated convergence; NUMERICAL | 1.7e-26 (fit stable 3e-22) | f1_04 |
| β_k explicit formula | PROVED | ≤ 6e-52 vs quadrature | f1_04 |
| ζ(2j+1) = 8⟨u, β_{2j} − 4β_{2j−2}⟩; ζ(3) single integral | PROVED | ≤ 1.6e-61 for ζ(3..11) | f1_05 |
| Hurwitz identity for the squircle | classical | 1e-29 (p = 1.2), 1e-12 (p = 1.5, K = 120) | f1_05 |

Scripts (in this folder): `f1_01_duality_check.py`, `f1_02_endpoint_coeffs.py`, `f1_03_diamond_fourier.py`,
`f1_04_circle_end.py`, `f1_05_odd_zeta_parseval.py`, `f1_06_comb_and_figure.py`; outputs `*.out`; figure
`f1_boundary_layer.png`. Runtime: under a minute each.

---

## 8. Open questions

1. **Single-zeta reduction for all k.** For k ≥ 1, are all τ^n coefficients c_n(k) single-zeta products? This is
   proved for n ≤ 4 and NUMERICAL for n = 5; k = 0 is D11.
   - Find the generating function in k. The k-dependence enters only through cos(2k·gd x), i.e. through the
     derivatives of cos 4kθ at the edge normals.
2. **A rigorous uniform-in-k bound** for Proposition 3b.
3. **Corner-angle family.** Does the same "comb at the edge normals with masses ζ(2)τ² + ζ(3)τ³" hold for other
   polygons rounded by an ℓ_q-type norm? Which constants replace √2/4 (they should lie in ℚ(cos(α/2)), linking to the
   ledger's corner-angle rule)?
4. **Exact closed form of a_k(v) for all k.** The numbers above suggest rationals built from odd harmonic sums.
5. **Higher orders at the circle.** Does the ε⁴ Hurwitz–Parseval term (ζ(3)λ, π⁴) have a similarly clean Euler-sum
   reading?
6. **Divergence in τ.** The τ-series diverges (R2's zero radius). Do the Borel data look simpler in τ than in d, given
   pure weight?

---

## 9. Suggested use in the paper (Non-linear Homotopies subsection, in Aritro's voice)

"Let h_p be the support function of the squircle. Cauchy's formula turns the perimeter into one Fourier coefficient,
P(p) = πa₀(h_p). The linear homotopy from the diamond to the square changes these coefficients linearly, and no zeta
value appears, exactly as with x² → x⁴. The squircle instead changes the *exponent*. Since ∂_p|x|^p = |x|^p log|x|,
it creates logarithms, whose Fourier coefficients 1/n have the parity my polynomials never had. As a result every
Fourier coefficient of h_p near the diamond is

a_k = (−1)^k(√2/π)[−4/(16k²−1) + ζ(2)τ² + ζ(3)τ³ + …],  τ = 1 − 1/p,

and for k = 0 this gives P/√2 = 4 + ζ(2)τ² + ζ(3)τ³ + …. The non-linear homotopy produces ζ(2) and then ζ(3),
one order apart."
