# PUR1: a purity criterion for layer integrals

**Outcome.** The conjecture holds as stated, after one repair in (e) (an extra pure term) and one sharpening in (g).

- (a) and (b) are PROVED, with an explicit remainder.
- The criterion is a one-line consequence of the Beta reflection identity: evaluate it on its fixed line λ = μ.
- For the squircle it reproves the paper's all-order pure-weight theorem. c₂…c₇ come out exactly; c₆ and c₇ agree with quadrature plus Richardson to 37 digits.
- For the soft-max it isolates 40γ₃₂ at τ⁶ exactly as SM6 found, and reproduces SM6's c₆ and c₇ to 51 and 47 digits.
- The designed square of (f) is a **convex body**, proved for every τ ∈ (0,1]. Its perimeter is 16 − 8Γ(1−τ/2)²/Γ(1−τ) plus an explicit exponentially small term (checked to 47–50 digits).
- What is *not* proved: that the alternating sums γ_{2j+1,2} really are not polynomials in zeta values. That is CONDITIONAL on standard independence conjectures. It is not even known unconditionally that γ₁₂ ∉ ℚζ(4).

Scripts in this folder (mpmath 1.3, sympy 1.14, throwaway venv, deleted at the end):

| script | output | what it does |
|---|---|---|
| `pur_core.py` | (library) | S-series, coefficient maps (a) and (b), numeric γ_ab |
| `pur_squircle.py` | `pur_squircle_output.txt` | (b) vs (a), squircle c₁…c₇ |
| `pur_quad.py` | `pur_quad_{squircle,softmax,test_e,rhombus}_output.txt` | exactness checks, Richardson extrapolation |
| `pur_softmax_e.py` | `pur_softmax_e_output.txt` | (d), (e), test example, rhombus |
| `pur_designed.py` | `pur_designed_output.txt` | (f), designed rectangle |

The small `.pkl` files pass predicted coefficients between scripts.

---

## 0. Notation

- g(s) = log(1+e^{−2s}).
- B(s,τ) := (1+e^{−2s})^τ − 1 = e^{τg} − 1.
- γ_{ab} = (1/(a!b!)) ∫₀^∞ s^a g^b ds.
- G(λ,μ) = ∫₀^∞ e^{λs}[(1+e^{−2s})^μ − 1] ds = Σ_{a≥0, b≥1} γ_{ab} λ^a μ^b.

G is analytic on {Re λ < 2} × ℂ, because |(1+x)^μ − 1| ≤ |μ| log(1+x) e^{|μ| log 2} ≤ |μ| 2^{|μ|} x.

**Lemma 5 of LAY1 in the variable x = λ − μ.** Put
Q(x,μ) := Γ(1−(μ+x)/2) Γ(1−(μ−x)/2) / Γ(1−μ),
so that
log Q = Σ_{k≥2} (ζ(k)/k) [((μ+x)/2)^k + ((μ−x)/2)^k − μ^k].
Euler's γ cancels. Then, with S(μ+x, μ) = 2μ(1 − Q)/(μ² − x²),

  **G(μ+x, μ) + G(μ−x, μ) = S(μ+x, μ).**  (R)

*Check of the rewriting.* Write Γ(−z) = −Γ(1−z)/z and 1/Γ(−μ) = −μ/Γ(1−μ). Then ½·Γ(−λ/2)Γ(λ/2−μ)/Γ(−μ) = −2μQ/(μ²−x²), and ½[2/λ + 2/(2μ−λ)] = 2μ/(μ²−x²). The rewriting therefore has no new content. At a test point, the S-series matches the Gamma form up to its truncation (`pur_squircle_output.txt`).

Q − 1 vanishes on x = ±μ, so (1−Q)/(μ²−x²) is a power series.

**Weight grading.** Give ζ(k) weight k and give x and μ weight −1. Then every monomial of log Q, and hence of Q, has total weight 0, and S has total weight 1.

---

## 1. (a) The expansion and its remainder: PROVED

**Hypotheses.**
- 0 < V ≤ ∞, and φ is C^{N−1} on [0, V₀] for some finite V₀ ∈ (0, V].
- |φ(v)| ≤ C e^{βv} on [0, V) for some β ≥ 0. This is automatic for V < ∞.

**Proposition A.** For 0 < τ ≤ min(1, 1/β),
  |F(τ) − Σ_{n=0}^{N} F_n τⁿ| ≤ C_N τ^{N+1} + C′ e^{−V/τ},
where F_n = Σ_{k+b=n−1, b≥1} φ_k k! γ_{kb}. Equivalently,
  F ~ Σ_k φ_k k! Σ_{b≥1} γ_{kb} τ^{k+1+b}.
The term C′e^{−V/τ} is absent when V = ∞. Since e^{−V/τ} = O(τ^{N+1}), the remainder is O(τ^{N+1}) in all cases.

*Proof.*
1. Substitute v = τs: F = τ ∫₀^{V/τ} φ(τs) B(s,τ) ds.
2. Bounds on B. Since 0 ≤ g ≤ min(log 2, e^{−2s}),
   - 0 ≤ B = e^{τg} − 1 ≤ τg e^{τg} ≤ 2τ e^{−2s};
   - |B − Σ_{b<L} τ^b g^b/b!| ≤ (τg)^L e^{τg}/L! ≤ 2τ^L e^{−2Ls}/L!.
3. Taylor remainder with K = N−1. Let P_K(v) = Σ_{k≤K−1} φ_k v^k and R_K = φ − P_K. On [0, V₀], Taylor gives |R_K| ≤ (sup|φ^{(K)}|/K!) v^K. On [V₀, V), |R_K| ≤ |φ| + |P_K| ≤ C″ e^{βv} (v/V₀)^K. So |R_K(v)| ≤ C_K v^K e^{βv} on [0, V).
4. The remainder piece:
   τ ∫ |R_K(τs)| B ds ≤ 2C_K τ^{K+2} ∫₀^∞ s^K e^{(βτ−2)s} ds ≤ 2C_K K! τ^{N+1},  using βτ ≤ 1.
5. The polynomial pieces. Each is φ_k τ^{k+1} [∫₀^∞ s^k B ds − ∫_{V/τ}^∞ s^k B ds].
   - The tail is ≤ 2τ ∫_{V/τ}^∞ s^k e^{−2s} ds = O(e^{−V/τ}).
   - By step 2 with L = N − k, ∫₀^∞ s^k B ds = Σ_{b=1}^{N−k−1} τ^b k! γ_{kb} + O(τ^{N−k}).
6. Collect powers of τ. ∎

For comparison, γ_{kb} ≤ 1/(b!(2b)^{k+1}).

---

## 2. (b) The criterion: PROVED

**Theorem B.** Write φ = e^v ψ and ψ(v) ~ Σ ψ_k v^k at 0. Under the hypotheses of (a), as formal power series in τ (each coefficient a finite sum),

  F ~ τ Σ_k ψ_k τ^k ∂_λ^k G(λ,τ)|_{λ=τ}.  (B1)

For even k,

  ∂_λ^{2k} G(λ,τ)|_{λ=τ} = ½ ∂_x^{2k} S(τ+x, τ)|_{x=0} = τ ∂_x^{2k}[(1−Q(x,τ))/(τ²−x²)]_{x=0}.  (B2)

So if ψ_k = 0 for every odd k (that is, e^{−v}φ is even as a Taylor series),

  **F ~ τ² Σ_k ψ_{2k} τ^{2k} ∂_x^{2k}[(1−Q(x,τ))/(τ²−x²)]_{x=0}.**  (B3)

Every coefficient of (B3) is a polynomial in ζ(2), ζ(3), … whose coefficients are ℚ-linear in the ψ_{2k}. If the ψ_{2k} are rational, the τⁿ coefficient is homogeneous of weight n. No log 2, Li_n(1/2) or alternating sum can occur.

*Proof.*
1. **(B1).** By (a), F ~ τ Σ_k φ_k τ^k ∂_λ^k G(0,τ), because ∂_λ^k G(0,μ) = Σ_b k! γ_{kb} μ^b.
   - Insert φ_k = Σ_{i+j=k} ψ_i/j!, which is the Cauchy product of e^v and ψ.
   - Then Σ_j (τ^j/j!) ∂_λ^{i+j} G(0,τ) = ∂_λ^i G(τ,τ). This is Taylor's theorem in λ, which converges since |τ| < 2.
   - The rearrangement is legitimate coefficientwise: ∂^{i+j}G(0,τ) = O(τ), so the τⁿ coefficient involves only i + j ≤ n − 2.
2. **(B2).** Differentiate (R) 2k times in x and set x = 0. The two terms on the left are equal at x = 0 and the sign is (−1)^{2k} = +1, so 2∂^{2k}G(μ,μ) = ∂_x^{2k}S|_{x=0}.
   - For odd derivatives the same step gives 0 = 0. This is why the reflection says nothing about odd ψ_k.
3. **Weight.** By §0, a weight-w monomial of S has degree w − 1 in (x, τ).
   - ∂_x^{2k} at x = 0 leaves degree w − 1 − 2k in τ.
   - Multiplying by τ^{1+2k} gives τ^w.
   - Hence the ζ-weight equals the power of τ. ∎

**Check, not proof.** (a) with γ_ab by quadrature and (B3) symbolically agree for the squircle to ≤ 3·10⁻⁷⁰, n = 1…7 (`pur_squircle_output.txt`), and for the soft-max to ≤ 1.5·10⁻⁷⁰.

**Special case ψ ≡ 1** (closed form): F ~ τG(τ,τ) = 1 − Γ(1−τ/2)²/Γ(1−τ).

---

## 3. (c) The squircle: PROVED (identity and all-order purity) + NUMERICAL

**Exact change of variables.**
- On θ = π/4 + u, 0 ≤ u ≤ π/4, LAY1 §3 gives h_p = m(1 + r^{1/τ})^τ, with m = (cos u + sin u)/√2 and r = tan(π/4 − u) = e^{−2v}, where v = artanh(tan u).
- dv/du = sec 2u.
- From tan u = tanh v: cos 2u = 1/cosh 2v, cos u = cosh v/√(cosh 2v), and cos u + sin u = e^v/√(cosh 2v).
- Hence m du = e^v/(√2 cosh^{3/2} 2v) dv exactly, and by D₄ symmetry
  P = 8∫₀^∞ φ(v)(1+e^{−2v/τ})^τ dv,  φ = e^v/(√2 cosh^{3/2}(2v)).
- ∫₀^∞ φ = ∫₀^{π/4} m du = 1/√2, so P = 4√2 + 8F with V = ∞.

ψ = 1/(√2 cosh^{3/2} 2v) is even as a function. Its coefficients ψ_{2k} lie in √2·ℚ: √2/2, −3√2/2, 13√2/4, −379√2/60, …. Theorem B therefore gives P/√2 to every order as a pure-weight polynomial in zeta values.

Coefficients from (B3) (sympy, symbolic; difference from the known values is 0 after ζ(2k) → π^{2k}):

| n | c_n = 8F_n | agrees with |
|---|---|---|
| 2 | √2 ζ(2) | paper |
| 3 | √2 ζ(3) | paper |
| 4 | −(33√2/16) ζ(4) | paper, LAY1 |
| 5 | −√2(¾ζ(5) + (7/4)ζ(2)ζ(3)) | paper, SM5 |
| 6 | √2(4039ζ(6)/256 − 7ζ(3)²/8) = 20.91147814643259204894245555535839107867 | SM6 |
| 7 | √2(129ζ(7)/16 + 141ζ(2)ζ(5)/16 + 1071ζ(3)ζ(4)/64) = 63.54451807765424743136703150673447368701 | SM6 |

**Numerics** (`pur_quad_squircle_output.txt`).
1. Exactness of the v-form against θ-quadrature of h_p, dps 50: P(v) − P(θ) = −1.1·10⁻⁵⁰ at τ = 0.5, 0.2, 0.1 and 0.05 (rounding level). At τ = ½, P = 2π to 50 digits, as it must.
2. Richardson. F(τ) by v-quadrature at dps 120 for τ = 1/(20·2^i), i = 0…14. F is computed directly, so the constant 4√2 is never cancelled.
   - (F − Σ_{n≤5}F_nτⁿ)/τ⁶ extrapolates to c₆/8. The difference from the exact c₆ falls to 7.8·10⁻³⁸ (order 14).
   - c₇ is reached to 2·10⁻³⁶.
   - So c₆ and c₇ are confirmed to ≥ 36 digits.

---

## 4. (d) The soft-max: PROVED (given SM6's exact reduction) + NUMERICAL

SM6 §2 gives the setup:
- I/8 − 1 = F + O(e^{−√2/τ}), with φ = 1 + v/√(1−v²) and V = 1/√2.
- ψ = e^{−v}φ = 1 − ½v² + (5/6)v³ − (5/8)v⁴ + (79/120)v⁵ − … (sympy).
- ψ₁ = 0 and ψ₃ = 5/6. The ψ₃ value comes from 1/3 (from e^{−v}(1+v)) plus 1/2 (from v³/2).

By (B1), c_n = 8F_n = 8·[even part via (B3)] + 8·[odd part]. The odd part is
  ψ₃τ⁴ ∂³G(τ,τ) + ψ₅τ⁶ ∂⁵G(τ,τ) + …,  with ∂^kG(τ,τ) = Σ_{a≥k,b≥1} γ_ab (a!/(a−k)!) τ^{a−k+b}.

**Odd-part contributions by order.**
- It first contributes at τ⁵, with 8ψ₃·3!γ₃₁ = 40γ₃₁ = 75ζ(5)/32, which is pure.
- At τ⁶ it gives 8ψ₃[4!γ₄₁ + 3!γ₃₂] = 160γ₄₁ + **40γ₃₂**. Here γ₄₁ = η(6)/32 is pure.
- Hence c₂ … c₅ are pure, and τ⁶ is the first order the symmetric part does not determine.

Symbolic results (`pur_softmax_e_output.txt`):
- c₂ = π²/3, c₃ = 2ζ(3), c₄ = ζ(4)/4, c₅ = 107ζ(5)/32 − ζ(2)ζ(3). These are LAY1's and SM5's values.
- **c₆ = 40γ₃₂ − ζ(3)²/2 + (11/16)ζ(6).** Symbolically identical to SM6's form.
- c₇ = 40γ₃₃ + 160γ₄₂ − ¾ζ(2)²ζ(3) − (5/2)ζ(2)ζ(5) − (7/4)ζ(3)ζ(4) + 7359ζ(7)/512.
  - This equals SM6's 40γ₃₃ + ¾ζ(3)ζ(4) + 2339ζ(7)/512 if and only if 160γ₄₂ is the pure number below.
  - Numerically the difference is 0 at dps 70.
  - PSLQ gives γ₄₂ = (−251ζ(7) + 64ζ(2)ζ(5) + 112ζ(3)ζ(4))/4096. This is NUMERICAL, consistent with the classical reducibility of odd-weight alternating linear Euler sums. I cite that reducibility from memory and have checked it only numerically here.
- Numeric c₆ = 0.0429508278802840339294338048351352655102967277: SM6's Richardson value − this = 2·10⁻⁵¹ (SM6 quoted 50 digits). c₇ agrees with SM6 to 2·10⁻⁴⁷.
- **Own Richardson** of the soft-max F (V = 1/√2, τ = 1/(200·2^i), i ≤ 13, dps 120; `pur_quad_softmax_output.txt`) gives c₆ = 0.0429508278802840339294338048351352655102967277. It agrees with the prediction in all 45 printed digits.

---

## 5. (e) First failure: PROVED (the identity), CONDITIONAL (the non-purity)

**Proposition E.** Suppose ψ_{2i+1} = 0 for i < j and ψ_{2j+1} = ω ≠ 0. Then:
- F_n is pure of weight n for n ≤ 2j+3;
- the next coefficient is

  **F_{2j+4} = [pure part from ψ_even] + ω·(2j+2)!·η(2j+4)/2^{2j+3} + ω·(2j+1)!·γ_{2j+1,2}.**

Here η(m) = (1−2^{1−m})ζ(m). The middle term is the extra pure term that the conjecture's statement omitted; it is harmless. These are per unit of F. For a perimeter I = c₀ + κF, multiply by κ (κ = 8 for the soft-max and for the squircle).

*Proof.* By (B1), odd k contributes τ^{k+1}ψ_k ∂^kG(τ,τ) = ψ_k[k!γ_{k1}τ^{k+2} + ((k+1)!γ_{k+1,1} + k!γ_{k2})τ^{k+3} + O(τ^{k+4})]. Take k = 2j+1, and note that k ≥ 2j+3 starts at τ^{2j+5}. The single moments are elementary: γ_{a1} = (1/a!)Σ_n (−1)^{n+1}(1/n)·a!/(2n)^{a+1} = η(a+2)/2^{a+1}. The interchange is justified by alternating partial sums bounded by x (LAY1 §5). The even-k part is pure by Theorem B. ∎

**The double moment as an Euler sum: PROVED.** We have log²(1+x) = Σ_{n≥2} (−1)^n (2H_{n−1}/n) xⁿ, since Σ_{i+j=n} 1/(ij) = 2H_{n−1}/n. For x ∈ [0,1] the series converges absolutely, as Σ H_{n−1}/n < ∞ times xⁿ. Integrating against s^a e^{−2ns} gives

  **γ_{a,2} = 2^{−(a+1)} Σ_{n≥2} (−1)^n H_{n−1}/n^{a+2} = 2^{−(a+1)} [η(a+3) − Σ_{n≥1}(−1)^{n+1}H_n/n^{a+2}].**

Checked against quadrature for a = 0…6 at dps 80, agreeing to working precision. Thus γ₃₂ = (1/16)Σ(−1)^nH_{n−1}/n⁵, as in SM6.

**j = 0 (weight 4).** The classical evaluation is
  Σ_{n≥1} (−1)^{n+1}H_n/n³ = (11/4)ζ(4) − (7/4)ζ(3)log 2 + ½ζ(2)log²2 − (1/12)log⁴2 − 2Li₄(½).
It appears for instance in Bailey–Borwein–Girgensohn, *Experimental evaluation of Euler sums*, Experiment. Math. 3 (1994), and in Flajolet–Salvy, *Euler sums and contour integral representations*, Experiment. Math. 7 (1998). I did not reopen the sources in this session. **Verified here to 2.9·10⁻¹⁰¹.** Hence

  γ₁₂ = −(15/32)ζ(4) + (7/16)ζ(3)log 2 − ⅛ζ(2)log²2 + (1/48)log⁴2 + ½Li₄(½)  (residual 6·10⁻¹⁰²).

PSLQ on {γ₁₂, ζ(4), ζ(3)log2, ζ(2)log²2, log⁴2, Li₄(½)} independently returns exactly this relation.

**Test of the normalisation at j = 0 (NUMERICAL).** Take ψ = (1+v)e^{−v²}, so ω = 1, V = ∞. Proposition E predicts F₄ = γ₁₂ + π⁴/640, using γ₂₁ = 7ζ(4)/64. Richardson of direct quadrature (τ = 1/(20·2^i), i ≤ 14, dps 120) gives
  F₄ = 0.174148122632792633691086343575542219881455036,
against the predicted 0.17414812263279263369108634357554221988145503563 (`pur_quad_test_e_output.txt`), about 45 digits.

**What is known about non-purity.**
- PROVED (the sums are not *formally* pure): for odd a = 2j+1, the γ_{a,2} are exactly the antisymmetric data that (R) leaves free (SM6 §2 made the same count).
- NUMERICAL: PSLQ at dps 100 (maxcoeff 10¹²) finds no relation between γ_{a,2} and the pure-zeta basis of weight a+3:
  - a = 1: {ζ(4)};
  - a = 3: {ζ(6), ζ(3)²};
  - a = 5: {ζ(8), ζ(3)ζ(5), ζ(2)ζ(3)²}.
  
  For even a = 2, 4 it does find one (γ₂₂ = (16ζ(2)ζ(3) − 29ζ(5))/256, and γ₄₂ above). This matches the parity: odd weight reduces, even weight does not.
- CONDITIONAL: "γ_{2j+1,2} is not a polynomial in ζ-values" follows from the expected ℚ-linear independence of the standard alternating-MZV basis at weight 2j+4. At weight 4 that basis is {ζ(4), ζ(3)log2, ζ(2)log²2, log⁴2, Li₄(½)}, in which γ₁₂ has Li₄(½)-coefficient ½. The independence itself would follow from Deligne's motivic theorem for N = 2 plus the period conjecture. I did not carry out the motivic argument.
- NOT KNOWN (do not claim): irrationality of Li₄(½) or of ζ(5); that γ₁₂/ζ(4) ∉ ℚ; any independence statement among these numbers. The unconditional results nearby are only:
  - ζ(3) ∉ ℚ (Apéry);
  - ζ(2k) ∈ ℚπ^{2k} is transcendental;
  - infinitely many ζ(2k+1) are irrational (Ball–Rivoal);
  - one of ζ(5), ζ(7), ζ(9), ζ(11) is irrational (Zudilin).

---

## 6. (f) The designed square: PROVED (identity, convexity) + NUMERICAL

**Definition.** Take the tent m = cos u + |sin u| at each edge normal jπ/2, |u| ≤ π/4. Let Λ(u) = log(2 + sin u − cos u) and h = m(1 + e^{−2Λ(|u|)/τ})^τ.

**Identity.**
- Λ(0) = 0 and Λ(π/4) = log 2.
- e^Λ Λ′ = cos u + sin u, so (cos u + sin u) du = e^v dv with ψ ≡ 1 and V = log 2.
- ∫₀^{2π} h = 8∫₀^{log 2} e^v(1+e^{−2v/τ})^τ dv = 8 + 8F.
- By §2 (special case), 8 + 8τG(τ,τ) = 16 − 8Γ(1−τ/2)²/Γ(1−τ) = 16 + 2τΓ(−τ/2)²/Γ(−τ).

**Exact exponentially small term.** Let T(V) := ∫_V^∞ e^v[(1+e^{−2v/τ})^τ − 1] dv. Substituting x = e^{−2v/τ} and integrating the binomial series termwise gives

  T(V) = τ e^V Σ_{n≥1} C(τ,n) e^{−2nV/τ}/(2n − τ),

so that

  **∫₀^{2π} h dθ = 16 − 8Γ(1−τ/2)²/Γ(1−τ) − 16τ Σ_{n≥1} C(τ,n) 4^{−n/τ}/(2n−τ).**

The last term is −16τ²·4^{−1/τ}/(2−τ)·(1 + O(4^{−1/τ})), valid for 0 < τ < 2. At τ = 1 the right side is 16 − 0 − 4 = 12, which is exactly ∫h.

**Numerical check** (dps 50, `pur_designed_output.txt`):
- τ = 1, 0.7, 0.5, 0.3, 0.2, 0.1 and 0.05: residual ≤ 4.3·10⁻⁵⁰;
- τ = 0.02: residual 1.4·10⁻⁴⁷;
- a full-circle quadrature without the symmetry reduction, at τ = 0.3, agrees exactly.

**Expansion through τ⁷** (log Γ(1−z) = γz + Σζ(k)z^k/k; pure weight n at τⁿ):

  ∫h = 8 + 2ζ(2)τ² + 2ζ(3)τ³ + (9/8)ζ(4)τ⁴ + (3/2ζ(5) − ½ζ(2)ζ(3))τ⁵ + (79π⁶/120960 − ¼ζ(3)²)τ⁶ + (9/8ζ(7) − (3/8)ζ(2)ζ(5) − (9/32)ζ(3)ζ(4))τ⁷ + O(τ⁸).

Here 79π⁶/120960 = (237/384)ζ(6). In the variables printed by sympy:
- τ⁴: π⁴/80;
- τ⁵: −π²ζ(3)/12 + 3ζ(5)/2;
- τ⁶: −ζ(3)²/4 + 79π⁶/120960;
- τ⁷: −π²ζ(5)/16 − π⁴ζ(3)/320 + 9ζ(7)/8.

As a check, closed form(0.01) − series₇(0.01) = 6·10⁻¹⁸ = O(τ⁸). The first terms c₂ = π²/3 and c₃ = 2ζ(3) coincide with the soft-max's, as they must (same tent, same g₁).

**Convexity: PROVED for every τ ∈ (0,1].** On 0 < u < π/4, write M = cos u + sin u, y = 2Λ/τ, σ = 1/(1+e^y), f = (1+e^{−y})^τ. Using e^ΛΛ′ = M:

  h + h″ = 2σ f M e^{−Λ}·[−3M′ + M²e^{−Λ}(1 + 2K)],  K = (1−σ)/τ + σ ≥ 1.

(This formula matches h + mp.diff(h,2) to 10⁻⁴⁸.) Hence the bracket is at least 3e^{−Λ}(M² − M′e^Λ). Since e^Λ = 2 − M′, M² = 1 + sin 2u and M′² = 1 − sin 2u,

  M² − M′e^Λ = 2(1 − M′) = 2(1 − cos u + sin u) ≥ 0.

So h + h″ ≥ 0 on (0, π/4). At u = 0 the bracket equals 1/τ − 1 ≥ 0.

- **At u = 0:** h′(0⁺) = 2^τ(1 − Λ′(0)) = 0, since Λ′(0) = 1. h is even there, so h ∈ C² across the edge normals (numerically h′(0⁺) ~ 10⁻⁵⁴).
- **At u = π/4:** m′ = 0 and f′ < 0, so h′(π/4⁻) = m f′ < 0. By the reflection u ↦ π/2 − u, the jump is h′(π/4⁺) − h′(π/4⁻) = −2h′(π/4⁻) > 0. It is exponentially small (≈ 3.8·10⁻⁶ at τ = 0.1). A positive jump is a positive point mass of h + h″, which means a short edge with normal π/4 + jπ/2.

Hence h + h″ is a nonnegative measure, so h is the support function of a convex body (Schneider §1.7 / the curvature-measure characterisation). ∫h is its perimeter by Cauchy's formula. Grid minima of h + h″ are ≥ 0 at τ = 1, 0.8, 0.5, 0.2, 0.1 and 0.03, with the minimum 0 attained only at τ = 1, u = 0. No repair was needed.

So: **the designed body is a convex curve whose perimeter is, to all orders, 16 − 8Γ(1−τ/2)²/Γ(1−τ). Every coefficient is pure, by construction, and the remainder is exponentially small and explicit.**

---

## 7. (g) Scope

**Does the criterion explain the squircle's pure-weight theorem?** It gives an independent complete proof, to every finite order, of thm:pw (§3). The mechanism is, however, the same as the paper's:
- The paper writes P as a full Beta integral ∫₀¹(t(1−t))^{η/2} √cosh(ηl) dt, with t ↔ 1−t symmetry and an even weight √cosh.
- (R) is the half-line version of that Beta integral.
- v ↦ −v is t ↔ 1−t.
- e^{−v}φ even ⟺ "the weight is even under t ↔ 1−t".

So for the squircle it **reproduces** the theorem rather than explaining something new. What it adds:
1. It isolates the single condition that matters. The body need not have a global Beta representation; any even ψ works (§6).
2. It locates the first failure and its constant exactly (§5).
3. It handles finite V, with exponentially small corrections.

It does not address the paper's n-dimensional surface-area extension (not examined).

**Rhombi.**

*(i) p-norm rhombi and every linear image of B_p: PROVED pure.*
- Let K = A·B_p, so h(n) = ‖Aᵀn‖_q. Put (x,y) = Aᵀn and t = y/x.
- Then x dy − y dx = det A·dθ, and |x| = Q(1,t)^{−1/2} with Q(ξ) = |A^{−T}ξ|². Hence m dθ = |x|³dt/det A.
- At a kink |x| = |y|, set t = ±e^{−2v} on one side and the reciprocal on the other. Homogeneity of Q gives the two-sided weight

  φ(v) = e^v[q(v) + q(−v)],  q(v) = 2/(det A·Q(e^v, ±e^{−v})^{3/2}).

- That is e^v times an **even** function. Each half-window runs over v ∈ [0, ∞) with no remainder.
- So all coefficients are pure of weight n, with coefficients in ℚ(entries)·ζ. This is consistent with the paper's "every rhombus".
- The two sides of a kink are individually *not* of the form e^v·even (the evenness is about v₀ = ½log(a/b)). Only their sum is.

*Check for diag(1,2), r = √5* (`pur_softmax_e_output.txt`, `pur_quad_rhombus_output.txt`):
- c₂ = 16√5ζ(2)/25 and c₃ = 16√5ζ(3)/25, equal to the paper's 4A²B²ζ(k)/r³ (symbolic difference 0);
- c₄ = 16√5(4ζ(4) − ζ(2)²)/125;
- c₅ = 16√5(3ζ(5) − 2ζ(2)ζ(3))/125;
- the exact v-form equals θ-quadrature to ≤ 2·10⁻⁵⁰;
- Richardson confirms c₄ and c₅ to 2.6·10⁻⁴⁰.

*(ii) What replaces e^v for a tent H cos u + T|sin u| (S = 0, e.g. a rectangle).* Nothing. The fixed line λ = μ of the reflection is universal, so the weight must always be e^v × (even). What changes is the normalisation of the layer coordinate:
- v must be chosen so that the layer reads (1+e^{−2v/τ})^τ;
- C¹ of h across the kink forces Λ′(0) = T/H, since h′(0±) = ±2^τ(T − HΛ′(0)).

A designed choice that works for every such tent is

  e^{Λ(u)} = 1 + (T/H) sin u + (T/H)²(1 − cos u),  (H cos u + T sin u) du = (H²/T) e^v dv,  ψ ≡ H²/T.

The convexity proof of §6 carries over:
- M² − c e^Λ M′ vanishes at 0 and has derivative M(c+T) > 0;
- at u = 0 the bracket is T(1/τ − 1) ≥ 0.

Checked for the rectangle A = 1, B = 2, splitting the quarter at the angle where the two Λ's agree:
- ∫h = 4(A+B) + 4(A²/B + B²/A)[1 − Γ(1−τ/2)²/Γ(1−τ) − T(V)] to ≤ 4.3·10⁻⁵⁰ at τ = 1, ½, 0.2, 0.1, 0.05;
- h + h″ > 0 on the grid, the jumps at the split are positive, and h′ = 0 at the edge normals (`pur_designed_output.txt`).

*(iii) The LAY1-type soft-max of a polygon is generally worse than the square's.*
- Near edge j the layer is e^{−2T sin u/τ}, so v = T sin u.
- The two sides give φ = 2H/T + 2 tan u. Here S cancels between the sides.
- e^{−v}φ has v¹-coefficient 2(1−H)/T. So **the soft-max is pure through τ³ always, and its τ⁴ coefficient contains γ₁₂ (hence Li₄(½), conditionally) unless every H_j = 1.**
- When H = 1, the v³-coefficient per side is 1/(3T) + 1/(2T³) > 0, so the τ⁶ failure is unavoidable.

The unit square (H = T = 1) is the case SM5/SM6 saw. This item rests on (a) plus the soft-max window reduction of SM6 §2. That reduction was not re-derived for general polygons: it is the same exponentially-small bookkeeping, so status PROVED modulo that sketch.

---

## 8. Numerics table

| quantity | method | agreement |
|---|---|---|
| (a) vs (b), squircle c₁…c₇ | γ_ab quadrature, dps 70 | ≤ 2.9·10⁻⁷⁰ |
| squircle c₁…c₇ vs paper / SM6 values | symbolic | exact (difference 0) |
| squircle P: v-form vs θ-quadrature | dps 50, 4 values of τ | 1.1·10⁻⁵⁰ |
| squircle c₆, c₇ | Richardson, dps 120 | 7.8·10⁻³⁸, 2.2·10⁻³⁶ |
| soft-max c₂…c₇: (a) vs criterion | dps 70 | ≤ 1.5·10⁻⁷⁰ |
| soft-max c₆ = 11/16ζ(6) − ζ(3)²/2 + 40γ₃₂ | symbolic vs SM6 | exact |
| soft-max c₆ | own Richardson, dps 120 | 45 digits |
| soft-max c₆, c₇ vs SM6 Richardson | numeric | 2·10⁻⁵¹, 2·10⁻⁴⁷ |
| γ_{a,2} = Euler sum, a = 0…6 | quad vs nsum, dps 80 | working precision |
| Σ(−1)^{n+1}H_n/n³ closed form | dps 100 | 2.9·10⁻¹⁰¹ |
| test ψ = (1+v)e^{−v²}: F₄ = γ₁₂ + π⁴/640 | Richardson, dps 120 | ~45 digits |
| designed square ∫h = closed form + E | dps 50, 8 values of τ | ≤ 1.4·10⁻⁴⁷ |
| designed square h + h″ formula | vs mp.diff | ≤ 8·10⁻⁴⁹ |
| designed rectangle ∫h | dps 50, 5 values of τ | ≤ 4.3·10⁻⁵⁰ |
| p-norm rhombus diag(1,2): v-form vs θ | dps 50 | ≤ 2·10⁻⁵⁰ |
| p-norm rhombus c₄, c₅ | Richardson, dps 120 | 2.6·10⁻⁴⁰ |

## 9. Status

| item | status |
|---|---|
| (a) expansion F_n = Σφ_k k!γ_{k,n−1−k}, remainder O(τ^{N+1}) + O(e^{−V/τ}) | **PROVED** |
| (b) F ~ τΣψ_kτ^k∂^kG(τ,τ); even derivatives = ½∂^{2k}S on λ = μ; pure weight n if ψ odd-free and ψ_{2k} ∈ ℚ | **PROVED** |
| (c) squircle is exactly of the form, φ = e^v/(√2cosh^{3/2}2v), V = ∞ | **PROVED** (+ 50-digit check) |
| (c) squircle c₂…c₇ from (b), equal to the known values | **PROVED** (symbolic) + **NUMERICAL** (≥ 36 digits) |
| (c) squircle pure weight to every order | **PROVED** (reproof of thm:pw) |
| (d) soft-max ψ₁ = 0, ψ₃ = 5/6; c₂…c₅ pure; c₆ = pure + 40γ₃₂ = SM6 | **PROVED** (given SM6's exponentially-small reduction) + **NUMERICAL** |
| (d) c₇ = SM6's form | **NUMERICAL** (needs γ₄₂ pure: PSLQ) |
| (e) through τ^{2j+3} pure; F_{2j+4} = pure + ω(2j+2)!η(2j+4)/2^{2j+3} + ω(2j+1)!γ_{2j+1,2} | **PROVED** (the conjecture missed the middle pure term) |
| (e) γ_{a,2} = 2^{−a−1}Σ(−1)^nH_{n−1}/n^{a+2} | **PROVED** + NUMERICAL |
| (e) γ₁₂ closed form with Li₄(½), log 2 | **NUMERICAL** (classical formula, cited; verified to 100 digits) |
| (e) γ_{2j+1,2} not a polynomial in ζ-values | **CONDITIONAL** (independence of alternating-MZV bases; numerical: no PSLQ relation for j = 0, 1, 2); unconditionally **OPEN** |
| (f) ∫h = 16 − 8Γ(1−τ/2)²/Γ(1−τ) − 16τΣC(τ,n)4^{−n/τ}/(2n−τ) | **PROVED** + NUMERICAL (≥ 47 digits) |
| (f) expansion through τ⁷ | **PROVED** (symbolic) |
| (f) h is a support function (convex), all τ ∈ (0,1] | **PROVED** + grid check |
| (g) linear images of B_p (all rhombi): pure to every order | **PROVED** (+ numerical for diag(1,2)) |
| (g) designed tents H cos u + T\|sin u\|: weight (H²/T)e^v; convex; closed form | **PROVED** (argument) + NUMERICAL (rectangle 1×2) |
| (g) polygon soft-max: τ⁴ impure unless all H_j = 1; τ⁶ always | **PROVED** modulo the window sketch; impurity itself **CONDITIONAL** |
| criterion "explains" thm:pw | **No**: same Beta mechanism; it reproduces it, and adds the condition, the designed examples and the failure law |

## 10. Five-sentence summary for a referee

1. For F(τ) = ∫₀^V φ(v)[(1+e^{−2v/τ})^τ − 1]dv, substituting v = τs and expanding gives F_n = Σ_k φ_k k! γ_{k,n−1−k} with an O(τ^{N+1}) remainder. Rewriting φ = e^vψ turns this into F ~ τΣψ_kτ^k∂_λ^kG(λ,τ)|_{λ=τ}.
2. Since λ = μ is the fixed line of λ ↦ 2μ−λ, Lemma 5 gives every even derivative there as ½∂^{2k} of the Gamma-ratio S. Therefore an even ψ yields coefficients that are rational polynomials in ζ(k) of weight n, with no alternating constants.
3. The squircle is exactly such an integral, with ψ = 1/(√2cosh^{3/2}2v). The criterion reproduces c₂…c₇ symbolically (c₆, c₇ also by 36-digit Richardson), and the same computation proves every linear image of the p-ball pure; the mechanism is the paper's Beta symmetry t ↔ 1−t, not a new one.
4. If e^{−v}φ has odd part ωv^{2j+1} + …, then F_{2j+4} = pure + ω(2j+2)!η(2j+4)/2^{2j+3} + ω(2j+1)!γ_{2j+1,2}, with γ_{a,2} = 2^{−a−1}Σ(−1)^nH_{n−1}/n^{a+2}. For the soft-max (ω = 5/6, j = 1) this gives exactly SM6's 40γ₃₂ at τ⁶, and γ₁₂ contains Li₄(½); that these sums are not pure is supported numerically and conditional on standard independence conjectures.
5. The designed square h = m(1+e^{−2Λ/τ})^τ with Λ = log(2 + sin u − cos u) is provably convex for τ ∈ (0,1]. Its perimeter is 16 − 8Γ(1−τ/2)²/Γ(1−τ) − 16τΣ_{n≥1}C(τ,n)4^{−n/τ}/(2n−τ), checked to 47–50 digits.

---

## Planner notes after the referee (2026-11-22)

Refereed by an independent Opus agent with its own code (`REFEREE.md`, `referee/`). (a), (b), (c), (d), (f) hold. Read these fixes into the text:
1. (a): for finite V assume φ bounded (or L¹), not "automatic".
2. (e): justify the Euler-sum step by Tonelli (Σ H_{n−1}/n diverges). Non-purity is CONDITIONAL everywhere; "not formally pure" in §5,
   "τ⁶ unavoidable" in §7(iii) and "γ₁₂ contains Li₄(½)" are remarks, not theorems.
3. (g)(i): the coefficient field for a linear image of B_p is ℚ(entries, √Q(1,±1)), not ℚ(entries).
4. (g)(ii): the jump sign at the split has a one-line proof in REFEREE.md.
5. (g)(iii) is FALSE as stated: the γ₁₂ coefficient of a polygon's soft-max is Σ_j 2(1−H_j)/T_j, which can vanish with H_j ≠ 1 (the
   rectangle 6/5 × 3/5 has c₄ = (5/16)ζ(4)). Do not use.
For the paper: the criterion, the squircle and soft-max cases, and the designed curve (closed form to ≥ 56 digits; convexity an inequality
chain for all τ in (0,1]; jump of h′ at π/4 equal to 4(1+4^{−1/τ})^τ/(1+4^{1/τ}) > 0).
