# LAY1: the Layer Lemma with a proved remainder

**Outcome.** The Layer Lemma is now a theorem (Theorem 1 below), with remainder O(τ^{N+1}) to every order N. Its hypotheses hold, with proofs, for the squircle at its diamond end and for the soft-max of the square, at N = 4. There is no gluing step to estimate: subtract m·A (not the polygon's m) on each window, and the window integrand is exactly m·A·(e^L − 1), which decays. The s-independent parts are then never separated from the bulk integral ∫ m·A, so there is nothing to recombine. Two new facts came out of the work:
1. the soft-max is convex for every τ in (0, 1] (proved);
2. the moments behind c₃ and c₄ follow from a single Beta-function reflection identity, which also proves AV5b's numerical J = ζ(4)/32 and gives c₄ for both examples.

Notation: T = ℝ/2πℤ, τ ∈ (0, τ₀], and I(τ) := ∫₀^{2π} h(θ; τ) dθ.

---

## 1. Theorem

**Data.**

(H0) **Polygon.** m is the support function of a convex polygon K₀. Its edge normals are θ₁, …, θ_J. Fix W > 0 smaller than half the least angular gap between consecutive θ_j, and set I_j = [θ_j − W, θ_j + W]. The windows are disjoint, and on I_j only the two vertices of edge j matter, so for |u| ≤ W
  m(θ_j + u) = H_j cos u + S_j sin u + T_j |sin u|,
where H_j = m(θ_j), T_j = ℓ_j/2 ≥ 0 (ℓ_j is the edge length), and S_j is the tangential coordinate of the edge midpoint. *Derivation:* for u > 0, m = v₊·n(θ_j+u) = H_j cos u + t₊ sin u; for u < 0, m = H_j cos u + t₋ sin u; t₊ ≥ t₋; then S = (t₊+t₋)/2 and T = (t₊−t₋)/2. The paper's β is T_j/H_j, and the paper's form assumes S_j = 0.

(H1) **Bulk factor.** A(θ, τ) > 0. There are bounded measurable a₁, …, a_N on T and C₁ with
  |log A(θ,τ) − Σ_{k=1}^N τ^k a_k(θ)| ≤ C₁ τ^{N+1} for all θ, τ.
In addition, on each window, a_k ∈ C^{N−k}(I_j) for 1 ≤ k ≤ N−1. For N = 3 this means a₁ ∈ C², a₂ ∈ C¹.

(H2) **Layer.** For each j there are measurable g_{j,1}, …, g_{j,N−1}: ℝ → ℝ and constants C₂, c > 0 such that for |s| ≤ W/τ
  h(θ_j + τs; τ) = m(θ_j+τs) A(θ_j+τs, τ) exp(L_j(s,τ)),  L_j = Σ_{k=1}^{N−1} τ^k g_{j,k}(s) + R_j(s,τ),
with |g_{j,k}(s)| ≤ C₂ e^{−c|s|} for all s ∈ ℝ, and |R_j(s,τ)| ≤ C₂ τ^N e^{−c|s|} for |s| ≤ W/τ.
- A polynomial weight (1+|s|)^M e^{−c|s|} is allowed; replace c by c/2.
- An error K e^{−κ/τ} that does not decay in s is allowed: on |s| ≤ W/τ it is ≤ K e^{−(κ−c′W)/τ} e^{−c′|s|} for any c′ < κ/W.

(H3) **Outside the windows.** |h(θ;τ) − m(θ)A(θ,τ)| ≤ C₃ e^{−κ/τ} for θ ∉ ∪_j I_j, and h(·;τ) is integrable.

**Coefficient polynomials.** Everything below is defined by generating functions in an auxiliary variable t.
- B_n(θ) is the tᵅ coefficient of exp(Σ_{k=1}^N t^k a_k(θ)), with n in place of α. So B₀ = 1, B₁ = a₁, B₂ = a₂ + a₁²/2, B₃ = a₃ + a₁a₂ + a₁³/6.
- M_{j,a}(s) is the t^a coefficient of H_j cos(ts) + S_j sin(ts) + T_j|sin(ts)|:
  - M₀ = H
  - M₁ = Ss + T|s|
  - M₂ = −Hs²/2
  - M₃ = −(Ss + T|s|)s²/6
- P_{j,b}(s) is the t^b coefficient of exp(Σ_{k=1}^{N−1} t^k Σ_{i=0}^{N−1−k} a_k^{(i)}(θ_j) (ts)^i/i!):
  - P₀ = 1
  - P₁ = a₁(θ_j)
  - P₂ = a₂(θ_j) + a₁(θ_j)²/2 + a₁′(θ_j) s
- G_{j,n}(s) is the tⁿ coefficient of exp(Σ_{k=1}^{N−1} t^k g_{j,k}(s)) − 1:
  - G₁ = g₁
  - G₂ = g₂ + g₁²/2
  - G₃ = g₃ + g₁g₂ + g₁³/6
- Φ_{j,n}(s) = Σ_{a+b+c=n, c≥1} M_{j,a}(s) P_{j,b}(s) G_{j,c}(s) for 1 ≤ n ≤ N−1, and Φ_{j,0} = 0.

**Theorem 1 (Layer Theorem).** Under (H0)–(H3), there are τ₁ > 0 and C such that for 0 < τ ≤ τ₁
  | I(τ) − Σ_{n=0}^N c_n τ^n | ≤ C τ^{N+1},  c_n = ∫_T m B_n dθ + Σ_j ∫_ℝ Φ_{j,n−1}(s) ds,
and every integral converges absolutely. Explicitly:
- c₀ = ∫ m dθ = perimeter of K₀
- c₁ = ∫ m a₁ dθ
- c₂ = ∫ m (a₂ + a₁²/2) dθ + Σ_j H_j ∫_ℝ g_{j,1} ds
- **c₃ = ∫ m (a₃ + a₁a₂ + a₁³/6) dθ + Σ_j ∫_ℝ [ H_j (g₂ + g₁²/2) + (S_j s + T_j|s| + H_j a₁(θ_j)) g₁ ] ds**
- c₄ = ∫ m B₄ dθ + Σ_j ∫_ℝ [ H G₃ + (Ss+T|s|+H P₁) G₂ + (−Hs²/2 + (Ss+T|s|)P₁ + H P₂) G₁ ] ds

For N = 3, with S_j = 0 (or with even g₁) and β_j = T_j/H_j, this is exactly the paper's statement. Each kink adds H∫g₁ to c₂ and H∫[g₂ + g₁²/2 + (β|s| + a₁)g₁] to c₃, and the remainder is O(τ⁴). The constant C depends only on N, C₁, C₂, C₃, c, κ, W, max|H_j|+|S_j|+T_j, and sup|a_k^{(i)}| on the windows. The proof gives it explicitly.

**Where convexity enters.** Theorem 1 is a statement about the integral I(τ) and uses no convexity. Convexity is used once, to read I(τ) as a perimeter. If h(·;τ) is the support function of a convex body K_τ, then I(τ) = per(K_τ) by Cauchy's formula. For C² functions h, being a support function is equivalent to h + h″ ≥ 0, and then per = ∫(h + h″) = ∫h. Both examples are support functions: the squircle classically, the soft-max by Proposition 4.

---

## 2. Proof of Theorem 1

Take τ ≤ τ₁ := min(τ₀, 1). Write C for constants of the kind listed above; their values change from line to line.

### Lemma 1 (exponential truncation with a majorant)

Let y₁, …, y_K ∈ ℂ, ω ≥ 1, Ω > 0, Y₀ ≥ 0 with |y_k| ≤ Y₀ ω^k. Define E_n by exp(Σ_k t^k y_k) = Σ_{n≥0} E_n tⁿ, and put S_Ω = Σ_{k=1}^K Ω^k and Z = Y₀ S_Ω e^{Y₀ S_Ω}. Then:
- |E_n| ≤ (ω/Ω)ⁿ Z for n ≥ 1;
- for M ≥ 1 and 0 ≤ τω ≤ Ω, |exp(Σ_k τ^k y_k) − Σ_{n<M} τⁿ E_n| ≤ (τω/Ω)^M Z.

*Proof.* Let Ẽ_n be the coefficients of exp(Σ t^k |y_k|) − 1. These are nonnegative, and |E_n| ≤ Ẽ_n for n ≥ 1, because every monomial in E_n has a nonnegative rational coefficient. Put r = Ω/ω. Then
  Σ_{n≥1} rⁿ Ẽ_n = e^{X} − 1, with X = Σ r^k |y_k| ≤ Y₀ Σ Ω^k = Y₀ S_Ω.
Since e^X − 1 ≤ X e^X, this sum is ≤ Z. Hence Ẽ_n ≤ Z r^{−n}, which is the first claim. For M ≥ 1 and τ ≤ r,
  |Σ_{n≥M} τⁿ E_n| ≤ (τ/r)^M Σ_{n≥M} rⁿ Ẽ_n ≤ (τ/r)^M Z,
and τ/r = τω/Ω. ∎

### Step 0: exact splitting (no gluing)

h and mA are integrable and the windows are disjoint, so
  I(τ) = ∫_T mA dθ + Σ_j ∫_{I_j} (h − mA) dθ + ∫_{T∖∪I_j} (h − mA) dθ.  (2.1)
By (H2), h − mA = m A (e^{L_j} − 1) on I_j. Substituting θ = θ_j + τs (dθ = τ ds) gives
  ∫_{I_j} (h − mA) dθ = τ ∫_{−W/τ}^{W/τ} Ψ_j(s,τ) ds,  Ψ_j := m(θ_j+τs) A(θ_j+τs,τ) (e^{L_j(s,τ)} − 1).  (2.2)
This identity is exact. Nothing has been expanded, and no window-edge matching is needed, because the bulk integral runs over all of T with the true m·A.

### Step 1: bulk

Apply Lemma 1 with y_k = a_k(θ), ω = Ω = 1, Y₀ = max_k sup|a_k|, M = N+1:
  |exp(Σ_{k≤N} τ^k a_k) − Σ_{n≤N} τⁿ B_n| ≤ τ^{N+1} Z.
By (H1), A = exp(Σ_{k≤N} τ^k a_k) · e^{ρ} with |ρ| ≤ C₁τ^{N+1} ≤ C₁. Also |e^ρ − 1| ≤ |ρ| e^{|ρ|} and |exp(Σ_{k≤N} τ^k a_k)| ≤ e^{Y₀N}. Therefore
  |A − Σ_{n≤N} τⁿ B_n| ≤ τ^{N+1} (Z + e^{Y₀N} C₁ e^{C₁})
uniformly in θ. Multiply by |m| ≤ max m and integrate over T:
  ∫_T mA dθ = Σ_{n=0}^N τⁿ ∫_T m B_n dθ + O(τ^{N+1}).  (2.3)

### Step 2: the three factors of Ψ_j on |s| ≤ W/τ

Note τ|s| ≤ W throughout.

**(a) The tent.** On each side of s = 0, the function τ ↦ m(θ_j+τs) is H cos(τs) + (S ± T) sin(τs). Taylor's theorem with Lagrange remainder, using |cos^{(N)}|, |sin^{(N)}| ≤ 1, gives
  m(θ_j+τs) = Σ_{a<N} τ^a M_a(s) + r_m,  |r_m| ≤ (|H|+|S|+T)(τ|s|)^N / N!.
Also |M_a(s)| ≤ (|H|+|S|+T)|s|^a/a!. So the truncation m_< := Σ_{a<N} τ^a M_a satisfies
  |m_<| ≤ (|H|+|S|+T) Σ_a (τ|s|)^a/a! ≤ (|H|+|S|+T) e^W.

**(b) The bulk factor on the window.** By (H1), and since |a_N| is bounded, log A(θ_j+τs, τ) = Σ_{k=1}^{N−1} τ^k a_k(θ_j+τs) + ρ₁ with |ρ₁| ≤ C τ^N. Because a_k ∈ C^{N−k}(I_j), Taylor's theorem gives
  a_k(θ_j+τs) = Σ_{i=0}^{N−1−k} a_k^{(i)}(θ_j)(τs)^i/i! + ρ_{k}',  |ρ_k'| ≤ sup_{I_j}|a_k^{(N−k)}| (τ|s|)^{N−k}/(N−k)!.
Hence log A = Σ_{n=1}^{N−1} τⁿ λ_n(s) + ρ₂, where
- λ_n(s) = Σ_{k+i=n} a_k^{(i)}(θ_j) s^i/i!, so |λ_n(s)| ≤ Y₀(1+|s|)ⁿ;
- |ρ₂| ≤ C τ^N (1+|s|)^N.

Apply Lemma 1 with ω = 1+|s| and Ω = 1+W. This is allowed because τω ≤ τ + W ≤ Ω. With M = N:
  |exp(Σ τⁿλ_n) − Σ_{n<N} τⁿ P_n(s)| ≤ Z (τ(1+|s|))^N/(1+W)^N,  |P_n(s)| ≤ Z (1+|s|)ⁿ.
Also |exp(Σ τⁿλ_n)| ≤ e^{Y₀ S_Ω}. Since |ρ₂| ≤ C(1+W)^N, we have |e^{ρ₂} − 1| ≤ C τ^N (1+|s|)^N. Therefore
- |A| ≤ C;
- |r_A| := |A − Σ_{n<N} τⁿ P_n| ≤ C τ^N (1+|s|)^N;
- |A_<| := |Σ_{n<N} τⁿ P_n| ≤ Z Σ (τ(1+|s|))ⁿ ≤ C.

**(c) The layer.** Apply Lemma 1 with y_k = g_k(s), ω = Ω = 1, Y₀ = C₂ e^{−c|s|} ≤ C₂, so Z ≤ C₂ S₁ e^{C₂ S₁} e^{−c|s|}:
- |G_n(s)| ≤ C e^{−c|s|} for all s ∈ ℝ;
- |exp(Σ τ^k g_k) − 1 − Σ_{n=1}^{N−1} τⁿ G_n| ≤ C τ^N e^{−c|s|}, by M = N;
- |exp(Σ τ^k g_k) − 1| ≤ C τ e^{−c|s|}, by M = 1.

With e^{L} = exp(Σ τ^k g_k) e^{R}, where |R| ≤ C₂ τ^N e^{−c|s|} ≤ C₂ so that |e^R − 1| ≤ C₂ e^{C₂} τ^N e^{−c|s|}, we get
- |E| := |e^{L} − 1| ≤ C e^{−c|s|};
- |r_E| := |e^{L} − 1 − E_<| ≤ C τ^N e^{−c|s|}, where E_< := Σ_{n=1}^{N−1} τⁿ G_n.

### Step 3: expansion of Ψ_j

  Ψ_j − m_< A_< E_< = r_m · A · E + m_< · r_A · E + m_< · A_< · r_E.
By Step 2, each of the three terms is bounded by C τ^N (1+|s|)^N e^{−c|s|}.

Expanding the product gives m_< A_< E_< = Σ_{n=1}^{3N−3} τⁿ Σ_{a+b+c=n; a,b≤N−1; 1≤c≤N−1} M_a P_b G_c. For n ≤ N−1 the constraints a, b, c ≤ N−1 hold automatically, so the inner sum is Φ_{j,n}. For n ≥ N each term is bounded by τⁿ C|s|^a (1+|s|)^b e^{−c|s|} ≤ C τ^N (1+|s|)^{2N−2} e^{−c|s|} (using τ ≤ 1). So
  |Ψ_j(s,τ) − Σ_{n=1}^{N−1} τⁿ Φ_{j,n}(s)| ≤ C τ^N (1+|s|)^{2N} e^{−c|s|}  for |s| ≤ W/τ.  (2.4)

### Step 4: integrate (2.4) and remove the window cut-off

Let K_N = sup_{x≥0} (1+x)^{2N} e^{−cx/2} < ∞. Then:
- Integrating (2.4) over |s| ≤ W/τ: τ ∫ |Ψ_j − Σ τⁿΦ_{j,n}| ds ≤ C τ^{N+1} ∫_ℝ (1+|s|)^{2N} e^{−c|s|} ds ≤ C τ^{N+1} · 4K_N/c.
- Each Φ_{j,n} is bounded by C(1+|s|)^{2N} e^{−c|s|} on all of ℝ (Step 2, where the bounds on M_a, P_b, G_c hold for every s). So it is absolutely integrable, and the cut-off tail satisfies
  |∫_{|s|>W/τ} Φ_{j,n}| ≤ C K_N ∫_{|s|>W/τ} e^{−c|s|/2} ds = (4CK_N/c) e^{−cW/(2τ)}.
- Since e^{x} ≥ x^{N+1}/(N+1)!, any e^{−κ′/τ} is ≤ (N+1)! τ^{N+1}/κ′^{N+1}.

Therefore
  ∫_{I_j} (h − mA) dθ = Σ_{n=1}^{N−1} τ^{n+1} ∫_ℝ Φ_{j,n} ds + O(τ^{N+1}).  (2.5)

### Step 5: the outside region and assembly

By (H3), |∫_{T∖∪I_j} (h − mA)| ≤ 2π C₃ e^{−κ/τ} = O(τ^{N+1}). Put (2.3), (2.5) and this bound into (2.1) and collect powers of τ; the window term contributes to c_n through Φ_{j,n−1}. This proves Theorem 1. The explicit c₂, c₃, c₄ follow by reading off Φ₁ = H G₁, Φ₂ = H G₂ + (M₁ + H P₁) G₁, and Φ₃ = H G₃ + (M₁ + HP₁) G₂ + (M₂ + M₁P₁ + HP₂) G₁. ∎

**Remarks.**
1. *Comparison with the formal derivation.* The formal derivation expands m·A·e^L on the window and then has to glue its s-independent part back to the bulk; that is the unestimated step. Step 0 removes it: the s-independent part is m·A itself, integrated over all of T, and only m·A·(e^L − 1) is localized.
2. *Profiles with a nonzero limit.* (H2) requires L_j → 0 as |s| grows. A nonzero limit, such as AV5b's hexagon profile g₁ = ln 2 + ½log(1+e^{−4α|s|}), must be moved into A. This works only if the limits on the two sides agree and the result satisfies (H1). The hexagon is not one of the two examples here, so this is not checked here.
3. *No logarithms.* Under (H0)–(H3) to order N, the expansion has no τ^k log τ terms. AV5b found d^k log d terms in its frozen-rest local models for m ≥ 4. If those terms are real, those models violate (H1)/(H2) at the relevant order. Not examined here.

---

## 3. Example (i): the squircle at the diamond end

Setup: p > 1, q = p/(p−1), τ = 1 − 1/p = 1/q ∈ (0, 1).

**Support function.** By Hölder, sup_{‖x‖_p≤1} x·n = ‖n‖_q, with equality at x_i = sgn(n_i)|n_i|^{q−1}/‖n‖_q^{q−1}. So
  h_p(θ) = (|cos θ|^q + |sin θ|^q)^{1/q},
which is the support function of the convex body B_p. Hence I(τ) = per(B_p) by Cauchy's formula.

**Polygon.** m(θ) = max(|cos θ|, |sin θ|) is the support function of the diamond conv{(±1,0), (0,±1)}, with per = 4√2. Its edge normals are θ_j = π/4 + jπ/2. All four kinks are alike under the dihedral symmetry D₄ of h_p and m, so take θ₀ = π/4.

Write r = min/max of (|cos θ|, |sin θ|). Then
  h_p = m (1 + r^q)^{1/q} = m · exp(τ log(1 + r^{1/τ})),  A ≡ 1, so all a_k = 0.

**Tent.** cos(π/4+u) = (cos u − sin u)/√2 and sin(π/4+u) = (cos u + sin u)/√2, so
  m(π/4+u) = (cos u + |sin u|)/√2,  giving H = T = 1/√2, S = 0, β = 1,
and r = (cos u − |sin u|)/(cos u + |sin u|) = tan(π/4 − |u|).

**Window: W = π/8.** Put x = |u| = τ|s| ≤ π/8 and Λ(x) := artanh(tan x). Since tan(π/4−x) = (1−tan x)/(1+tan x), we have log r = −2Λ(x) and
  L(s,τ) = τ F(2Λ(x)/τ),  F(y) := log(1 + e^{−y}).
Λ′(x) = sec²x/(1 − tan²x) = 1/cos 2x = sec 2x, and Λ(0) = 0. On y ∈ [0, π/4], sec y − 1 = (1 − cos y)/cos y ≤ (y²/2)/(1/√2) = y²/√2. Hence
  0 ≤ Λ(x) − x = ∫₀^x (sec 2t − 1) dt ≤ (2√2/3) x³  (0 ≤ x ≤ π/8).
So 2Λ/τ = 2|s| + δ with 0 ≤ δ ≤ (4√2/3) τ²|s|³. The derivative satisfies |F′(y)| = 1/(1+e^y) ≤ e^{−y}. By the mean value theorem on [2|s|, 2|s|+δ],
  |L − τ g₁(s)| ≤ τ δ e^{−2|s|} ≤ (4√2/3) τ³ |s|³ e^{−2|s|} ≤ (4√2/3)(27/e³) τ³ e^{−|s|},  g₁(s) = log(1+e^{−2|s|}),
using |s|³ e^{−|s|} ≤ 27e^{−3}. This gives (H2) with N = 3, g₂ = 0 and c = 1. The decay bound on g₁ itself holds because log(1+y) ≤ y, so 0 ≤ g₁ ≤ e^{−2|s|}.

**Order N = 4.** sec has nonnegative Taylor coefficients, so (sec y − 1 − y²/2)/y⁴ is increasing on [0, π/2). Hence for y ≤ π/4,
  sec y − 1 − y²/2 ≤ K y⁴,  K = (√2 − 1 − π²/32)/(π/4)⁴ < 0.28.
So 2Λ/τ = 2|s| + (4/3)τ²|s|³ + δ′ with 0 ≤ δ′ ≤ (32K/5) τ⁴|s|⁵. Write δ_tot := (4/3)τ²|s|³ + δ′ for the total shift. Second-order Taylor of F, with |F″(y)| = e^y/(1+e^y)² ≤ e^{−y} for y ≥ 0, gives
  L = τ g₁ + τ³ g₃ + R,  g₃(s) = −(4/3)|s|³/(1+e^{2|s|}),
  |R| ≤ τ e^{−2|s|} [(32K/5)τ⁴|s|⁵ + ½ δ_tot²] ≤ C τ⁵ (1+|s|)⁶ e^{−2|s|},
using δ_tot ≤ (4√2/3)τ²|s|³. So (H2) holds with N = 4, g₂ = 0, and g₃ as above.

**Outside the windows.** The set |θ − θ_j| ≥ π/8 for all j consists, by D₄ symmetry, of copies of |θ| ≤ π/8. There r = |tan θ| ≤ tan(π/8) = √2 − 1, and (1+y)^τ − 1 ≤ τy for τ ∈ (0,1] (Bernoulli). Therefore
  0 ≤ h_p − m ≤ τ (√2−1)^{1/τ} = τ e^{−log(1+√2)/τ},
which is (H3) with κ = log(1+√2).

**Status: (H0)–(H3) PROVED for N = 3 and N = 4.**

---

## 4. Example (ii): the soft-max of the square

h(θ) = m [Σ_v e^{(v·n − m)/τ}]^τ, where v ranges over (±1, ±1), n = (cos θ, sin θ), and m = |cos θ| + |sin θ| is the support function of [−1,1]², with per = 8. Edge normals θ_j = jπ/2; H = 1, T = 1 (edge length 2), S = 0, β = 1. A ≡ 1.

**Window: W = π/4.** The four windows tile T, so (H3) is vacuous. For 0 ≤ u ≤ π/4 the maximizing vertex is (1,1), and the other three exponents v·n − m are −2 sin u, −2 cos u and −2(cos u + sin u). The case −π/4 ≤ u ≤ 0 is the mirror image, with the same values in |u|. Exactly, then,
  L(s,τ) = τ log(1 + e^{−2 sin(τ|s|)/τ} + b),  b = e^{−2 cos u/τ}(1 + e^{−2|sin u|/τ}) ∈ (0, 2e^{−√2/τ}].

Put x = τ|s| ≤ π/4. Then 2 sin x/τ = 2|s| − δ with 0 ≤ δ = 2(x − sin x)/τ ≤ x³/(3τ) = τ²|s|³/3. Also δ ≤ |s| x²/3 ≤ 0.206|s|, and 2 sin x/τ ≥ 0. Three facts give the bound:
- the mean value theorem with |F′(y)| ≤ e^{−y} on [2|s| − δ, 2|s|];
- |log(1+a+b) − log(1+a)| ≤ b for a, b ≥ 0;
- the (H2) remark on non-decaying e^{−κ/τ} errors.
Together:
  |L − τ g₁(s)| ≤ (τ³|s|³/3) e^{−1.79|s|} + 2τ e^{−√2/τ},  g₁ = log(1+e^{−2|s|}).
On |s| ≤ π/(4τ), e^{−√2/τ} ≤ e^{−(√2−π/4)/τ} e^{−|s|}. So this is (H2) with N = 3, g₂ = 0, c = 1.

**Order N = 4.** Using x³/6 − x⁵/120 ≤ x − sin x ≤ x³/6 and second-order Taylor of F (ξ ≥ 1.79|s|, so |F″(ξ)| ≤ e^{−1.79|s|}):
  L = τ g₁ + τ³ g₃ + O(τ⁵ (1+|s|)⁶ e^{−1.79|s|}) + O(τ e^{−√2/τ}),  g₃(s) = |s|³/(3(1+e^{2|s|})) = (|s|³/6)(1 − tanh|s|),
which agrees with AV5b's profile.

**Proposition 4 (convexity, PROVED for all τ ∈ (0,1]).** h is C², h + h″ > 0, so h is the support function of a convex body, and I(τ) is its perimeter.

*Proof.* On [0, π/2] the maximizing vertex is (1,1), so there h = Φ, where
  Φ(θ) = M f,  M = cos θ + sin θ,  f = (1+ε)^τ,  ε = Σ_{i=2}^4 E_i,  E_i = e^{ψ_i/τ},
  ψ₂ = −2 sin θ,  ψ₃ = −2 cos θ,  ψ₄ = −2M.

*C² at the kinks.* h is even about 0 by the symmetry θ ↦ −θ of the vertex set, so near 0, h(θ) = Φ(|θ|). At 0: M = M′ = 1, 1 + ε = S := 2 + 2e^{−2/τ}, and ε′ = −S/τ. Hence
  Φ′(0) = S^τ + τ S^{τ−1}(−S/τ) = 0,
so h′ and h″ = Φ″(|θ|) are continuous at 0. The same holds at every jπ/2 by D₄ symmetry.

*Sign of h + h″.* Because M″ = −M, Φ + Φ″ = 2M′f′ + Mf″. With ℓ = log(1+ε) this is
  Φ + Φ″ = τ f Q,  Q = 2M′ℓ′ + Mℓ″ + τMℓ′².
The functions ψ_i are linear in n, so ψ_i″ = −ψ_i. Hence ε′ = Σ(ψ_i′/τ)E_i and ε″ = Σ(ψ_i′²/τ² − ψ_i/τ)E_i, and
  τ²(1+ε)² Q = (1+ε) Σ_i E_i [τ(2M′ψ_i′ − Mψ_i) + Mψ_i′²] − (1−τ) M (Σ_i ψ_i′ E_i)².
Cauchy–Schwarz gives (Σψ_i′E_i)² ≤ ε Σψ_i′²E_i. Since (1−τ)M ≥ 0,
  τ²(1+ε)² Q ≥ Σ_i E_i [τ(1+ε)(2M′ψ_i′ − Mψ_i) + Mψ_i′²(1+τε)].

On [0, π/4] we have M ≥ 1, M′ = cos − sin ≥ 0, and ψ_i ≤ 0. Check each term:
- **i = 2** (ψ₂′ = −2cos): 2M′ψ₂′ − Mψ₂ = −4cos² + 6 sin cos + 2 sin². The term is ≥ 4cos²[M(1+τε) − τ(1+ε)] ≥ 4cos²(1−τ) ≥ 0.
- **i = 3** (ψ₃′ = 2 sin): 2M′ψ₃′ = 4(cos − sin)sin ≥ 0 and −Mψ₃ = 2M cos > 0. The term is > 0.
- **i = 4** (ψ₄′ = −2M′): 2M′ψ₄′ − Mψ₄ = 2M² − 4M′². The term is ≥ 4M′²(1−τ) + 2τ(1+ε)M² ≥ 0.

So Q > 0 on [0, π/4]. By the symmetry θ ↦ π/2 − θ the same holds on [π/4, π/2], and by D₄ on all of T. Then h + h″ = τ f Q > 0.

A C² periodic h with h + h″ ≥ 0 is a support function: its 1-homogeneous extension has Hessian eigenvalues 0 and (h+h″)/|x| off the origin, so it is convex along every segment that avoids 0, and by continuity along all segments; a convex 1-homogeneous function is sublinear (Schneider, *Convex Bodies*, §1.7). The boundary curve x(θ) = h n + h′ n′ has |x′| = h + h″, so per = ∫(h + h″) = ∫h. ∎

AV5b's numerical value min(h + h″) ≈ −2·10⁻⁵⁶ at dps 55 was rounding noise. With dps 400 the minimum is positive, for example 6.1·10⁻¹¹ at τ = 0.05 (`hypcheck_output.txt`), consistent with Q > 0.

**Status: (H0)–(H3) PROVED for N = 3 and N = 4; convexity PROVED.**

---

## 5. Profile moments, from one identity

Let g(s) = log(1+e^{−2s}) for s ≥ 0, and define γ_{ab} := (1/(a! b!)) ∫₀^∞ s^a g^b ds.

**Lemma 5 (reflection).** Let G(λ, μ) := ∫₀^∞ e^{λs}[(1+e^{−2s})^μ − 1] ds, which is analytic for |λ| < 2 and has Taylor coefficients γ_{ab}. Then, near (0,0),
  G(λ, μ) + G(2μ − λ, μ) = ½ [ Γ(−λ/2) Γ(λ/2 − μ)/Γ(−μ) + 2/λ + 2/(2μ − λ) ].

*Proof.* Substitute u = e^{−2s}: G = ½ K(−λ/2, μ), with K(a, μ) = ∫₀¹ u^{a−1}((1+u)^μ − 1) du, which is analytic for Re a > −1. For Re a > 0 and Re(a+μ) < 0, split ∫₀^∞ u^{a−1}(1+u)^μ du = B(a, −a−μ) at u = 1 and map (1, ∞) to (0, 1) by u ↦ 1/u. This gives
  K(a,μ) + 1/a + K(−a−μ, μ) − 1/(a+μ) = Γ(a)Γ(−a−μ)/Γ(−μ).
Both sides are analytic on the connected domain {|a| < 1, |a+μ| < 1}; the poles on the right cancel. The identity theorem extends the equation to that whole domain. Finally set a = −λ/2. ∎

**Coefficients.** At λ = 0 the bracket equals −ψ(1−μ) − γ_E = Σ_{k≥1} ζ(k+1) μ^k. On the left, the μ² and μ³ coefficients at λ = 0 are 2(γ₀₂ + γ₁₁) and 2γ₀₃ + 2γ₁₂ + 4γ₂₁. Hence:
- **γ₀₂ + γ₁₁ = ζ(3)/4**, that is, ∫₀^∞ (g²/2 + s g) ds = ζ(3)/4;
- 2γ₀₃ + 2γ₁₂ + 4γ₂₁ = ζ(4)/2.

The single moments come from log(1+x) = Σ(−1)^{n+1}xⁿ/n; the alternating partial sums are bounded by x, which justifies the interchange. This gives ∫₀^∞ s^j g = j! η(j+2)/2^{j+1}, so γ₀₁ = ζ(2)/4 and γ₂₁ = η(4)/8 = 7ζ(4)/64. Therefore
  **J = γ₁₂ + γ₀₃ = ζ(4)/4 − 2γ₂₁ = ζ(4)/32.**
AV5b had J only numerically; it is now PROVED. Similarly ∫₀^∞ s³/(1+e^{2s}) ds = (3/8)η(4) = 21ζ(4)/64.

By the same identity at λ = μ, the frozen-tent layer has the closed form G(τ,τ) = ¼[Γ(−τ/2)²/Γ(−τ) + 4/τ]. Every coefficient of this model is a ζ-value. The diagonal is the fixed line of λ ↦ 2μ − λ.

**The coefficients.** Both examples have four kinks, A ≡ 1, S = 0, β = 1, g₂ = 0, and an even profile, so each kink contributes 2H times a half-line integral.

| | squircle (H = 1/√2) | soft-max (H = 1) |
|---|---|---|
| c₀ | 4√2 | 8 |
| c₁ | 0 | 0 |
| c₂ = 4H∫_ℝ g₁ = 8Hγ₀₁ | √2 ζ(2) | π²/3 |
| c₃ = 4H∫_ℝ(g₁²/2 + \|s\|g₁) = 8H(γ₀₂+γ₁₁) | **√2 ζ(3)** | **2ζ(3)** |
| c₄ = 8H[∫₀^∞ g₃ + J − γ₂₁] | (4√2)[−7/16 + 1/32 − 7/64]ζ(4) = **−(33√2/16) ζ(4)** | 8[7/64 + 1/32 − 7/64]ζ(4) = **ζ(4)/4** |

All entries are PROVED via Theorem 1 at N = 4, giving P = Σ_{n≤4} c_n τⁿ + O(τ⁵). The squircle c₄ = −(33√2/16)ζ(4) is a new result. The soft-max c₄ = ζ(4)/4 was numerical in the paper.

---

## 6. Numerical sanity checks (mpmath, dps 50; own code)

P(τ) = 8∫₀^{π/4}(…) dθ by tanh-sinh quadrature, with breakpoints doubling in units of τ toward the layer. The reported quadrature error is ≤ 2·10⁻⁵¹ throughout. Source: `lay1_numerics.py`; output: `numerics_output.txt`.

**Squircle** (predicted c₄ = −3.156937154220705612…)

| 1/τ | (P − c₀ − c₂τ² − c₃τ³)/τ⁴ | (P − … − c₄τ⁴)/τ⁵ |
|---|---|---|
| 10 | −3.517957136452290404 | −3.610 |
| 40 | −3.292837787595325946 | −5.436 |
| 160 | −3.193564047823720326 | −5.860 |
| 640 | −3.166250557978090184 | −5.961 |
| 2560 | −3.159275134043971507 | −5.985 |

**Soft-max** (predicted c₄ = ζ(4)/4 = 0.270580808427784548…)

| 1/τ | (P − c₀ − c₂τ² − c₃τ³)/τ⁴ | (P − … − c₄τ⁴)/τ⁵ |
|---|---|---|
| 10 | 0.426193483767337137 | 1.5561 |
| 40 | 0.307943501167550406 | 1.4945 |
| 160 | 0.279895868125351200 | 1.4904 |
| 640 | 0.272908939020196591 | 1.49000 |
| 2560 | 0.271162816420343349 | 1.48994 |

In both cases R₃/τ⁴ stays bounded and converges to the proved c₄, and R₄/τ⁵ stays bounded. The limits of R₄/τ⁵ (≈ −5.99 and ≈ 1.4899) are c₅, which I did not identify.

Other checks:
- `lay1_moments.py` → `moments_output.txt`:
  - Lemma 5 residual at (λ, μ) = (0.3, −0.45): 1.9·10⁻⁵⁷.
  - Moments and J: residuals ≤ 10⁻⁶¹.
  - Diagonal closed form: residual ≤ 6·10⁻⁵¹.
- `lay1_hypcheck.py` → `hypcheck_output.txt`, a grid check at dps 400 of the explicit window bounds in §3 and §4. The largest ratio of actual to claimed bound is 0.707 (squircle) and 0.715 (soft-max); the bound requires ≤ 1. At dps 40 one spurious ratio of 9.9 appeared at s ≈ 47, τ = 10⁻³, where e^{−2s} is below working precision; this is noted so nobody is misled by it.

---

## 7. Status

| item | status |
|---|---|
| Theorem 1 (Layer Theorem, any order N, remainder O(τ^{N+1})) | **PROVED** |
| Gluing / recombination of s-independent parts | **PROVED**: none needed (Step 0, exact splitting against m·A) |
| Convexity needed? | Only to read ∫h as a perimeter (Cauchy). Theorem 1 itself is about ∫h. |
| Squircle: exact h_p, (H0)–(H3) for N = 3 and 4 | **PROVED** |
| Squircle: c₂ = √2ζ(2), c₃ = √2ζ(3), remainder O(τ⁴) | **PROVED** |
| Squircle: c₄ = −(33√2/16)ζ(4), remainder O(τ⁵) | **PROVED** (new) |
| Soft-max: (H0)–(H3) for N = 3 and 4 | **PROVED** |
| Soft-max: convexity, all τ ∈ (0,1] | **PROVED** (was numerical) |
| Soft-max: c₂ = π²/3, c₃ = 2ζ(3), c₄ = ζ(4)/4, remainder O(τ⁵) | **PROVED** (c₄ was numerical) |
| J = ∫₀^∞(sg²/2 + g³/6) = ζ(4)/32 | **PROVED** (Lemma 5; was numerical) |
| Remainders bounded as τ → 0 (both examples, 1/τ up to 2560) | **NUMERICAL**, consistent |
| c₅ of either example | **OPEN** (numerically ≈ −5.99 and ≈ 1.4899) |
| Hexagon / m ≥ 3 cases of AV5b under these hypotheses | **NOT EXAMINED** (see the remark on nonzero limits of L) |

No step failed, so no extra assumption was needed. The only hypotheses beyond the formal lemma are the decay and remainder bounds in (H2), the smoothness of a_k on the windows in (H1), and the exponential bound (H3). All of them were verified for both examples.

## 8. Summary for a referee

Write ∫h = ∫m·A + Σ_j ∫_{window j} m·A·(e^{L_j} − 1) + ∫_outside (h − m·A). The first term is a uniform τ-expansion of the bulk. In the second, the integrand after θ = θ_j + τs is a product of three factors (tent, bulk, layer), each expanded with an explicit Taylor or Cauchy-majorant remainder (Lemma 1); this gives an error C τ^N (1+|s|)^{2N} e^{−c|s|}, which integrates to O(τ^{N+1}). The tails beyond W/τ and the outside region are O(e^{−const/τ}). The hypotheses hold for the squircle (h_p = ‖n‖_q, log r = −2 artanh tan|u|, Λ′ = sec 2x) and for the soft-max (exact four-term window formula, x − sin x ≤ x³/6, and a Cauchy–Schwarz proof that h + h″ > 0). The resulting moments all follow from the Beta reflection identity G(λ,μ) + G(2μ−λ,μ) = ½[Γ(−λ/2)Γ(λ/2−μ)/Γ(−μ) + 2/λ + 2/(2μ−λ)], which gives c₃ = √2ζ(3) and 2ζ(3), and c₄ = −(33√2/16)ζ(4) and ζ(4)/4, all confirmed at dps 50.

---

## 9. Planner notes after the referee (2026-11-22)

Refereed by an independent Opus agent with its own code (`REFEREE.md`, `referee/`): all five items hold. Fixes to read into the text above:
1. (H0): read "W at most half the least angular gap, with the tent formula on all of I_j". The windows may then share endpoints (a null set), which Step 0 allows. This is what the soft-max uses (W = π/4).
2. Lemma 5, pole cancellation: Γ(a)Γ(−a−μ)/Γ(−μ) − 1/a + 1/(a+μ) = μ/(a(a+μ))·[Γ(1+a)Γ(1−a−μ)/Γ(1−μ) − 1]; the bracket vanishes on a = 0 and on a + μ = 0, so the right side is analytic on the whole domain.
3. Section 5, "Coefficients": the μ^k coefficient on the right is ζ(k+1)/2 (the factor ½ of the identity).
4. Define Φ_{j,−1} := 0; in the (H2) remark, c becomes min(c, c′).
5. Not new: the squircle's c₄ = −(33√2/16)ζ(4) is the τ⁴ term of the paper's pure-weight theorem (thm:pw), and that theorem's τ⁵ term −√2(¾ζ(5) + (7/4)ζ(2)ζ(3)) = −5.99341 is the "open" c₅ ≈ −5.99 of Section 6. Both are cross-checks.
6. The soft-max c₅ is in `../SM5_softmax_c5/result.md`: (107/32)ζ(5) − ζ(2)ζ(3), 55 digits against quadrature.
